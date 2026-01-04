"""
Compute WEIGHTED RWR influence matrix M using STRING confidence scores.

Enhancement based on Köhler et al. 2008:
    "Nodes connected by multiple paths receive higher similarity than nodes 
    connected by only one path."

This script:
1. Loads the STRING network WITH confidence scores (900-999)
2. Builds a weighted transition matrix (edges weighted by confidence)
3. Computes M = α(I - (1-α)W)^{-1} with weighted W
4. Saves alongside the original unweighted matrix for comparison

Output: research/ecnp-closed-form/data/influence_matrix_900_weighted.npz

Author: Research Session 2026-01-02
Status: Implementation of Algorithm Improvement Area 1 (Edge Weights)
"""

import numpy as np
import pandas as pd
import networkx as nx
from pathlib import Path
import time

# Paths - script is in ecnp-closed-form/src/precompute/, so 4 levels up to project root
PROJECT_ROOT = Path(__file__).resolve().parents[4]
DATA_DIR = PROJECT_ROOT / "data" / "processed"
OUTPUT_DIR = PROJECT_ROOT / "research" / "ecnp-closed-form" / "data"

# RWR parameters
ALPHA = 0.15  # Restart probability


def load_weighted_network():
    """
    Load STRING network with confidence scores, filter to liver LCC.
    
    The processed network_900.parquet has weights, but the liver_lcc version doesn't.
    We need to:
    1. Load the weighted network
    2. Load liver genes
    3. Filter and extract LCC while preserving weights
    """
    # Load weighted network
    network_path = DATA_DIR / "network_900.parquet"
    edges = pd.read_parquet(network_path)
    print(f"Loaded weighted network: {len(edges)} edges")
    print(f"  Weight range: {edges['weight'].min()} - {edges['weight'].max()}")
    
    # Standardize column names
    edges = edges.rename(columns={'protein1': 'gene1', 'protein2': 'gene2'})
    
    # Load liver genes
    liver_df = pd.read_csv(DATA_DIR / "liver_proteome.csv")
    liver_genes = set(liver_df['gene_symbol'])
    print(f"Liver genes: {len(liver_genes)}")
    
    # Filter to liver-expressed genes
    edges_liver = edges[
        edges['gene1'].isin(liver_genes) & 
        edges['gene2'].isin(liver_genes)
    ].copy()
    print(f"After liver filter: {len(edges_liver)} edges")
    
    # Build graph and extract LCC
    G = nx.from_pandas_edgelist(edges_liver, 'gene1', 'gene2', edge_attr='weight')
    lcc_nodes = max(nx.connected_components(G), key=len)
    G_lcc = G.subgraph(lcc_nodes).copy()
    print(f"LCC: {G_lcc.number_of_nodes()} nodes, {G_lcc.number_of_edges()} edges")
    
    # Convert back to edgelist with weights
    edges_lcc = []
    for u, v, data in G_lcc.edges(data=True):
        edges_lcc.append({
            'gene1': u,
            'gene2': v,
            'weight': data.get('weight', 900)  # Default to minimum if missing
        })
    
    df_lcc = pd.DataFrame(edges_lcc)
    return df_lcc


def build_weighted_transition_matrix(edges, weight_transform='normalized'):
    """
    Build column-stochastic transition matrix W from weighted edge list.
    
    Weight transforms:
    - 'binary': Ignore weights (for comparison baseline)
    - 'normalized': w_ij = score / max_score (0.9 to 1.0)
    - 'log': w_ij = log(score / 900 + 1) (emphasizes differences)
    - 'linear': w_ij = (score - 900) / 100 + 1 (1.0 to 2.0 range)
    
    Args:
        edges: DataFrame with 'gene1', 'gene2', 'weight' columns
        weight_transform: How to transform STRING scores
    
    Returns:
        W: numpy array (n x n) transition matrix
        node_list: list of gene symbols
        degrees: degree array
    """
    # Get unique nodes
    all_nodes = pd.concat([edges['gene1'], edges['gene2']]).unique()
    node_list = sorted(all_nodes)
    node_to_idx = {g: i for i, g in enumerate(node_list)}
    n = len(node_list)
    
    print(f"Building weighted transition matrix for {n} nodes...")
    print(f"  Transform: {weight_transform}")
    
    # Build weighted adjacency matrix
    A = np.zeros((n, n))
    
    for _, row in edges.iterrows():
        i = node_to_idx[row['gene1']]
        j = node_to_idx[row['gene2']]
        score = row['weight']
        
        # Apply weight transform
        if weight_transform == 'binary':
            w = 1.0
        elif weight_transform == 'normalized':
            w = score / 1000.0  # 0.9 to 0.999
        elif weight_transform == 'log':
            w = np.log(score / 900 + 1)  # ~0 to ~0.1
        elif weight_transform == 'linear':
            w = (score - 900) / 100 + 1.0  # 1.0 to 2.0
        else:
            raise ValueError(f"Unknown weight transform: {weight_transform}")
        
        A[i, j] = w
        A[j, i] = w  # Symmetric
    
    # Compute weighted degrees
    degrees = A.sum(axis=1)
    
    # Column-normalize to get transition matrix
    col_sums = A.sum(axis=0)
    col_sums[col_sums == 0] = 1  # Avoid division by zero
    W = A / col_sums[np.newaxis, :]
    
    print(f"  Weight stats: mean={A[A>0].mean():.4f}, std={A[A>0].std():.4f}")
    print(f"  Degree stats: mean={degrees.mean():.1f}, max={degrees.max():.0f}")
    
    return W, node_list, degrees


def compute_influence_matrix(W, alpha=ALPHA):
    """
    Compute M = α(I - (1-α)W)^{-1}
    
    For n ≈ 8000, dense inversion takes ~30-60 seconds.
    """
    n = W.shape[0]
    I = np.eye(n)
    
    # (I - (1-α)W)
    A = I - (1 - alpha) * W
    
    print(f"Computing influence matrix M for n={n}...")
    print("This may take 30-60 seconds...")
    
    start = time.time()
    A_inv = np.linalg.inv(A)
    M = alpha * A_inv
    elapsed = time.time() - start
    
    print(f"Computed M in {elapsed:.1f}s")
    print(f"M shape: {M.shape}")
    print(f"M stats: min={M.min():.6f}, max={M.max():.6f}, mean={M.mean():.6f}")
    
    return M


def apply_hub_correction(M, degrees, beta=0.5):
    """
    Apply hub correction: down-weight influence from high-degree nodes.
    
    Köhler et al. 2008 insight:
        "Proteins x and y connected via a hub... global similarity is LESS 
        than via a low-degree intermediate"
    
    Implementation: Scale columns by inverse degree factor.
    M_corrected[i,j] = M[i,j] * (1 / (1 + normalized_degree[j]))^beta
    
    Args:
        M: Original influence matrix
        degrees: Array of node degrees
        beta: Correction strength (0=none, 0.5=moderate, 1=strong)
    
    Returns:
        M_corrected: Hub-corrected influence matrix
    """
    print(f"\nApplying hub correction (beta={beta})...")
    
    # Normalize degrees to [0, 1]
    deg_normalized = degrees / degrees.max()
    
    # Inverse degree factor
    inv_deg_factor = (1.0 / (1 + deg_normalized)) ** beta
    
    # Apply column-wise (influence RECEIVED is penalized by source degree)
    M_corrected = M * inv_deg_factor.reshape(1, -1)
    
    # Re-normalize to preserve scale
    scale = M.mean() / M_corrected.mean()
    M_corrected *= scale
    
    print(f"  M_corrected stats: min={M_corrected.min():.6f}, max={M_corrected.max():.6f}")
    
    return M_corrected


def compute_dili_influence_vector(M, node_list, dili_genes_path):
    """
    Compute m_j = Σ_{d ∈ DILI} M[j, d] for each node j.
    
    This is the per-node influence on the DILI module.
    """
    # Load DILI genes
    dili_df = pd.read_csv(dili_genes_path)
    dili_genes = set(dili_df['gene_name'])  # Column is 'gene_name' not 'gene'
    
    node_to_idx = {g: i for i, g in enumerate(node_list)}
    dili_indices = [node_to_idx[g] for g in dili_genes if g in node_to_idx]
    
    print(f"\nComputing DILI influence vector...")
    print(f"  DILI genes in network: {len(dili_indices)}")
    
    # Sum influence to DILI nodes for each source node
    M_dili = M[:, dili_indices]
    m_vector = M_dili.sum(axis=1)
    
    print(f"  m_vector stats: min={m_vector.min():.6f}, max={m_vector.max():.6f}")
    
    # Create DataFrame
    m_df = pd.DataFrame({
        'gene': node_list,
        'dili_influence': m_vector
    })
    
    return m_df


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Load weighted network
    edges = load_weighted_network()
    
    # Save the weighted liver LCC network for reference
    weighted_lcc_path = DATA_DIR / "network_900_liver_lcc_weighted.parquet"
    edges.to_parquet(weighted_lcc_path, index=False)
    print(f"\nSaved weighted LCC network: {weighted_lcc_path}")
    
    # Build weighted transition matrix
    W, node_list, degrees = build_weighted_transition_matrix(edges, weight_transform='normalized')
    
    # Compute influence matrix
    M = compute_influence_matrix(W)
    
    # Apply hub correction (Algorithm Improvement Area 4)
    M_corrected = apply_hub_correction(M, degrees, beta=0.5)
    
    # Compute DILI influence vectors
    dili_path = DATA_DIR / "dili_900_lcc.csv"
    m_df = compute_dili_influence_vector(M, node_list, dili_path)
    m_df_corrected = compute_dili_influence_vector(M_corrected, node_list, dili_path)
    
    # Save weighted influence matrix
    output_path = OUTPUT_DIR / "influence_matrix_900_weighted.npz"
    np.savez_compressed(
        output_path,
        M=M,
        M_corrected=M_corrected,
        node_list=np.array(node_list),
        degrees=degrees
    )
    print(f"\nSaved weighted influence matrix: {output_path}")
    
    # Save DILI influence vectors
    m_df.to_csv(OUTPUT_DIR / "dili_influence_vector_900_weighted.csv", index=False)
    print(f"Saved: {OUTPUT_DIR / 'dili_influence_vector_900_weighted.csv'}")
    
    m_df_corrected.to_csv(OUTPUT_DIR / "dili_influence_vector_900_corrected.csv", index=False)
    print(f"Saved: {OUTPUT_DIR / 'dili_influence_vector_900_corrected.csv'}")
    
    # Print comparison stats
    print("\n" + "="*60)
    print("COMPARISON: Original vs Weighted Influence")
    print("="*60)
    
    # Load original for comparison
    orig = np.load(OUTPUT_DIR / "influence_matrix_900.npz", allow_pickle=True)
    M_orig = orig['M']
    
    print(f"\nOriginal (binary edges):")
    print(f"  M range: [{M_orig.min():.6f}, {M_orig.max():.6f}]")
    print(f"  M mean:  {M_orig.mean():.6f}")
    
    print(f"\nWeighted (STRING scores):")
    print(f"  M range: [{M.min():.6f}, {M.max():.6f}]")
    print(f"  M mean:  {M.mean():.6f}")
    
    print(f"\nWeighted + Hub Corrected:")
    print(f"  M range: [{M_corrected.min():.6f}, {M_corrected.max():.6f}]")
    print(f"  M mean:  {M_corrected.mean():.6f}")
    
    # Correlation between matrices
    corr_weighted = np.corrcoef(M_orig.flatten(), M.flatten())[0, 1]
    corr_corrected = np.corrcoef(M_orig.flatten(), M_corrected.flatten())[0, 1]
    print(f"\nCorrelation with original:")
    print(f"  Weighted:  r = {corr_weighted:.4f}")
    print(f"  Corrected: r = {corr_corrected:.4f}")


if __name__ == "__main__":
    main()

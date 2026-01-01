"""
Compute the RWR influence matrix M = α(I - (1-α)W)^{-1}

This is the fundamental matrix for closed-form ECNP derivation.
Computing M once allows O(1) influence lookups for any node pair.

NOTE FOR SCALABILITY:
For production use on very large networks (n > 20,000), M should NOT be
materialized explicitly. Instead, use sparse linear solves or truncated
Krylov methods. This script uses dense inversion for validation on the
liver LCC (~8,000 nodes), which is tractable.

Output: research/ecnp-closed-form/data/influence_matrix_900.npz
"""

import numpy as np
import pandas as pd
from pathlib import Path
import time

# Paths - script is at research/ecnp-closed-form/scripts/, so 3 levels up
PROJECT_ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = PROJECT_ROOT / "data" / "processed"
OUTPUT_DIR = PROJECT_ROOT / "research" / "ecnp-closed-form" / "data"

# RWR parameters
ALPHA = 0.15  # Restart probability


def load_network():
    """Load the liver-expressed LCC network (STRING ≥900)."""
    network_path = DATA_DIR / "network_900_liver_lcc.parquet"
    edges = pd.read_parquet(network_path)
    print(f"Loaded network: {len(edges)} edges")
    return edges


def build_transition_matrix(edges):
    """
    Build column-stochastic transition matrix W from edge list.
    
    Args:
        edges: DataFrame with 'gene1' and 'gene2' columns
    
    Returns:
        W: numpy array (n x n) transition matrix
        node_list: list of gene symbols
    """
    # Get unique nodes
    all_nodes = pd.concat([edges['gene1'], edges['gene2']]).unique()
    node_list = sorted(all_nodes)
    node_to_idx = {g: i for i, g in enumerate(node_list)}
    n = len(node_list)
    
    print(f"Building transition matrix for {n} nodes...")
    
    # Build adjacency matrix (dense for small n)
    A = np.zeros((n, n))
    for _, row in edges.iterrows():
        i = node_to_idx[row['gene1']]
        j = node_to_idx[row['gene2']]
        A[i, j] = 1
        A[j, i] = 1  # Symmetric (undirected)
    
    # Column-normalize to get transition matrix
    col_sums = A.sum(axis=0)
    col_sums[col_sums == 0] = 1  # Avoid division by zero
    W = A / col_sums[np.newaxis, :]
    
    print(f"Built transition matrix: {W.shape}")
    return W, node_list


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
    print(f"M shape: {M.shape}, density: {(np.abs(M) > 1e-10).sum() / M.size:.2%}")
    
    return M


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Load network
    edges = load_network()
    
    # Build transition matrix
    W, node_list = build_transition_matrix(edges)
    
    # Compute influence matrix
    M = compute_influence_matrix(W)
    
    # Save
    output_path = OUTPUT_DIR / "influence_matrix_900.npz"
    np.savez_compressed(
        output_path,
        M=M,
        node_list=np.array(node_list)
    )
    print(f"\nSaved: {output_path}")
    
    # Also save node list separately for convenience
    pd.DataFrame({'gene': node_list}).to_csv(
        OUTPUT_DIR / "node_list_900.csv", index=False
    )
    print(f"Saved: {OUTPUT_DIR / 'node_list_900.csv'}")


if __name__ == "__main__":
    main()

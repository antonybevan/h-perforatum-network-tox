"""
Compute per-node DILI-influence vector: m_j = Σ_{i∈D} M_{ij}

This vector allows O(1) lookup of how much influence any node has on the DILI module.
Key insight: I(T) = Σ_{j∈T} m_j (influence is additive due to RWR linearity).

Output: research/ecnp-closed-form/data/dili_influence_vector_900.csv
"""

import numpy as np
import pandas as pd
from pathlib import Path

# Paths - script is at research/ecnp-closed-form/scripts/, so 3 levels up
PROJECT_ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = PROJECT_ROOT / "data" / "processed"
RESEARCH_DATA_DIR = PROJECT_ROOT / "research" / "ecnp-closed-form" / "data"


def load_influence_matrix():
    """Load precomputed M matrix."""
    npz_path = RESEARCH_DATA_DIR / "influence_matrix_900.npz"
    data = np.load(npz_path, allow_pickle=True)
    M = data['M']
    node_list = data['node_list'].tolist()
    print(f"Loaded M: {M.shape}")
    return M, node_list


def load_dili_genes():
    """Load DILI gene set (in LCC)."""
    dili_path = DATA_DIR / "dili_900_lcc.csv"
    dili = pd.read_csv(dili_path)
    dili_genes = set(dili['gene_name'].tolist())
    print(f"Loaded {len(dili_genes)} DILI genes")
    return dili_genes


def compute_dili_influence_vector(M, node_list, dili_genes):
    """
    Compute m_j = Σ_{i∈D} M_{ij} for all nodes j.
    
    This is the sum of column j over rows in the DILI set.
    """
    node_to_idx = {g: i for i, g in enumerate(node_list)}
    
    # Get indices of DILI genes
    dili_indices = [node_to_idx[g] for g in dili_genes if g in node_to_idx]
    print(f"DILI genes in network: {len(dili_indices)}")
    
    # Sum over DILI rows for each column
    # m_j = Σ_{i∈D} M[i, j]
    m_vector = M[dili_indices, :].sum(axis=0)
    
    print(f"DILI-influence vector: min={m_vector.min():.6f}, max={m_vector.max():.6f}, mean={m_vector.mean():.6f}")
    
    return m_vector


def main():
    # Load data
    M, node_list = load_influence_matrix()
    dili_genes = load_dili_genes()
    
    # Compute per-node DILI-influence
    m_vector = compute_dili_influence_vector(M, node_list, dili_genes)
    
    # Save as DataFrame
    df = pd.DataFrame({
        'gene': node_list,
        'dili_influence': m_vector
    })
    
    output_path = RESEARCH_DATA_DIR / "dili_influence_vector_900.csv"
    df.to_csv(output_path, index=False)
    print(f"\nSaved: {output_path}")
    
    # Print top influencers
    print("\nTop 10 nodes by DILI-influence:")
    print(df.nlargest(10, 'dili_influence').to_string(index=False))


if __name__ == "__main__":
    main()

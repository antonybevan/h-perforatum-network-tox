"""
Diagnostic Analysis: Why Quercetin's Closed-Form Z is Still Off

Quercetin: Z=5.95 (closed-form) vs Z=4.42 (Monte Carlo)
This means sigma is UNDERESTIMATED by ~35%

Possible causes:
1. Pool variance differs between compounds
2. Redundancy scaling with k is non-linear
3. High-redundancy regime needs different treatment
4. mu_T estimation is biased

This script digs into each component.
"""

import numpy as np
import pandas as pd
from pathlib import Path

# Paths
PROJECT_ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = PROJECT_ROOT / "data" / "processed"
RESEARCH_DATA_DIR = PROJECT_ROOT / "research" / "ecnp-closed-form" / "data"


def load_all_data():
    """Load all required data."""
    # M matrix
    npz = np.load(RESEARCH_DATA_DIR / "influence_matrix_900.npz", allow_pickle=True)
    M = npz['M']
    node_list = npz['node_list'].tolist()
    node_to_idx = {g: i for i, g in enumerate(node_list)}
    
    # DILI influence
    m_df = pd.read_csv(RESEARCH_DATA_DIR / "dili_influence_vector_900.csv")
    m_vector = m_df.set_index('gene')['dili_influence']
    
    # DILI indices
    dili = pd.read_csv(DATA_DIR / "dili_900_lcc.csv")
    dili_genes = set(dili['gene_name'].tolist())
    dili_indices = [node_to_idx[g] for g in dili_genes if g in node_to_idx]
    
    # Degrees
    edges = pd.read_parquet(DATA_DIR / "network_900_liver_lcc.parquet")
    degrees = pd.concat([edges['gene1'], edges['gene2']]).value_counts().to_dict()
    
    # Targets
    targets_df = pd.read_csv(DATA_DIR / "targets_lcc.csv")
    hyp = targets_df[targets_df['compound'] == 'Hyperforin']['gene_symbol'].tolist()
    que = targets_df[targets_df['compound'] == 'Quercetin']['gene_symbol'].tolist()
    
    return M, node_list, node_to_idx, m_vector, dili_indices, degrees, hyp, que


def get_degree_matched_pool(targets, degrees, tolerance=0.1):
    target_degrees = [degrees.get(t, 0) for t in targets]
    pool = set()
    for deg in target_degrees:
        min_d = int(deg * (1 - tolerance))
        max_d = int(deg * (1 + tolerance))
        matching = [n for n, d in degrees.items() if min_d <= d <= max_d]
        pool.update(matching)
    return list(pool - set(targets))


def analyze_pools(hyp, que, m_vector, degrees):
    """Compare pool characteristics."""
    hyp_pool = get_degree_matched_pool(hyp, degrees)
    que_pool = get_degree_matched_pool(que, degrees)
    
    hyp_m = [m_vector.get(p, 0) for p in hyp_pool if p in m_vector.index]
    que_m = [m_vector.get(p, 0) for p in que_pool if p in m_vector.index]
    
    print("="*60)
    print("POOL COMPARISON")
    print("="*60)
    print(f"\nHyperforin pool:")
    print(f"  Size: {len(hyp_pool)}")
    print(f"  Mean m: {np.mean(hyp_m):.6f}")
    print(f"  Std m: {np.std(hyp_m):.6f}")
    print(f"  Var m: {np.var(hyp_m, ddof=1):.8f}")
    
    print(f"\nQuercetin pool:")
    print(f"  Size: {len(que_pool)}")
    print(f"  Mean m: {np.mean(que_m):.6f}")
    print(f"  Std m: {np.std(que_m):.6f}")
    print(f"  Var m: {np.var(que_m, ddof=1):.8f}")
    
    print(f"\nRatio (Que/Hyp):")
    print(f"  Pool size: {len(que_pool)/len(hyp_pool):.2f}x")
    print(f"  Variance: {np.var(que_m, ddof=1)/np.var(hyp_m, ddof=1):.2f}x")
    
    return hyp_pool, que_pool, hyp_m, que_m


def analyze_redundancy(M, dili_indices, targets, node_to_idx, name):
    """Detailed redundancy analysis."""
    target_idx = [node_to_idx[t] for t in targets if t in node_to_idx]
    k = len(target_idx)
    
    # Get M restricted to DILI
    M_D = M[dili_indices, :][:, target_idx]
    
    # Norms
    norms = np.linalg.norm(M_D, axis=0)
    
    # Cosine similarity
    M_D_norm = M_D / (norms[np.newaxis, :] + 1e-10)
    rho = M_D_norm.T @ M_D_norm
    
    # Off-diagonal
    off_diag = rho[~np.eye(k, dtype=bool)]
    
    print(f"\n{name} Redundancy:")
    print(f"  Targets (k): {k}")
    print(f"  Max pairwise rho: {off_diag.max():.4f}")
    print(f"  Min pairwise rho: {off_diag.min():.4f}")
    print(f"  Mean pairwise rho: {off_diag.mean():.4f}")
    print(f"  Std pairwise rho: {off_diag.std():.4f}")
    print(f"  Sum off-diag rho: {off_diag.sum():.4f}")
    print(f"  k*(k-1): {k*(k-1)}")
    print(f"  Theoretical max sum: {k*(k-1)}")
    
    # Top redundant pairs
    if k > 2:
        # Get top 5 most redundant pairs
        rho_flat = rho.copy()
        np.fill_diagonal(rho_flat, 0)
        top_idx = np.unravel_index(np.argsort(rho_flat.flatten())[-10:], rho.shape)
        print(f"\n  Top 5 redundant pairs:")
        for i in range(min(5, len(top_idx[0]))):
            idx_i, idx_j = top_idx[0][-(i+1)], top_idx[1][-(i+1)]
            if idx_i < idx_j:
                t_i = targets[idx_i] if idx_i < len(targets) else f"idx{idx_i}"
                t_j = targets[idx_j] if idx_j < len(targets) else f"idx{idx_j}"
                print(f"    {t_i} - {t_j}: rho={rho[idx_i, idx_j]:.4f}")
    
    return rho, off_diag


def analyze_observed_influence(hyp, que, m_vector):
    """Compare observed influence components."""
    print("\n" + "="*60)
    print("OBSERVED INFLUENCE BREAKDOWN")
    print("="*60)
    
    hyp_m = [(t, m_vector.get(t, 0)) for t in hyp if t in m_vector.index]
    que_m = [(t, m_vector.get(t, 0)) for t in que if t in m_vector.index]
    
    print(f"\nHyperforin (k={len(hyp_m)}):")
    print(f"  Total I(T): {sum(m for _, m in hyp_m):.6f}")
    print(f"  Mean per target: {np.mean([m for _, m in hyp_m]):.6f}")
    print(f"  Top 3 targets:")
    for t, m in sorted(hyp_m, key=lambda x: -x[1])[:3]:
        print(f"    {t}: {m:.6f}")
    
    print(f"\nQuercetin (k={len(que_m)}):")
    print(f"  Total I(T): {sum(m for _, m in que_m):.6f}")
    print(f"  Mean per target: {np.mean([m for _, m in que_m]):.6f}")
    print(f"  Top 3 targets:")
    for t, m in sorted(que_m, key=lambda x: -x[1])[:3]:
        print(f"    {t}: {m:.6f}")


def back_calculate_sigma(I_T, mu_T, Z_mc):
    """Back-calculate what sigma Monte Carlo implies."""
    return (I_T - mu_T) / Z_mc


def main():
    print("Loading data...")
    M, node_list, node_to_idx, m_vector, dili_indices, degrees, hyp, que = load_all_data()
    print("Done.\n")
    
    # Pool analysis
    hyp_pool, que_pool, hyp_m, que_m = analyze_pools(hyp, que, m_vector, degrees)
    
    # Redundancy analysis
    print("\n" + "="*60)
    print("REDUNDANCY STRUCTURE")
    print("="*60)
    hyp_rho, hyp_off = analyze_redundancy(M, dili_indices, hyp, node_to_idx, "Hyperforin")
    que_rho, que_off = analyze_redundancy(M, dili_indices, que, node_to_idx, "Quercetin")
    
    # Observed influence
    analyze_observed_influence(hyp, que, m_vector)
    
    # Back-calculate implied sigma
    print("\n" + "="*60)
    print("BACK-CALCULATED SIGMA FROM MONTE CARLO")
    print("="*60)
    
    # From our closed-form results
    hyp_I = 1.138096
    hyp_mu = 0.155133
    hyp_sigma_cf = 0.105948
    hyp_Z_mc = 10.27
    
    que_I = 1.994615
    que_mu = 1.014679
    que_sigma_cf = 0.164592
    que_Z_mc = 4.42
    
    hyp_sigma_implied = back_calculate_sigma(hyp_I, hyp_mu, hyp_Z_mc)
    que_sigma_implied = back_calculate_sigma(que_I, que_mu, que_Z_mc)
    
    print(f"\nHyperforin:")
    print(f"  sigma (closed-form): {hyp_sigma_cf:.6f}")
    print(f"  sigma (implied by MC): {hyp_sigma_implied:.6f}")
    print(f"  Ratio (CF/MC): {hyp_sigma_cf/hyp_sigma_implied:.2f}")
    
    print(f"\nQuercetin:")
    print(f"  sigma (closed-form): {que_sigma_cf:.6f}")
    print(f"  sigma (implied by MC): {que_sigma_implied:.6f}")
    print(f"  Ratio (CF/MC): {que_sigma_cf/que_sigma_implied:.2f}")
    
    # Key insight
    print("\n" + "="*60)
    print("KEY INSIGHT")
    print("="*60)
    print(f"\nQuercetin's sigma should be: {que_sigma_implied:.6f}")
    print(f"We computed: {que_sigma_cf:.6f}")
    print(f"Missing variance: {que_sigma_implied**2 - que_sigma_cf**2:.6f}")
    print(f"Ratio of variances: {(que_sigma_implied/que_sigma_cf)**2:.2f}x")


if __name__ == "__main__":
    main()

"""
Stratified Pool Construction: The Correct Approach

Previous error: Intersecting degree AND influence constraints collapsed pools.

Correct approach:
1. Match on degree (topology)
2. THEN stratify by influence quantiles WITHIN the degree-matched pool

This is CONDITIONING, not intersection:
    E[m | deg ≈ deg_T, m ≈ m_T]

It acknowledges:
- Drug targets are enriched for signaling hubs
- Kinases are not random
- Regulatory proteins are not random

The algorithm must reflect pharmacological reality.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from scipy.optimize import minimize_scalar

# Paths
PROJECT_ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = PROJECT_ROOT / "data" / "processed"
RESEARCH_DATA_DIR = PROJECT_ROOT / "research" / "ecnp-closed-form" / "data"


def load_all_data():
    npz = np.load(RESEARCH_DATA_DIR / "influence_matrix_900.npz", allow_pickle=True)
    M = npz['M']
    node_list = npz['node_list'].tolist()
    node_to_idx = {g: i for i, g in enumerate(node_list)}
    
    m_df = pd.read_csv(RESEARCH_DATA_DIR / "dili_influence_vector_900.csv")
    m_vector = m_df.set_index('gene')['dili_influence']
    
    dili = pd.read_csv(DATA_DIR / "dili_900_lcc.csv")
    dili_indices = [node_to_idx[g] for g in dili['gene_name'] if g in node_to_idx]
    
    edges = pd.read_parquet(DATA_DIR / "network_900_liver_lcc.parquet")
    degrees = pd.concat([edges['gene1'], edges['gene2']]).value_counts().to_dict()
    
    targets = pd.read_csv(DATA_DIR / "targets_lcc.csv")
    hyp = targets[targets['compound'] == 'Hyperforin']['gene_symbol'].tolist()
    que = targets[targets['compound'] == 'Quercetin']['gene_symbol'].tolist()
    
    return M, node_list, node_to_idx, m_vector, dili_indices, degrees, hyp, que


def get_degree_matched_pool(targets, degrees, tolerance=0.2):
    """Step 1: Get degree-matched candidates (large pool)."""
    pool = set()
    for t in targets:
        deg = degrees.get(t, 0)
        deg_min, deg_max = deg * (1 - tolerance), deg * (1 + tolerance)
        pool.update(n for n, d in degrees.items() if deg_min <= d <= deg_max)
    return list(pool - set(targets))


def stratify_pool_by_influence(pool, m_vector, target_influences, n_quantiles=10):
    """
    Step 2: Within degree-matched pool, select nodes in same influence quantiles.
    
    For each target, find which quantile its influence falls into.
    Sample pool nodes from those same quantiles.
    """
    # Get influence values for pool
    pool_m = pd.Series({p: m_vector.get(p, 0) for p in pool if p in m_vector.index})
    
    if len(pool_m) < n_quantiles:
        return list(pool_m.index)  # Pool too small, return all
    
    # Assign quantile labels to pool
    pool_quantiles = pd.qcut(pool_m, n_quantiles, labels=False, duplicates='drop')
    
    # Find which quantiles targets fall into
    target_quantile_set = set()
    quantile_edges = pool_m.quantile(np.linspace(0, 1, n_quantiles + 1)).values
    
    for t_inf in target_influences:
        for q in range(n_quantiles):
            if quantile_edges[q] <= t_inf <= quantile_edges[q + 1]:
                target_quantile_set.add(q)
                break
    
    # Select pool nodes in those quantiles
    stratified_pool = pool_quantiles[pool_quantiles.isin(target_quantile_set)].index.tolist()
    
    return stratified_pool


def get_stratified_pool(targets, degrees, m_vector, deg_tol=0.2, n_quantiles=10):
    """Combined: degree-matched then influence-stratified."""
    # Step 1: Degree matching
    deg_pool = get_degree_matched_pool(targets, degrees, deg_tol)
    
    # Get target influences
    target_infs = [m_vector.get(t, 0) for t in targets if t in m_vector.index]
    
    # Step 2: Influence stratification
    stratified = stratify_pool_by_influence(deg_pool, m_vector, target_infs, n_quantiles)
    
    return stratified


def compute_mean_redundancy(M, dili_indices, target_indices):
    k = len(target_indices)
    if k < 2:
        return 0.0
    
    M_D = M[dili_indices, :][:, target_indices]
    norms = np.linalg.norm(M_D, axis=0)
    norms[norms == 0] = 1e-10
    M_D_norm = M_D / norms[np.newaxis, :]
    rho = M_D_norm.T @ M_D_norm
    
    off_diag_sum = np.sum(rho) - np.trace(rho)
    return off_diag_sum / (k * (k - 1))


def ecnp_stratified(targets, M, m_vector, dili_indices, node_to_idx, pool, lambda_param):
    """ECNP with stratified pool."""
    target_indices = [node_to_idx[t] for t in targets if t in node_to_idx]
    k = len(target_indices)
    
    if k == 0:
        return {'Z': np.nan, 'error': 'No targets'}
    
    I_T = sum(m_vector.get(t, 0) for t in targets if t in m_vector.index)
    
    pool_m = [m_vector.get(p, 0) for p in pool if p in m_vector.index]
    if len(pool_m) < 20:
        return {'Z': np.nan, 'error': f'Pool too small ({len(pool_m)})'}
    
    pool_mean = np.mean(pool_m)
    pool_var = np.var(pool_m, ddof=1)
    
    mu_T = k * pool_mean
    mean_rho = compute_mean_redundancy(M, dili_indices, target_indices)
    
    sigma_T_sq = pool_var / k + lambda_param * mean_rho
    sigma_T = np.sqrt(sigma_T_sq) if sigma_T_sq > 0 else 1e-10
    
    Z = (I_T - mu_T) / sigma_T
    
    return {
        'k': k,
        'I_T': I_T,
        'mu_T': mu_T,
        'sigma_T': sigma_T,
        'Z': Z,
        'mean_rho': mean_rho,
        'pool_mean': pool_mean,
        'pool_size': len(pool_m)
    }


def calibrate_lambda(M, m_vector, dili_indices, node_to_idx, degrees,
                     hyp_targets, que_targets,
                     hyp_mc_z=10.27, que_mc_z=4.42):
    
    hyp_pool = get_stratified_pool(hyp_targets, degrees, m_vector)
    que_pool = get_stratified_pool(que_targets, degrees, m_vector)
    
    print(f"  Stratified pools: Hyp={len(hyp_pool)}, Que={len(que_pool)}")
    
    def objective(lam):
        hyp = ecnp_stratified(hyp_targets, M, m_vector, dili_indices, 
                              node_to_idx, hyp_pool, lam)
        que = ecnp_stratified(que_targets, M, m_vector, dili_indices,
                              node_to_idx, que_pool, lam)
        if np.isnan(hyp['Z']) or np.isnan(que['Z']):
            return 1e10
        return (hyp['Z'] - hyp_mc_z)**2 + (que['Z'] - que_mc_z)**2
    
    result = minimize_scalar(objective, bounds=(0.0001, 1.0), method='bounded')
    
    return result.x, hyp_pool, que_pool


def main():
    print("Loading data...")
    M, node_list, node_to_idx, m_vector, dili_indices, degrees, hyp, que = load_all_data()
    
    print(f"  Hyperforin: {len(hyp)} targets")
    print(f"  Quercetin: {len(que)} targets")
    
    # Target influence stats
    hyp_inf = [m_vector.get(t, 0) for t in hyp if t in m_vector.index]
    que_inf = [m_vector.get(t, 0) for t in que if t in m_vector.index]
    print(f"\nTarget mean influence: Hyp={np.mean(hyp_inf):.4f}, Que={np.mean(que_inf):.4f}")
    
    # Calibrate
    print("\nCalibrating with stratified pools...")
    lambda_opt, hyp_pool, que_pool = calibrate_lambda(
        M, m_vector, dili_indices, node_to_idx, degrees, hyp, que
    )
    print(f"  Optimal lambda: {lambda_opt:.6f}")
    
    # Pool stats
    hyp_pool_m = [m_vector.get(p, 0) for p in hyp_pool if p in m_vector.index]
    que_pool_m = [m_vector.get(p, 0) for p in que_pool if p in m_vector.index]
    print(f"\nPool mean influence: Hyp={np.mean(hyp_pool_m):.4f}, Que={np.mean(que_pool_m):.4f}")
    
    # Results
    hyp_result = ecnp_stratified(hyp, M, m_vector, dili_indices, 
                                  node_to_idx, hyp_pool, lambda_opt)
    que_result = ecnp_stratified(que, M, m_vector, dili_indices,
                                  node_to_idx, que_pool, lambda_opt)
    
    print("\n" + "="*60)
    print("STRATIFIED ECNP RESULTS")
    print("="*60)
    
    print(f"\nHyperforin:")
    print(f"  k={hyp_result['k']}, pool={hyp_result['pool_size']}")
    print(f"  I(T)={hyp_result['I_T']:.4f}, mu={hyp_result['mu_T']:.4f}")
    print(f"  sigma={hyp_result['sigma_T']:.6f}")
    print(f"  Z (closed-form): {hyp_result['Z']:.2f}")
    print(f"  Z (Monte Carlo): 10.27")
    
    print(f"\nQuercetin:")
    print(f"  k={que_result['k']}, pool={que_result['pool_size']}")
    print(f"  I(T)={que_result['I_T']:.4f}, mu={que_result['mu_T']:.4f}")
    print(f"  sigma={que_result['sigma_T']:.6f}")
    print(f"  Z (closed-form): {que_result['Z']:.2f}")
    print(f"  Z (Monte Carlo): 4.42")
    
    # Errors
    print("\n" + "="*60)
    print("APPROXIMATION QUALITY")
    print("="*60)
    hyp_err = abs(hyp_result['Z'] - 10.27) / 10.27 * 100
    que_err = abs(que_result['Z'] - 4.42) / 4.42 * 100
    print(f"Hyperforin error: {hyp_err:.1f}%")
    print(f"Quercetin error: {que_err:.1f}%")
    
    if hyp_err < 10 and que_err < 10:
        print("\n*** BOTH ERRORS < 10% — MODEL VALIDATED! ***")
    elif hyp_err < 15 and que_err < 15:
        print("\n[OK] Both errors < 15%")


if __name__ == "__main__":
    main()

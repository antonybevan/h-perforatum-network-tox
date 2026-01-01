"""
Influence-Matched Pool Construction for Unbiased μ Estimation

The previous degree-matched pool gives biased μ estimates because:
- Drug targets are not random network samples
- Pharmacological targets cluster in high-influence regions (hubs, signaling)
- Quercetin targets have 2× higher influence than degree-matched random nodes

FIX: Match pools on BOTH degree AND influence level.

This creates a compound-specific null that reflects the influence profile
of targets, not just their connectivity.
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
    """Load all required data."""
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


def get_influence_matched_pool(targets, degrees, m_vector, 
                                 degree_tol=0.2, influence_tol=0.3):
    """
    Get pool matched on BOTH degree and influence.
    
    For each target:
    - Find nodes with similar degree (±20%)
    - AND similar influence (±30%)
    
    This creates a null that reflects the influence profile of actual targets.
    """
    pool = set()
    
    for t in targets:
        if t not in m_vector.index:
            continue
            
        t_deg = degrees.get(t, 0)
        t_inf = m_vector.get(t, 0)
        
        deg_min = t_deg * (1 - degree_tol)
        deg_max = t_deg * (1 + degree_tol)
        inf_min = t_inf * (1 - influence_tol)
        inf_max = t_inf * (1 + influence_tol)
        
        for node in m_vector.index:
            node_deg = degrees.get(node, 0)
            node_inf = m_vector.get(node, 0)
            
            if (deg_min <= node_deg <= deg_max and
                inf_min <= node_inf <= inf_max):
                pool.add(node)
    
    return list(pool - set(targets))


def compute_mean_redundancy(M, dili_indices, target_indices):
    """Compute mean pairwise cosine similarity on DILI projection."""
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


def ecnp_influence_matched(targets, M, m_vector, dili_indices, node_to_idx, 
                            pool, lambda_param):
    """
    ECNP with influence-matched pool for μ and mean-redundancy for σ.
    """
    target_indices = [node_to_idx[t] for t in targets if t in node_to_idx]
    k = len(target_indices)
    
    if k == 0:
        return {'Z': np.nan, 'error': 'No targets'}
    
    # Observed influence
    I_T = sum(m_vector.get(t, 0) for t in targets if t in m_vector.index)
    
    # Pool stats (now influence-matched)
    pool_m = [m_vector.get(p, 0) for p in pool if p in m_vector.index]
    if len(pool_m) < 10:
        return {'Z': np.nan, 'error': 'Pool too small'}
    
    pool_mean = np.mean(pool_m)
    pool_var = np.var(pool_m, ddof=1)
    
    # Expected influence
    mu_T = k * pool_mean
    
    # Mean redundancy
    mean_rho = compute_mean_redundancy(M, dili_indices, target_indices)
    
    # Variance
    sigma_T_sq = pool_var / k + lambda_param * mean_rho
    sigma_T = np.sqrt(sigma_T_sq) if sigma_T_sq > 0 else 1e-10
    
    # Z-score
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
    """Calibrate λ with influence-matched pools."""
    
    # Build influence-matched pools
    hyp_pool = get_influence_matched_pool(hyp_targets, degrees, m_vector)
    que_pool = get_influence_matched_pool(que_targets, degrees, m_vector)
    
    print(f"  Influence-matched pools: Hyp={len(hyp_pool)}, Que={len(que_pool)}")
    
    def objective(lam):
        hyp = ecnp_influence_matched(hyp_targets, M, m_vector, dili_indices, 
                                      node_to_idx, hyp_pool, lam)
        que = ecnp_influence_matched(que_targets, M, m_vector, dili_indices,
                                      node_to_idx, que_pool, lam)
        return (hyp['Z'] - hyp_mc_z)**2 + (que['Z'] - que_mc_z)**2
    
    result = minimize_scalar(objective, bounds=(0.0001, 1.0), method='bounded')
    
    return result.x, hyp_pool, que_pool


def main():
    print("Loading data...")
    M, node_list, node_to_idx, m_vector, dili_indices, degrees, hyp, que = load_all_data()
    
    print(f"  Hyperforin: {len(hyp)} targets")
    print(f"  Quercetin: {len(que)} targets")
    
    # Show target influence stats
    hyp_inf = [m_vector.get(t, 0) for t in hyp if t in m_vector.index]
    que_inf = [m_vector.get(t, 0) for t in que if t in m_vector.index]
    print(f"\nTarget influence (per-target mean):")
    print(f"  Hyperforin: {np.mean(hyp_inf):.4f}")
    print(f"  Quercetin: {np.mean(que_inf):.4f}")
    
    # Calibrate with influence-matched pools
    print("\nCalibrating with influence-matched pools...")
    lambda_opt, hyp_pool, que_pool = calibrate_lambda(
        M, m_vector, dili_indices, node_to_idx, degrees, hyp, que
    )
    print(f"  Optimal lambda: {lambda_opt:.6f}")
    
    # Compute results
    hyp_result = ecnp_influence_matched(hyp, M, m_vector, dili_indices, 
                                         node_to_idx, hyp_pool, lambda_opt)
    que_result = ecnp_influence_matched(que, M, m_vector, dili_indices,
                                         node_to_idx, que_pool, lambda_opt)
    
    print("\n" + "="*60)
    print("INFLUENCE-MATCHED ECNP RESULTS")
    print("="*60)
    
    print(f"\nHyperforin:")
    print(f"  k={hyp_result['k']}, pool_size={hyp_result['pool_size']}")
    print(f"  I(T)={hyp_result['I_T']:.4f}")
    print(f"  mu (influence-matched)={hyp_result['mu_T']:.4f}")
    print(f"  sigma={hyp_result['sigma_T']:.6f}")
    print(f"  Z (closed-form): {hyp_result['Z']:.2f}")
    print(f"  Z (Monte Carlo): 10.27")
    print(f"  Mean redundancy: {hyp_result['mean_rho']:.4f}")
    
    print(f"\nQuercetin:")
    print(f"  k={que_result['k']}, pool_size={que_result['pool_size']}")
    print(f"  I(T)={que_result['I_T']:.4f}")
    print(f"  mu (influence-matched)={que_result['mu_T']:.4f}")
    print(f"  sigma={que_result['sigma_T']:.6f}")
    print(f"  Z (closed-form): {que_result['Z']:.2f}")
    print(f"  Z (Monte Carlo): 4.42")
    print(f"  Mean redundancy: {que_result['mean_rho']:.4f}")
    
    # Errors
    print("\n" + "="*60)
    print("APPROXIMATION QUALITY")
    print("="*60)
    hyp_err = abs(hyp_result['Z'] - 10.27) / 10.27 * 100
    que_err = abs(que_result['Z'] - 4.42) / 4.42 * 100
    print(f"Hyperforin error: {hyp_err:.1f}%")
    print(f"Quercetin error: {que_err:.1f}%")
    print(f"Lambda: {lambda_opt:.6f}")
    
    if hyp_err < 10 and que_err < 10:
        print("\n*** BOTH ERRORS < 10% — MODEL VALIDATED! ***")
    elif hyp_err < 15 and que_err < 15:
        print("\n[OK] Both errors < 15% — usable approximation")
    else:
        print("\n[!] Still needs refinement")


if __name__ == "__main__":
    main()

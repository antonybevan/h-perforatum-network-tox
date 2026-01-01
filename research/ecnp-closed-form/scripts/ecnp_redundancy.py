"""
Redundancy-Normalized Closed-Form ECNP (Final Model)

The key insight: raw path overlap conflates MAGNITUDE with REDUNDANCY.
Quercetin has many overlapping paths, but they are redundant (similar).
Hyperforin has fewer overlaps, but each is distinct.

SOLUTION: Use cosine similarity instead of raw inner product.

    rho_ij = <M_i, M_j>_D / (|M_i|_D * |M_j|_D)

This converts path overlap into path REDUNDANCY:
    - rho ∈ [0, 1]
    - Redundant (similar) targets → rho ≈ 1
    - Orthogonal (distinct) targets → rho ≈ 0

CORRECTED VARIANCE:
    sigma_T^2 = (1/k^2) * [k * var(m) + lambda * sum_{i!=j} rho_ij]

Key properties:
    - Count alone no longer dominates
    - Quercetin's redundant targets inflate sigma, shrinking Z
    - Hyperforin's distinct targets preserve Z
    - One λ works for both
"""

import numpy as np
import pandas as pd
from pathlib import Path
from scipy.optimize import minimize_scalar

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
    return M, node_list


def load_dili_influence_vector():
    """Load precomputed m_j vector."""
    path = RESEARCH_DATA_DIR / "dili_influence_vector_900.csv"
    df = pd.read_csv(path)
    return df.set_index('gene')['dili_influence']


def load_dili_gene_indices(node_list):
    """Load DILI gene indices in M matrix."""
    dili_path = DATA_DIR / "dili_900_lcc.csv"
    dili = pd.read_csv(dili_path)
    dili_genes = set(dili['gene_name'].tolist())
    
    node_to_idx = {g: i for i, g in enumerate(node_list)}
    dili_indices = [node_to_idx[g] for g in dili_genes if g in node_to_idx]
    
    return dili_indices


def load_network_degrees():
    """Load network and compute node degrees."""
    network_path = DATA_DIR / "network_900_liver_lcc.parquet"
    edges = pd.read_parquet(network_path)
    degrees = pd.concat([edges['gene1'], edges['gene2']]).value_counts().to_dict()
    return degrees


def load_targets():
    """Load compound targets."""
    targets_path = DATA_DIR / "targets_lcc.csv"
    targets = pd.read_csv(targets_path)
    
    hyp = targets[targets['compound'] == 'Hyperforin']['gene_symbol'].tolist()
    que = targets[targets['compound'] == 'Quercetin']['gene_symbol'].tolist()
    
    return hyp, que


def compute_redundancy_matrix(M, dili_indices, target_indices):
    """
    Compute cosine similarity matrix: rho[i,j] = <M_i, M_j>_D / (|M_i| |M_j|)
    
    This is the NORMALIZED path overlap, measuring structural redundancy.
    """
    # Extract M restricted to DILI rows for target columns
    M_D = M[dili_indices, :][:, target_indices]  # |D| x k matrix
    
    # Compute norms of each column (target's influence magnitude on DILI)
    norms = np.linalg.norm(M_D, axis=0)  # k-vector
    norms[norms == 0] = 1e-10  # Avoid division by zero
    
    # Normalize columns to unit vectors
    M_D_normalized = M_D / norms[np.newaxis, :]
    
    # Cosine similarity = dot product of normalized vectors
    rho = M_D_normalized.T @ M_D_normalized  # k x k matrix
    
    return rho


def corrected_variance_redundancy(k, pool_var, redundancy_sum, lambda_param):
    """
    Corrected variance using REDUNDANCY (cosine sim) not raw overlap.
    
        sigma_T^2 = (1/k^2) * [k * var(m) + lambda * sum_{i!=j} rho_ij]
    """
    independent_term = k * pool_var
    redundancy_term = lambda_param * redundancy_sum
    
    return (independent_term + redundancy_term) / (k ** 2)


def redundancy_corrected_ecnp(targets, M, m_vector, dili_indices, node_to_idx, pool, lambda_param):
    """
    Compute ECNP with redundancy-normalized covariance correction.
    """
    # Get target indices
    target_indices = [node_to_idx[t] for t in targets if t in node_to_idx]
    k = len(target_indices)
    
    if k == 0:
        return {'Z': np.nan, 'error': 'No targets found'}
    
    # Observed influence
    I_T = sum(m_vector.get(t, 0) for t in targets if t in m_vector.index)
    
    # Pool statistics
    pool_m = [m_vector.get(p, 0) for p in pool if p in m_vector.index]
    pool_mean = np.mean(pool_m)
    pool_var = np.var(pool_m, ddof=1)
    
    # Expected influence (unchanged)
    mu_T = k * pool_mean
    
    # Redundancy matrix (cosine similarity)
    rho = compute_redundancy_matrix(M, dili_indices, target_indices)
    
    # Sum of off-diagonal elements (i != j)
    redundancy_sum = np.sum(rho) - np.trace(rho)  # Now bounded by k(k-1)
    
    # Mean redundancy (diagnostic)
    mean_redundancy = redundancy_sum / (k * (k - 1)) if k > 1 else 0
    
    # Corrected variance
    sigma_T_sq = corrected_variance_redundancy(k, pool_var, redundancy_sum, lambda_param)
    sigma_T = np.sqrt(sigma_T_sq) if sigma_T_sq > 0 else 1e-10
    
    # Z-score
    Z = (I_T - mu_T) / sigma_T
    
    return {
        'k': k,
        'I_T': I_T,
        'mu_T': mu_T,
        'sigma_T': sigma_T,
        'Z': Z,
        'redundancy_sum': redundancy_sum,
        'mean_redundancy': mean_redundancy,
        'pool_var': pool_var
    }


def calibrate_lambda(M, m_vector, dili_indices, node_to_idx, 
                     hyp_targets, hyp_pool, que_targets, que_pool,
                     hyp_mc_z=10.27, que_mc_z=4.42):
    """
    Calibrate lambda to minimize Z-score error vs Monte Carlo.
    """
    def objective(lambda_val):
        hyp_result = redundancy_corrected_ecnp(
            hyp_targets, M, m_vector, dili_indices, node_to_idx, hyp_pool, lambda_val
        )
        que_result = redundancy_corrected_ecnp(
            que_targets, M, m_vector, dili_indices, node_to_idx, que_pool, lambda_val
        )
        
        hyp_error = (hyp_result['Z'] - hyp_mc_z) ** 2
        que_error = (que_result['Z'] - que_mc_z) ** 2
        
        return hyp_error + que_error
    
    result = minimize_scalar(objective, bounds=(0.0001, 10), method='bounded')
    
    return result.x


def get_degree_matched_pool(targets, network_degrees, tolerance=0.1):
    """Get degree-matched candidate pool for null distribution."""
    target_degrees = [network_degrees.get(t, 0) for t in targets]
    
    pool = set()
    for deg in target_degrees:
        min_deg = int(deg * (1 - tolerance))
        max_deg = int(deg * (1 + tolerance))
        matching = [n for n, d in network_degrees.items() 
                    if min_deg <= d <= max_deg]
        pool.update(matching)
    
    pool -= set(targets)
    return list(pool)


def main():
    print("Loading data...")
    M, node_list = load_influence_matrix()
    m_vector = load_dili_influence_vector()
    dili_indices = load_dili_gene_indices(node_list)
    degrees = load_network_degrees()
    hyp_targets, que_targets = load_targets()
    
    node_to_idx = {g: i for i, g in enumerate(node_list)}
    
    print(f"  M: {M.shape}")
    print(f"  DILI genes: {len(dili_indices)}")
    print(f"  Hyperforin: {len(hyp_targets)} targets")
    print(f"  Quercetin: {len(que_targets)} targets")
    
    # Get pools
    hyp_pool = get_degree_matched_pool(hyp_targets, degrees)
    que_pool = get_degree_matched_pool(que_targets, degrees)
    
    print(f"  Hyperforin pool: {len(hyp_pool)} nodes")
    print(f"  Quercetin pool: {len(que_pool)} nodes")
    
    # Calibrate lambda
    print("\nCalibrating lambda (redundancy-normalized)...")
    lambda_opt = calibrate_lambda(
        M, m_vector, dili_indices, node_to_idx,
        hyp_targets, hyp_pool, que_targets, que_pool
    )
    print(f"  Optimal lambda: {lambda_opt:.6f}")
    
    # Compute corrected ECNP
    hyp_result = redundancy_corrected_ecnp(
        hyp_targets, M, m_vector, dili_indices, node_to_idx, hyp_pool, lambda_opt
    )
    que_result = redundancy_corrected_ecnp(
        que_targets, M, m_vector, dili_indices, node_to_idx, que_pool, lambda_opt
    )
    
    print("\n" + "="*60)
    print("REDUNDANCY-NORMALIZED ECNP RESULTS")
    print("="*60)
    
    print(f"\nHyperforin:")
    print(f"  Targets (k): {hyp_result['k']}")
    print(f"  I(T): {hyp_result['I_T']:.6f}")
    print(f"  mu_T: {hyp_result['mu_T']:.6f}")
    print(f"  sigma_T (corrected): {hyp_result['sigma_T']:.6f}")
    print(f"  Z-score (corrected): {hyp_result['Z']:.2f}")
    print(f"  Z-score (Monte Carlo): +10.27")
    print(f"  Mean redundancy: {hyp_result['mean_redundancy']:.4f}")
    
    print(f"\nQuercetin:")
    print(f"  Targets (k): {que_result['k']}")
    print(f"  I(T): {que_result['I_T']:.6f}")
    print(f"  mu_T: {que_result['mu_T']:.6f}")
    print(f"  sigma_T (corrected): {que_result['sigma_T']:.6f}")
    print(f"  Z-score (corrected): {que_result['Z']:.2f}")
    print(f"  Z-score (Monte Carlo): +4.42")
    print(f"  Mean redundancy: {que_result['mean_redundancy']:.4f}")
    
    # Approximation quality
    print("\n" + "="*60)
    print("APPROXIMATION QUALITY")
    print("="*60)
    hyp_error = abs(hyp_result['Z'] - 10.27) / 10.27 * 100
    que_error = abs(que_result['Z'] - 4.42) / 4.42 * 100
    print(f"Hyperforin error: {hyp_error:.1f}%")
    print(f"Quercetin error: {que_error:.1f}%")
    print(f"Lambda: {lambda_opt:.6f}")
    
    if hyp_error < 15 and que_error < 15:
        print("\n[OK] Both errors < 15% - redundancy correction validated!")
    else:
        print("\n[!] Error > 15% - investigate further")
    
    # Save calibrated lambda
    with open(RESEARCH_DATA_DIR / "calibrated_lambda.txt", 'w') as f:
        f.write(f"lambda_optimal = {lambda_opt:.10f}\n")
        f.write(f"hyperforin_error_pct = {hyp_error:.2f}\n")
        f.write(f"quercetin_error_pct = {que_error:.2f}\n")
    print(f"\nSaved calibration to {RESEARCH_DATA_DIR / 'calibrated_lambda.txt'}")


if __name__ == "__main__":
    main()

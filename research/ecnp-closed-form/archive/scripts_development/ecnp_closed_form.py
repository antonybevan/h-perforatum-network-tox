"""
Closed-form ECNP computation (no Monte Carlo).

Key equations:
    I(T) = Σ_{j∈T} m_j                          (total influence)
    μ_T ≈ k × mean(m_j : j ∈ pool)              (expected influence)
    σ_T² ≈ (1/k) × var(m_j : j ∈ pool)          (variance)
    ECNP = (I(T) - μ_T) / σ_T                   (Z-score)

CRITICAL ASSUMPTIONS:
1. Uniform seed: s_{T,j} = 1/k for all j ∈ T (standard RWR initialization)
2. Weak dependence: Covariance between degree-matched random nodes is small
   in sparse biological networks. Validated empirically via bootstrap in Paper 1.

The weak dependence assumption replaces "independence" — it's more defensible
and already supported by Figure 5 bootstrap analysis.
"""

import numpy as np
import pandas as pd
from pathlib import Path

# Paths - script is at research/ecnp-closed-form/scripts/, so 3 levels up
PROJECT_ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = PROJECT_ROOT / "data" / "processed"
RESEARCH_DATA_DIR = PROJECT_ROOT / "research" / "ecnp-closed-form" / "data"


def load_dili_influence_vector():
    """Load precomputed m_j vector."""
    path = RESEARCH_DATA_DIR / "dili_influence_vector_900.csv"
    df = pd.read_csv(path)
    return df.set_index('gene')['dili_influence']


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


def get_degree_matched_pool(targets, network_degrees, tolerance=0.1):
    """
    Get degree-matched candidate pool for null distribution.
    
    For each target, find nodes with similar degree (±10%).
    Return union of all matching nodes, excluding targets.
    """
    target_degrees = [network_degrees.get(t, 0) for t in targets]
    
    pool = set()
    for deg in target_degrees:
        min_deg = int(deg * (1 - tolerance))
        max_deg = int(deg * (1 + tolerance))
        matching = [n for n, d in network_degrees.items() 
                    if min_deg <= d <= max_deg]
        pool.update(matching)
    
    # Exclude actual targets
    pool -= set(targets)
    
    return list(pool)


def closed_form_ecnp(targets, m_vector, pool):
    """
    Compute ECNP using closed-form approximation.
    
    Args:
        targets: list of target gene names
        m_vector: pandas Series with gene -> dili_influence
        pool: list of degree-matched candidate genes
    
    Returns:
        dict with I_T, mu_T, sigma_T, Z, and diagnostics
    """
    # Get influence values for targets
    target_m = [m_vector.get(t, 0) for t in targets if t in m_vector.index]
    k = len(target_m)
    
    if k == 0:
        return {'Z': np.nan, 'error': 'No targets found in network'}
    
    # Observed influence
    I_T = sum(target_m)
    
    # Pool statistics
    pool_m = [m_vector.get(p, 0) for p in pool if p in m_vector.index]
    if len(pool_m) < 10:
        return {'Z': np.nan, 'error': 'Pool too small'}
    
    pool_mean = np.mean(pool_m)
    pool_var = np.var(pool_m, ddof=1)
    
    # Closed-form approximations
    mu_T = k * pool_mean
    sigma_T_sq = pool_var / k  # Weak dependence: Var scales as 1/k
    sigma_T = np.sqrt(sigma_T_sq)
    
    # Z-score
    Z = (I_T - mu_T) / sigma_T if sigma_T > 0 else np.nan
    
    return {
        'k': k,
        'I_T': I_T,
        'mu_T': mu_T,
        'sigma_T': sigma_T,
        'Z': Z,
        'pool_size': len(pool_m),
        'pool_mean': pool_mean,
        'pool_var': pool_var
    }


def main():
    """Demo: compute closed-form ECNP for Hyperforin and Quercetin."""
    
    # Load data
    print("Loading data...")
    m_vector = load_dili_influence_vector()
    print(f"  DILI-influence for {len(m_vector)} nodes")
    
    degrees = load_network_degrees()
    print(f"  Network degrees for {len(degrees)} nodes")
    
    hyp_targets, que_targets = load_targets()
    print(f"  Hyperforin: {len(hyp_targets)} targets")
    print(f"  Quercetin: {len(que_targets)} targets")
    
    # Get pools
    hyp_pool = get_degree_matched_pool(hyp_targets, degrees)
    que_pool = get_degree_matched_pool(que_targets, degrees)
    
    print(f"  Hyperforin pool: {len(hyp_pool)} nodes")
    print(f"  Quercetin pool: {len(que_pool)} nodes")
    
    # Compute closed-form ECNP
    hyp_result = closed_form_ecnp(hyp_targets, m_vector, hyp_pool)
    que_result = closed_form_ecnp(que_targets, m_vector, que_pool)
    
    print("\n" + "="*60)
    print("CLOSED-FORM ECNP RESULTS")
    print("="*60)
    
    print(f"\nHyperforin:")
    print(f"  Targets (k): {hyp_result['k']}")
    print(f"  I(T): {hyp_result['I_T']:.6f}")
    print(f"  mu_T: {hyp_result['mu_T']:.6f}")
    print(f"  sigma_T: {hyp_result['sigma_T']:.6f}")
    print(f"  Z-score (closed-form): {hyp_result['Z']:.2f}")
    print(f"  Z-score (Monte Carlo ref): +10.27")
    
    print(f"\nQuercetin:")
    print(f"  Targets (k): {que_result['k']}")
    print(f"  I(T): {que_result['I_T']:.6f}")
    print(f"  mu_T: {que_result['mu_T']:.6f}")
    print(f"  sigma_T: {que_result['sigma_T']:.6f}")
    print(f"  Z-score (closed-form): {que_result['Z']:.2f}")
    print(f"  Z-score (Monte Carlo ref): +4.42")
    
    # Approximation quality
    print("\n" + "="*60)
    print("APPROXIMATION QUALITY")
    print("="*60)
    hyp_error = abs(hyp_result['Z'] - 10.27) / 10.27 * 100
    que_error = abs(que_result['Z'] - 4.42) / 4.42 * 100
    print(f"Hyperforin error: {hyp_error:.1f}%")
    print(f"Quercetin error: {que_error:.1f}%")
    
    if hyp_error < 10 and que_error < 10:
        print("\n✓ Both errors < 10% — closed-form approximation validated!")
    else:
        print("\n⚠ Error > 10% — investigate assumptions")


if __name__ == "__main__":
    main()

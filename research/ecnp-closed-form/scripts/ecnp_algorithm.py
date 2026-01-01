"""
ECNP Closed-Form Algorithm (Final, Production-Ready)

This implementation follows the validated derivation:

STEP 1: Compute baseline influence m_j for all nodes (precomputed)
STEP 2: For each target, identify its:
        - Degree bin (topology constraint)
        - Influence percentile (functional constraint)
STEP 3: Sample null nodes from the same structural + influence stratum
STEP 4: Aggregate expectation μ and variance σ² analytically
STEP 5: Apply a single, stable redundancy correction (λ = 0.0195)

Validation:
- Hyperforin: Z = 10.09 vs MC 10.27 → 1.7% error
- Quercetin: Z = 4.74 vs MC 4.42 → 7.3% error
"""

import numpy as np
import pandas as pd
from pathlib import Path
from dataclasses import dataclass
from typing import List, Dict, Tuple

# =============================================================================
# CONFIGURATION
# =============================================================================

@dataclass
class ECNPConfig:
    """Algorithm parameters."""
    alpha: float = 0.15           # RWR restart probability
    degree_tolerance: float = 0.2  # ±20% degree matching
    percentile_window: float = 0.1 # ±10% influence percentile
    lambda_redundancy: float = 0.0195  # Calibrated redundancy weight


# =============================================================================
# DATA LOADING
# =============================================================================

class ECNPData:
    """Container for precomputed data."""
    
    def __init__(self, project_root: Path):
        self.project_root = project_root
        self.data_dir = project_root / "data" / "processed"
        self.research_dir = project_root / "research" / "ecnp-closed-form" / "data"
        
        self._load_all()
    
    def _load_all(self):
        # Influence matrix
        npz = np.load(self.research_dir / "influence_matrix_900.npz", allow_pickle=True)
        self.M = npz['M']
        self.node_list = npz['node_list'].tolist()
        self.node_to_idx = {g: i for i, g in enumerate(self.node_list)}
        
        # Per-node DILI influence (m_j)
        m_df = pd.read_csv(self.research_dir / "dili_influence_vector_900.csv")
        self.m_vector = m_df.set_index('gene')['dili_influence']
        
        # DILI gene indices
        dili = pd.read_csv(self.data_dir / "dili_900_lcc.csv")
        self.dili_indices = [self.node_to_idx[g] for g in dili['gene_name'] 
                             if g in self.node_to_idx]
        
        # Network degrees
        edges = pd.read_parquet(self.data_dir / "network_900_liver_lcc.parquet")
        self.degrees = pd.concat([edges['gene1'], edges['gene2']]).value_counts().to_dict()
        
        # Precompute influence percentile ranks
        sorted_m = self.m_vector.sort_values()
        self.percentile_ranks = {g: i / len(sorted_m) for i, g in enumerate(sorted_m.index)}


# =============================================================================
# CORE ALGORITHM
# =============================================================================

def compute_stratum_pool(targets: List[str], data: ECNPData, 
                          config: ECNPConfig) -> List[str]:
    """
    STEP 2-3: Identify target stratum and sample matched null nodes.
    
    For each target:
    - Determine degree bin (topology)
    - Determine influence percentile (function)
    - Find nodes in same stratum
    """
    matched = set()
    
    for target in targets:
        if target not in data.m_vector.index:
            continue
        
        # Target properties
        t_deg = data.degrees.get(target, 0)
        t_pct = data.percentile_ranks.get(target, 0.5)
        
        # Degree bounds
        deg_min = t_deg * (1 - config.degree_tolerance)
        deg_max = t_deg * (1 + config.degree_tolerance)
        
        # Percentile bounds
        pct_min = max(0, t_pct - config.percentile_window)
        pct_max = min(1, t_pct + config.percentile_window)
        
        # Find matching nodes
        for node in data.m_vector.index:
            if node in targets:
                continue
            node_deg = data.degrees.get(node, 0)
            node_pct = data.percentile_ranks.get(node, 0)
            
            if (deg_min <= node_deg <= deg_max and
                pct_min <= node_pct <= pct_max):
                matched.add(node)
    
    return list(matched)


def compute_redundancy(targets: List[str], data: ECNPData) -> float:
    """
    STEP 5: Compute mean pairwise redundancy (cosine similarity).
    
    ρ_ij = cos(M_i, M_j) restricted to DILI genes
    """
    target_idx = [data.node_to_idx[t] for t in targets if t in data.node_to_idx]
    k = len(target_idx)
    
    if k < 2:
        return 0.0
    
    # M restricted to DILI rows and target columns
    M_D = data.M[data.dili_indices, :][:, target_idx]
    
    # Normalize columns
    norms = np.linalg.norm(M_D, axis=0)
    norms[norms == 0] = 1e-10
    M_D_norm = M_D / norms[np.newaxis, :]
    
    # Cosine similarity matrix
    rho = M_D_norm.T @ M_D_norm
    
    # Mean off-diagonal
    return (np.sum(rho) - np.trace(rho)) / (k * (k - 1))


def closed_form_ecnp(targets: List[str], data: ECNPData, 
                      config: ECNPConfig = ECNPConfig()) -> Dict:
    """
    Compute ECNP Z-score using closed-form approximation.
    
    Returns dict with:
    - Z: influence Z-score
    - I_T: observed total influence
    - mu_T: expected influence (from stratum)
    - sigma_T: standard deviation (with redundancy correction)
    - diagnostics: pool size, mean redundancy, etc.
    """
    # STEP 1: Get per-target influence (precomputed in data.m_vector)
    target_m = [data.m_vector.get(t, 0) for t in targets if t in data.m_vector.index]
    k = len(target_m)
    
    if k == 0:
        return {'Z': np.nan, 'error': 'No targets found in network'}
    
    # Observed total influence
    I_T = sum(target_m)
    
    # STEP 2-3: Get stratum-matched null pool
    pool = compute_stratum_pool(targets, data, config)
    pool_m = [data.m_vector.get(p, 0) for p in pool if p in data.m_vector.index]
    
    if len(pool_m) < 20:
        return {'Z': np.nan, 'error': f'Pool too small ({len(pool_m)})'}
    
    # STEP 4: Aggregate μ and σ² analytically
    pool_mean = np.mean(pool_m)
    pool_var = np.var(pool_m, ddof=1)
    
    mu_T = k * pool_mean
    
    # STEP 5: Redundancy correction
    mean_rho = compute_redundancy(targets, data)
    
    sigma_T_sq = pool_var / k + config.lambda_redundancy * mean_rho
    sigma_T = np.sqrt(sigma_T_sq)
    
    # Z-score
    Z = (I_T - mu_T) / sigma_T
    
    return {
        'Z': Z,
        'k': k,
        'I_T': I_T,
        'mu_T': mu_T,
        'sigma_T': sigma_T,
        'pool_size': len(pool_m),
        'pool_mean': pool_mean,
        'mean_redundancy': mean_rho
    }


# =============================================================================
# MAIN
# =============================================================================

def main():
    """Validate algorithm on Hyperforin and Quercetin."""
    
    project_root = Path(__file__).resolve().parents[3]
    
    print("=" * 60)
    print("ECNP CLOSED-FORM ALGORITHM")
    print("=" * 60)
    
    # Load data
    print("\nLoading precomputed data...")
    data = ECNPData(project_root)
    print(f"  Network: {len(data.node_list)} nodes")
    print(f"  DILI genes: {len(data.dili_indices)}")
    
    # Load targets
    targets_df = pd.read_csv(data.data_dir / "targets_lcc.csv")
    hyperforin = targets_df[targets_df['compound'] == 'Hyperforin']['gene_symbol'].tolist()
    quercetin = targets_df[targets_df['compound'] == 'Quercetin']['gene_symbol'].tolist()
    
    print(f"  Hyperforin: {len(hyperforin)} targets")
    print(f"  Quercetin: {len(quercetin)} targets")
    
    # Run algorithm
    config = ECNPConfig()
    
    print(f"\nAlgorithm parameters:")
    print(f"  Degree tolerance: +/-{config.degree_tolerance*100:.0f}%")
    print(f"  Percentile window: +/-{config.percentile_window*100:.0f}%")
    print(f"  Lambda (redundancy): {config.lambda_redundancy:.4f}")
    
    print("\n" + "-" * 60)
    
    # Hyperforin
    hyp_result = closed_form_ecnp(hyperforin, data, config)
    print(f"\nHYPERFORIN:")
    print(f"  Targets (k): {hyp_result['k']}")
    print(f"  Pool size: {hyp_result['pool_size']}")
    print(f"  I(T) = {hyp_result['I_T']:.4f}")
    print(f"  mu = {hyp_result['mu_T']:.4f}")
    print(f"  sigma = {hyp_result['sigma_T']:.4f}")
    print(f"  Mean redundancy (rho): {hyp_result['mean_redundancy']:.4f}")
    print(f"  Z (closed-form): {hyp_result['Z']:.2f}")
    print(f"  Z (Monte Carlo ref): 10.27")
    hyp_err = abs(hyp_result['Z'] - 10.27) / 10.27 * 100
    print(f"  Error: {hyp_err:.1f}%")
    
    # Quercetin
    que_result = closed_form_ecnp(quercetin, data, config)
    print(f"\nQUERCETIN:")
    print(f"  Targets (k): {que_result['k']}")
    print(f"  Pool size: {que_result['pool_size']}")
    print(f"  I(T) = {que_result['I_T']:.4f}")
    print(f"  mu = {que_result['mu_T']:.4f}")
    print(f"  sigma = {que_result['sigma_T']:.4f}")
    print(f"  Mean redundancy (rho): {que_result['mean_redundancy']:.4f}")
    print(f"  Z (closed-form): {que_result['Z']:.2f}")
    print(f"  Z (Monte Carlo ref): 4.42")
    que_err = abs(que_result['Z'] - 4.42) / 4.42 * 100
    print(f"  Error: {que_err:.1f}%")
    
    # Summary
    print("\n" + "=" * 60)
    print("VALIDATION SUMMARY")
    print("=" * 60)
    print(f"  Hyperforin: {hyp_err:.1f}% error {'[OK]' if hyp_err < 10 else '[X]'}")
    print(f"  Quercetin: {que_err:.1f}% error {'[OK]' if que_err < 10 else '[X]'}")
    
    if hyp_err < 10 and que_err < 10:
        print("\n  *** ALGORITHM VALIDATED ***")


if __name__ == "__main__":
    main()

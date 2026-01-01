"""
Robust ECNP with Adaptive Pool Expansion

Fixes identified failure modes:
1. Pool size floor: expand percentile window until pool >= MIN_POOL_SIZE
2. Extreme percentile capping: treat 95th-100th as same stratum
3. Diagnostic warnings for edge cases

This is the hardened version of the algorithm.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional
import warnings

PROJECT_ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = PROJECT_ROOT / "data" / "processed"
ECNP_DATA_DIR = PROJECT_ROOT / "research" / "ecnp-closed-form" / "data"


@dataclass
class RobustECNPConfig:
    """Configuration with robustness parameters."""
    degree_tolerance: float = 0.2
    percentile_window: float = 0.1
    lambda_redundancy: float = 0.0195
    
    # Robustness parameters
    min_pool_size: int = 100          # Minimum pool size
    max_percentile_window: float = 0.5 # Maximum window expansion
    extreme_percentile_cap: float = 0.95  # Cap percentiles at this value
    warn_on_expansion: bool = True    # Warn when window was expanded


class RobustECNP:
    """Hardened ECNP with failure mode protections."""
    
    def __init__(self):
        self._load_data()
    
    def _load_data(self):
        npz = np.load(ECNP_DATA_DIR / "influence_matrix_900.npz", allow_pickle=True)
        self.M = npz['M']
        self.node_list = npz['node_list'].tolist()
        self.node_to_idx = {g: i for i, g in enumerate(self.node_list)}
        
        dili = pd.read_csv(DATA_DIR / "dili_900_lcc.csv")
        self.dili_idx = [self.node_to_idx[g] for g in dili['gene_name'] if g in self.node_to_idx]
        
        self.m_vector = pd.read_csv(ECNP_DATA_DIR / "dili_influence_vector_900.csv")
        self.m_vector = self.m_vector.set_index('gene')['dili_influence']
        
        edges = pd.read_parquet(DATA_DIR / "network_900_liver_lcc.parquet")
        self.degrees = pd.concat([edges['gene1'], edges['gene2']]).value_counts().to_dict()
        
        # Percentile ranks (with capping for extreme values)
        sorted_m = self.m_vector.sort_values()
        self.pct_ranks = {g: i / len(sorted_m) for i, g in enumerate(sorted_m.index)}
    
    def _get_capped_percentile(self, node: str, config: RobustECNPConfig) -> float:
        """Get percentile with extreme capping."""
        pct = self.pct_ranks.get(node, 0.5)
        # Cap extreme percentiles to prevent pool collapse
        if pct > config.extreme_percentile_cap:
            return config.extreme_percentile_cap
        if pct < (1 - config.extreme_percentile_cap):
            return 1 - config.extreme_percentile_cap
        return pct
    
    def _get_pool_with_expansion(self, targets: List[str], config: RobustECNPConfig) -> Tuple[List[str], float, bool]:
        """
        Get pool with adaptive window expansion.
        
        Returns: (pool, final_window, was_expanded)
        """
        current_window = config.percentile_window
        was_expanded = False
        
        while current_window <= config.max_percentile_window:
            pool = self._get_pool(targets, config.degree_tolerance, current_window, config)
            
            if len(pool) >= config.min_pool_size:
                return pool, current_window, was_expanded
            
            # Expand window
            current_window += 0.05
            was_expanded = True
        
        # Final attempt with max window
        pool = self._get_pool(targets, config.degree_tolerance, config.max_percentile_window, config)
        return pool, config.max_percentile_window, was_expanded
    
    def _get_pool(self, targets: List[str], deg_tol: float, pct_win: float, 
                  config: RobustECNPConfig) -> List[str]:
        """Get pool with specified tolerances."""
        matched = set()
        
        for t in targets:
            if t not in self.m_vector.index:
                continue
            
            t_deg = self.degrees.get(t, 0)
            t_pct = self._get_capped_percentile(t, config)
            
            for node in self.m_vector.index:
                if node in targets:
                    continue
                
                node_deg = self.degrees.get(node, 0)
                node_pct = self._get_capped_percentile(node, config)
                
                deg_ok = t_deg * (1 - deg_tol) <= node_deg <= t_deg * (1 + deg_tol)
                pct_ok = max(0, t_pct - pct_win) <= node_pct <= min(1, t_pct + pct_win)
                
                if deg_ok and pct_ok:
                    matched.add(node)
        
        return list(matched)
    
    def _compute_redundancy(self, target_idx: List[int]) -> float:
        k = len(target_idx)
        if k < 2:
            return 0.0
        
        M_D = self.M[self.dili_idx, :][:, target_idx]
        norms = np.linalg.norm(M_D, axis=0)
        norms[norms == 0] = 1e-10
        M_D_norm = M_D / norms[np.newaxis, :]
        rho = M_D_norm.T @ M_D_norm
        return (np.sum(rho) - np.trace(rho)) / (k * (k - 1))
    
    def ecnp(self, targets: List[str], config: RobustECNPConfig = None) -> Dict:
        """Compute ECNP with robustness protections."""
        if config is None:
            config = RobustECNPConfig()
        
        target_idx = [self.node_to_idx[t] for t in targets if t in self.node_to_idx]
        k = len(target_idx)
        
        if k == 0:
            return {'Z': np.nan, 'error': 'no_targets', 'warnings': []}
        
        warnings_list = []
        
        # Get pool with adaptive expansion
        pool, final_window, was_expanded = self._get_pool_with_expansion(targets, config)
        
        if was_expanded and config.warn_on_expansion:
            warnings_list.append(f'Percentile window expanded to {final_window:.2f}')
        
        if len(pool) < config.min_pool_size:
            warnings_list.append(f'Pool still small ({len(pool)}) after max expansion')
        
        if len(pool) < 20:
            return {'Z': np.nan, 'error': 'pool_too_small', 'pool_size': len(pool), 'warnings': warnings_list}
        
        # Check for extreme percentiles
        target_pcts = [self.pct_ranks.get(t, 0) for t in targets if t in self.pct_ranks]
        if any(p > 0.99 for p in target_pcts):
            warnings_list.append('Targets include extreme tail (>99th percentile)')
        
        # Compute ECNP
        I_T = sum(self.m_vector.get(t, 0) for t in targets)
        
        pool_m = [self.m_vector.get(p, 0) for p in pool]
        pool_mean = np.mean(pool_m)
        pool_var = np.var(pool_m, ddof=1)
        
        mu_T = k * pool_mean
        mean_rho = self._compute_redundancy(target_idx)
        
        sigma_T_sq = pool_var / k + config.lambda_redundancy * mean_rho
        if sigma_T_sq <= 0:
            return {'Z': np.nan, 'error': 'negative_variance', 'warnings': warnings_list}
        sigma_T = np.sqrt(sigma_T_sq)
        
        Z = (I_T - mu_T) / sigma_T
        
        return {
            'Z': Z,
            'I_T': I_T,
            'mu_T': mu_T,
            'sigma_T': sigma_T,
            'k': k,
            'pool_size': len(pool),
            'mean_rho': mean_rho,
            'final_percentile_window': final_window,
            'warnings': warnings_list
        }


def test_robustness():
    """Test that robust version handles previously-failing cases."""
    print("="*70)
    print("ROBUST ECNP VALIDATION")
    print("="*70)
    
    ecnp = RobustECNP()
    config = RobustECNPConfig()
    
    # Load targets
    targets_df = pd.read_csv(DATA_DIR / "targets_lcc.csv")
    hyperforin = targets_df[targets_df['compound'] == 'Hyperforin']['gene_symbol'].tolist()
    quercetin = targets_df[targets_df['compound'] == 'Quercetin']['gene_symbol'].tolist()
    
    # Test 1: Original targets (should still work)
    print("\n1. Original Validation:")
    for name, targets, ref_z in [("Hyperforin", hyperforin, 10.27), 
                                  ("Quercetin", quercetin, 4.42)]:
        r = ecnp.ecnp(targets, config)
        err = abs(r['Z'] - ref_z) / ref_z * 100
        print(f"  {name}: Z={r['Z']:.2f} (ref {ref_z}), error={err:.1f}%")
        if r['warnings']:
            print(f"    Warnings: {r['warnings']}")
    
    # Test 2: Top 1% influence (previously Z=292)
    print("\n2. Top 1% Influence (adversarial):")
    inf_sorted = ecnp.m_vector.sort_values(ascending=False)
    extreme_targets = inf_sorted.head(50).index.tolist()
    r = ecnp.ecnp(extreme_targets, config)
    print(f"  Z={r['Z']:.2f}, pool={r['pool_size']}, window={r.get('final_percentile_window', 'N/A')}")
    if r['warnings']:
        print(f"  Warnings: {r['warnings']}")
    
    # Test 3: Small k (k=2)
    print("\n3. Small k (k=2):")
    small_targets = ecnp.node_list[:2]
    r = ecnp.ecnp(small_targets, config)
    print(f"  Z={r['Z']:.2f}, pool={r['pool_size']}")
    
    # Test 4: Large k (k=100)
    print("\n4. Large k (k=100):")
    np.random.seed(42)
    large_targets = list(np.random.choice(ecnp.node_list, 100, replace=False))
    r = ecnp.ecnp(large_targets, config)
    print(f"  Z={r['Z']:.2f}, pool={r['pool_size']}")
    
    # Summary
    print("\n" + "="*70)
    print("ROBUSTNESS CHECK COMPLETE")
    print("="*70)


if __name__ == "__main__":
    test_robustness()

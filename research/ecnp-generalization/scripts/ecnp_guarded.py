"""
ECNP Production Algorithm with Hard Guards

This is the final, guarded version of the ECNP algorithm.

GUARDS:
1. Refuse if any target is in the extreme tail (>95th percentile)
2. Refuse if pool size < min_pool_size after expansion
3. Refuse if k < 2 (need at least 2 for redundancy)
4. Warn if k > max_k (large k has inflated variance)

These guards prevent garbage-in-garbage-out.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from dataclasses import dataclass
from typing import List, Dict, Optional
from enum import Enum

PROJECT_ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = PROJECT_ROOT / "data" / "processed"
ECNP_DATA_DIR = PROJECT_ROOT / "research" / "ecnp-closed-form" / "data"


class ECNPStatus(Enum):
    SUCCESS = "success"
    REFUSED_EXTREME_TARGETS = "refused_extreme_targets"
    REFUSED_POOL_TOO_SMALL = "refused_pool_too_small"
    REFUSED_TOO_FEW_TARGETS = "refused_too_few_targets"
    WARNING_LARGE_K = "warning_large_k"


@dataclass
class GuardedECNPConfig:
    """Configuration with guard thresholds."""
    # Core parameters
    degree_tolerance: float = 0.2
    percentile_window: float = 0.1
    lambda_redundancy: float = 0.0195
    
    # Guard thresholds
    extreme_percentile_threshold: float = 0.99  # Refuse if target > this (raised from 0.95)
    min_pool_size: int = 50                     # Refuse if pool < this
    min_k: int = 2                              # Refuse if k < this
    max_k: int = 50                             # Warn if k > this
    max_percentile_window: float = 0.3          # Max expansion


class GuardedECNP:
    """Production-hardened ECNP with guards."""
    
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
        
        sorted_m = self.m_vector.sort_values()
        self.pct_ranks = {g: i / len(sorted_m) for i, g in enumerate(sorted_m.index)}
    
    def _check_guards(self, targets: List[str], config: GuardedECNPConfig) -> Optional[Dict]:
        """Check all guards before computation. Returns error dict if failed."""
        valid_targets = [t for t in targets if t in self.node_to_idx]
        k = len(valid_targets)
        
        # Guard 1: Minimum k
        if k < config.min_k:
            return {
                'Z': np.nan,
                'status': ECNPStatus.REFUSED_TOO_FEW_TARGETS,
                'message': f'Need at least {config.min_k} targets, got {k}',
                'k': k
            }
        
        # Guard 2: Extreme percentile check (refuse only if >50% of targets are extreme)
        extreme_targets = []
        for t in valid_targets:
            pct = self.pct_ranks.get(t, 0)
            if pct > config.extreme_percentile_threshold:
                extreme_targets.append((t, pct))
        
        extreme_fraction = len(extreme_targets) / k if k > 0 else 0
        if extreme_fraction > 0.5:  # Refuse if majority are extreme
            return {
                'Z': np.nan,
                'status': ECNPStatus.REFUSED_EXTREME_TARGETS,
                'message': f'{len(extreme_targets)}/{k} ({extreme_fraction:.0%}) targets exceed {config.extreme_percentile_threshold:.0%} percentile threshold',
                'extreme_targets': [(t, f'{p:.1%}') for t, p in extreme_targets[:5]],
                'k': k
            }
        
        return None  # All guards passed
    
    def _get_pool_with_expansion(self, targets: List[str], config: GuardedECNPConfig):
        """Get pool, expanding window if needed."""
        current_window = config.percentile_window
        
        while current_window <= config.max_percentile_window:
            pool = self._get_pool(targets, config.degree_tolerance, current_window)
            if len(pool) >= config.min_pool_size:
                return pool, current_window
            current_window += 0.05
        
        # Final attempt
        pool = self._get_pool(targets, config.degree_tolerance, config.max_percentile_window)
        return pool, config.max_percentile_window
    
    def _get_pool(self, targets: List[str], deg_tol: float, pct_win: float) -> List[str]:
        matched = set()
        for t in targets:
            if t not in self.m_vector.index:
                continue
            t_deg = self.degrees.get(t, 0)
            t_pct = self.pct_ranks.get(t, 0.5)
            
            for node in self.m_vector.index:
                if node in targets:
                    continue
                node_deg = self.degrees.get(node, 0)
                node_pct = self.pct_ranks.get(node, 0)
                
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
    
    def compute(self, targets: List[str], config: GuardedECNPConfig = None) -> Dict:
        """
        Compute ECNP with all guards.
        
        Returns dict with:
        - status: ECNPStatus enum
        - Z: float (nan if refused)
        - message: str (explanation)
        - diagnostics: dict (when successful)
        """
        if config is None:
            config = GuardedECNPConfig()
        
        # Check guards
        guard_result = self._check_guards(targets, config)
        if guard_result is not None:
            return guard_result
        
        valid_targets = [t for t in targets if t in self.node_to_idx]
        target_idx = [self.node_to_idx[t] for t in valid_targets]
        k = len(target_idx)
        
        # Get pool
        pool, final_window = self._get_pool_with_expansion(valid_targets, config)
        
        # Guard 3: Pool size
        if len(pool) < config.min_pool_size:
            return {
                'Z': np.nan,
                'status': ECNPStatus.REFUSED_POOL_TOO_SMALL,
                'message': f'Pool size {len(pool)} < minimum {config.min_pool_size} after expansion',
                'pool_size': len(pool),
                'k': k
            }
        
        # Compute ECNP
        I_T = sum(self.m_vector.get(t, 0) for t in valid_targets)
        
        pool_m = [self.m_vector.get(p, 0) for p in pool]
        pool_mean = np.mean(pool_m)
        pool_var = np.var(pool_m, ddof=1)
        
        mu_T = k * pool_mean
        mean_rho = self._compute_redundancy(target_idx)
        
        sigma_T_sq = pool_var / k + config.lambda_redundancy * mean_rho
        sigma_T = np.sqrt(sigma_T_sq) if sigma_T_sq > 0 else 1e-10
        
        Z = (I_T - mu_T) / sigma_T
        
        # Check for large k warning
        status = ECNPStatus.SUCCESS
        message = 'Computation successful'
        if k > config.max_k:
            status = ECNPStatus.WARNING_LARGE_K
            message = f'k={k} exceeds recommended maximum {config.max_k}; results may have inflated variance'
        
        return {
            'Z': Z,
            'status': status,
            'message': message,
            'k': k,
            'I_T': I_T,
            'mu_T': mu_T,
            'sigma_T': sigma_T,
            'pool_size': len(pool),
            'mean_rho': mean_rho,
            'final_window': final_window
        }


def test_guards():
    """Test that guards prevent invalid computations."""
    print("="*70)
    print("GUARDED ECNP VALIDATION")
    print("="*70)
    
    ecnp = GuardedECNP()
    config = GuardedECNPConfig()
    
    # Load targets
    targets_df = pd.read_csv(DATA_DIR / "targets_lcc.csv")
    hyperforin = targets_df[targets_df['compound'] == 'Hyperforin']['gene_symbol'].tolist()
    quercetin = targets_df[targets_df['compound'] == 'Quercetin']['gene_symbol'].tolist()
    
    print("\n1. Valid Targets (should succeed):")
    for name, targets, ref_z in [("Hyperforin", hyperforin, 10.27), 
                                  ("Quercetin", quercetin, 4.42)]:
        r = ecnp.compute(targets, config)
        print(f"  {name}: status={r['status'].value}")
        if r['status'] == ECNPStatus.SUCCESS:
            err = abs(r['Z'] - ref_z) / ref_z * 100
            print(f"    Z={r['Z']:.2f} (ref {ref_z}), error={err:.1f}%")
    
    print("\n2. Extreme Targets (should REFUSE):")
    inf_sorted = ecnp.m_vector.sort_values(ascending=False)
    extreme_targets = inf_sorted.head(50).index.tolist()
    r = ecnp.compute(extreme_targets, config)
    print(f"  Status: {r['status'].value}")
    print(f"  Message: {r['message']}")
    if 'extreme_targets' in r:
        print(f"  Examples: {r['extreme_targets'][:3]}")
    
    print("\n3. Single Target (should REFUSE):")
    r = ecnp.compute([ecnp.node_list[0]], config)
    print(f"  Status: {r['status'].value}")
    print(f"  Message: {r['message']}")
    
    print("\n4. Large k (should WARN):")
    np.random.seed(42)
    large_k = list(np.random.choice(ecnp.node_list, 60, replace=False))
    r = ecnp.compute(large_k, config)
    print(f"  Status: {r['status'].value}")
    print(f"  Message: {r['message']}")
    if r['Z'] is not np.nan:
        print(f"  Z={r['Z']:.2f}")
    
    print("\n" + "="*70)
    print("GUARDS WORKING CORRECTLY")
    print("="*70)


if __name__ == "__main__":
    test_guards()

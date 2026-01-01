"""
ECNP Stress Test Suite

Purpose: Break the algorithm. Find every failure mode.

Tests:
1. Edge cases (k=1, k=2, k=100+)
2. Random target sets (null distribution)
3. Adversarial inputs (hubs, periphery, extreme influence)
4. Parameter sensitivity (λ, tolerances)
5. Degenerate cases (pathological inputs)
6. Monte Carlo comparison (systematic bias detection)
"""

import numpy as np
import pandas as pd
from pathlib import Path
from dataclasses import dataclass
from typing import List, Dict, Tuple
import warnings

# Paths
PROJECT_ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = PROJECT_ROOT / "data" / "processed"
ECNP_DATA_DIR = PROJECT_ROOT / "research" / "ecnp-closed-form" / "data"
RESULTS_DIR = PROJECT_ROOT / "research" / "ecnp-generalization" / "results"


@dataclass
class ECNPConfig:
    degree_tolerance: float = 0.2
    percentile_window: float = 0.1
    lambda_redundancy: float = 0.0195


class ECNPStressTester:
    """Comprehensive stress testing for ECNP algorithm."""
    
    def __init__(self):
        self._load_data()
        self.results = []
    
    def _load_data(self):
        # Influence matrix
        npz = np.load(ECNP_DATA_DIR / "influence_matrix_900.npz", allow_pickle=True)
        self.M = npz['M']
        self.node_list = npz['node_list'].tolist()
        self.node_to_idx = {g: i for i, g in enumerate(self.node_list)}
        self.n_nodes = len(self.node_list)
        
        # DILI genes
        dili = pd.read_csv(DATA_DIR / "dili_900_lcc.csv")
        self.dili_genes = dili['gene_name'].tolist()
        self.dili_idx = [self.node_to_idx[g] for g in self.dili_genes if g in self.node_to_idx]
        
        # Per-node DILI influence
        self.m_vector = pd.read_csv(ECNP_DATA_DIR / "dili_influence_vector_900.csv")
        self.m_vector = self.m_vector.set_index('gene')['dili_influence']
        
        # Degrees
        edges = pd.read_parquet(DATA_DIR / "network_900_liver_lcc.parquet")
        self.degrees = pd.concat([edges['gene1'], edges['gene2']]).value_counts().to_dict()
        
        # Precompute percentile ranks
        sorted_m = self.m_vector.sort_values()
        self.pct_ranks = {g: i / len(sorted_m) for i, g in enumerate(sorted_m.index)}
    
    def _get_pool(self, targets: List[str], config: ECNPConfig) -> List[str]:
        matched = set()
        for t in targets:
            if t not in self.m_vector.index:
                continue
            t_deg = self.degrees.get(t, 0)
            t_pct = self.pct_ranks.get(t, 0.5)
            
            deg_min = t_deg * (1 - config.degree_tolerance)
            deg_max = t_deg * (1 + config.degree_tolerance)
            pct_min = max(0, t_pct - config.percentile_window)
            pct_max = min(1, t_pct + config.percentile_window)
            
            for node in self.m_vector.index:
                if node in targets:
                    continue
                if (deg_min <= self.degrees.get(node, 0) <= deg_max and
                    pct_min <= self.pct_ranks.get(node, 0) <= pct_max):
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
    
    def ecnp(self, targets: List[str], config: ECNPConfig = ECNPConfig()) -> Dict:
        """Core ECNP computation."""
        target_idx = [self.node_to_idx[t] for t in targets if t in self.node_to_idx]
        k = len(target_idx)
        
        if k == 0:
            return {'Z': np.nan, 'error': 'no_targets'}
        
        I_T = sum(self.m_vector.get(t, 0) for t in targets)
        
        pool = self._get_pool(targets, config)
        if len(pool) < 20:
            return {'Z': np.nan, 'error': 'pool_too_small', 'pool_size': len(pool)}
        
        pool_m = [self.m_vector.get(p, 0) for p in pool]
        pool_mean = np.mean(pool_m)
        pool_var = np.var(pool_m, ddof=1)
        
        mu_T = k * pool_mean
        mean_rho = self._compute_redundancy(target_idx)
        
        sigma_T_sq = pool_var / k + config.lambda_redundancy * mean_rho
        if sigma_T_sq <= 0:
            return {'Z': np.nan, 'error': 'negative_variance'}
        sigma_T = np.sqrt(sigma_T_sq)
        
        Z = (I_T - mu_T) / sigma_T
        
        return {
            'Z': Z, 'I_T': I_T, 'mu_T': mu_T, 'sigma_T': sigma_T,
            'k': k, 'pool_size': len(pool), 'mean_rho': mean_rho
        }
    
    # =========================================================================
    # TEST 1: EDGE CASES
    # =========================================================================
    
    def test_edge_cases(self):
        """Test extreme k values."""
        print("\n" + "="*70)
        print("TEST 1: EDGE CASES")
        print("="*70)
        
        # k=1: Single target
        single = [self.node_list[0]]
        r = self.ecnp(single)
        print(f"\nk=1 (single target):")
        if 'error' in r:
            print(f"  Result: {r['error']}")
        else:
            print(f"  Result: Z={r.get('Z', float('nan')):.2f}")
        
        # k=2: Minimal pair
        pair = self.node_list[:2]
        r = self.ecnp(pair)
        print(f"\nk=2 (minimal pair):")
        print(f"  Result: Z={r.get('Z', np.nan):.2f}, rho={r.get('mean_rho', 0):.3f}")
        
        # k=100: Large set
        large = self.node_list[:100]
        r = self.ecnp(large)
        print(f"\nk=100 (large set):")
        if not np.isnan(r.get('Z', np.nan)):
            print(f"  Z={r['Z']:.2f}, pool={r['pool_size']}, rho={r['mean_rho']:.3f}")
        else:
            print(f"  FAILED: {r.get('error')}")
        
        # k=500: Very large
        very_large = self.node_list[:500]
        r = self.ecnp(very_large)
        print(f"\nk=500 (very large set):")
        if not np.isnan(r.get('Z', np.nan)):
            print(f"  Z={r['Z']:.2f}, pool={r['pool_size']}, rho={r['mean_rho']:.3f}")
        else:
            print(f"  FAILED: {r.get('error')}")
    
    # =========================================================================
    # TEST 2: RANDOM TARGET SETS (NULL DISTRIBUTION)
    # =========================================================================
    
    def test_null_distribution(self, n_samples=100, k=10):
        """Random target sets should give Z ~ N(0,1) if properly calibrated."""
        print("\n" + "="*70)
        print(f"TEST 2: NULL DISTRIBUTION (n={n_samples}, k={k})")
        print("="*70)
        
        np.random.seed(42)
        z_values = []
        failures = 0
        
        for _ in range(n_samples):
            targets = list(np.random.choice(self.node_list, k, replace=False))
            r = self.ecnp(targets)
            if not np.isnan(r.get('Z', np.nan)):
                z_values.append(r['Z'])
            else:
                failures += 1
        
        z_values = np.array(z_values)
        
        print(f"\nResults:")
        print(f"  Successful: {len(z_values)}/{n_samples}")
        print(f"  Failures: {failures}")
        print(f"  Mean Z: {z_values.mean():.2f} (expected: 0)")
        print(f"  Std Z: {z_values.std():.2f} (expected: 1)")
        print(f"  Min Z: {z_values.min():.2f}")
        print(f"  Max Z: {z_values.max():.2f}")
        
        # Check if approximately N(0,1)
        mean_ok = abs(z_values.mean()) < 0.5
        std_ok = 0.5 < z_values.std() < 2.0
        
        if mean_ok and std_ok:
            print("  [OK] Null distribution approximately correct")
        else:
            print("  [!] Null distribution may be miscalibrated")
        
        return z_values
    
    # =========================================================================
    # TEST 3: ADVERSARIAL INPUTS
    # =========================================================================
    
    def test_adversarial(self):
        """Test pathological target selections."""
        print("\n" + "="*70)
        print("TEST 3: ADVERSARIAL INPUTS")
        print("="*70)
        
        # Sort by degree and influence
        deg_sorted = sorted(self.degrees.items(), key=lambda x: -x[1])
        inf_sorted = self.m_vector.sort_values(ascending=False)
        
        tests = [
            ("Top 1% degree (hubs)", [g for g, _ in deg_sorted[:int(0.01*len(deg_sorted))]]),
            ("Bottom 10% degree (periphery)", [g for g, _ in deg_sorted[-int(0.1*len(deg_sorted)):]]),
            ("Top 1% influence", inf_sorted.head(int(0.01*len(inf_sorted))).index.tolist()),
            ("Bottom 10% influence", inf_sorted.tail(int(0.1*len(inf_sorted))).index.tolist()),
        ]
        
        for name, targets in tests:
            targets = [t for t in targets if t in self.node_list][:50]  # Cap at 50
            r = self.ecnp(targets)
            print(f"\n{name} (k={len(targets)}):")
            if not np.isnan(r.get('Z', np.nan)):
                print(f"  Z={r['Z']:.2f}, pool={r['pool_size']}, rho={r['mean_rho']:.3f}")
            else:
                print(f"  FAILED: {r.get('error')}, pool_size={r.get('pool_size', 'N/A')}")
    
    # =========================================================================
    # TEST 4: PARAMETER SENSITIVITY
    # =========================================================================
    
    def test_lambda_sensitivity(self, targets_name="Hyperforin"):
        """How sensitive is Z to lambda?"""
        print("\n" + "="*70)
        print("TEST 4: LAMBDA SENSITIVITY")
        print("="*70)
        
        # Use real targets
        targets_df = pd.read_csv(DATA_DIR / "targets_lcc.csv")
        targets = targets_df[targets_df['compound'] == targets_name]['gene_symbol'].tolist()
        
        lambdas = [0.001, 0.005, 0.01, 0.0195, 0.05, 0.1, 0.5, 1.0]
        
        print(f"\n{targets_name} (k={len(targets)}):")
        for lam in lambdas:
            config = ECNPConfig(lambda_redundancy=lam)
            r = self.ecnp(targets, config)
            if not np.isnan(r.get('Z', np.nan)):
                print(f"  lambda={lam:.4f}: Z={r['Z']:.2f}")
            else:
                print(f"  lambda={lam:.4f}: FAILED")
    
    # =========================================================================
    # TEST 5: TOLERANCE SENSITIVITY
    # =========================================================================
    
    def test_tolerance_sensitivity(self):
        """How sensitive is the algorithm to tolerance parameters?"""
        print("\n" + "="*70)
        print("TEST 5: TOLERANCE SENSITIVITY")
        print("="*70)
        
        targets_df = pd.read_csv(DATA_DIR / "targets_lcc.csv")
        targets = targets_df[targets_df['compound'] == 'Hyperforin']['gene_symbol'].tolist()
        
        # Degree tolerance sweep
        print("\nDegree tolerance sweep (percentile=0.1 fixed):")
        for deg_tol in [0.05, 0.1, 0.2, 0.3, 0.5]:
            config = ECNPConfig(degree_tolerance=deg_tol, percentile_window=0.1)
            r = self.ecnp(targets, config)
            if not np.isnan(r.get('Z', np.nan)):
                print(f"  deg_tol={deg_tol:.2f}: Z={r['Z']:.2f}, pool={r['pool_size']}")
            else:
                print(f"  deg_tol={deg_tol:.2f}: FAILED (pool={r.get('pool_size', 'N/A')})")
        
        # Percentile window sweep
        print("\nPercentile window sweep (deg_tol=0.2 fixed):")
        for pct_win in [0.01, 0.05, 0.1, 0.2, 0.5]:
            config = ECNPConfig(degree_tolerance=0.2, percentile_window=pct_win)
            r = self.ecnp(targets, config)
            if not np.isnan(r.get('Z', np.nan)):
                print(f"  pct_win={pct_win:.2f}: Z={r['Z']:.2f}, pool={r['pool_size']}")
            else:
                print(f"  pct_win={pct_win:.2f}: FAILED (pool={r.get('pool_size', 'N/A')})")
    
    # =========================================================================
    # TEST 6: DEGENERATE DISEASE MODULES
    # =========================================================================
    
    def test_degenerate_disease(self):
        """Test with pathological disease definitions."""
        print("\n" + "="*70)
        print("TEST 6: DEGENERATE DISEASE MODULES")
        print("="*70)
        
        targets_df = pd.read_csv(DATA_DIR / "targets_lcc.csv")
        targets = targets_df[targets_df['compound'] == 'Hyperforin']['gene_symbol'].tolist()
        
        # Single gene disease
        single_disease = [self.dili_genes[0]]
        m_single = self.M[[self.node_to_idx[single_disease[0]]], :].sum(axis=0)
        m_single = pd.Series(m_single, index=self.node_list)
        
        # Store original and swap
        orig_m = self.m_vector.copy()
        
        print("\nSingle-gene disease module:")
        # Would need to recompute with different m_vector
        print("  [SKIP] Requires full recomputation")
        
        # Targets overlap with disease
        overlap_targets = self.dili_genes[:10]
        r = self.ecnp(overlap_targets)
        print(f"\nTargets ⊂ Disease (k=10):")
        if not np.isnan(r.get('Z', np.nan)):
            print(f"  Z={r['Z']:.2f}, pool={r['pool_size']}, rho={r['mean_rho']:.3f}")
        else:
            print(f"  FAILED: {r.get('error')}")
    
    def run_all(self):
        """Run complete stress test suite."""
        print("\n" + "="*70)
        print("ECNP STRESS TEST SUITE")
        print("="*70)
        
        self.test_edge_cases()
        self.test_null_distribution(n_samples=100, k=10)
        self.test_adversarial()
        self.test_lambda_sensitivity()
        self.test_tolerance_sensitivity()
        self.test_degenerate_disease()
        
        print("\n" + "="*70)
        print("STRESS TESTING COMPLETE")
        print("="*70)


def main():
    tester = ECNPStressTester()
    tester.run_all()


if __name__ == "__main__":
    main()

"""
Deep Failure Mode Analysis

Investigate the root causes of ECNP failures:
1. Top 1% influence → Z=292 (exploding Z)
2. Percentile window sensitivity (Z=-3 to Z=18)
3. Large k → inflated Z (k=500 → Z=70)

Goal: Understand whether these are fundamental limitations or fixable.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from dataclasses import dataclass
from typing import List, Dict

PROJECT_ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = PROJECT_ROOT / "data" / "processed"
ECNP_DATA_DIR = PROJECT_ROOT / "research" / "ecnp-closed-form" / "data"


@dataclass
class ECNPConfig:
    degree_tolerance: float = 0.2
    percentile_window: float = 0.1
    lambda_redundancy: float = 0.0195


class FailureModeAnalyzer:
    
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
    
    def analyze_top_influence_failure(self):
        """Why does selecting top-influence nodes give Z=292?"""
        print("\n" + "="*70)
        print("FAILURE MODE 1: TOP INFLUENCE TARGETS")
        print("="*70)
        
        # Get top 1% influence nodes
        inf_sorted = self.m_vector.sort_values(ascending=False)
        top_pct = inf_sorted.head(int(0.01 * len(inf_sorted)))
        targets = top_pct.index.tolist()[:50]
        
        print(f"\nTarget set: top {len(targets)} nodes by influence")
        print(f"  Min target influence: {self.m_vector[targets].min():.4f}")
        print(f"  Max target influence: {self.m_vector[targets].max():.4f}")
        print(f"  Mean target influence: {self.m_vector[targets].mean():.4f}")
        
        # What percentile are these?
        target_pcts = [self.pct_ranks[t] for t in targets]
        print(f"  Percentile range: [{min(target_pcts):.2%}, {max(target_pcts):.2%}]")
        
        # Get pool
        config = ECNPConfig()
        pool = self._get_pool(targets, config)
        
        print(f"\nPool (degree+percentile matched):")
        print(f"  Size: {len(pool)}")
        
        if len(pool) > 0:
            pool_m = self.m_vector[pool]
            print(f"  Pool mean influence: {pool_m.mean():.4f}")
            print(f"  Pool std influence: {pool_m.std():.4f}")
            pool_pcts = [self.pct_ranks[p] for p in pool]
            print(f"  Pool percentile range: [{min(pool_pcts):.2%}, {max(pool_pcts):.2%}]")
            
            # The problem
            print(f"\n  GAP:")
            print(f"    Target mean: {self.m_vector[targets].mean():.4f}")
            print(f"    Pool mean: {pool_m.mean():.4f}")
            print(f"    Ratio: {self.m_vector[targets].mean() / pool_m.mean():.1f}x")
            
            # Why the gap?
            print(f"\n  ROOT CAUSE ANALYSIS:")
            print(f"    Targets are at the EXTREME tail (99th+ percentile)")
            print(f"    Pool must match degree AND percentile")
            print(f"    Very few nodes have BOTH high degree AND high percentile")
            print(f"    Pool ends up sampling from lower percentiles")
            print(f"    Result: mu << I(T), Z explodes")
    
    def _get_pool(self, targets, config):
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
                
                deg_ok = t_deg * (1-config.degree_tolerance) <= node_deg <= t_deg * (1+config.degree_tolerance)
                pct_ok = max(0, t_pct - config.percentile_window) <= node_pct <= min(1, t_pct + config.percentile_window)
                
                if deg_ok and pct_ok:
                    matched.add(node)
        return list(matched)
    
    def analyze_percentile_sensitivity(self):
        """Why does percentile window cause Z to swing -3 to +18?"""
        print("\n" + "="*70)
        print("FAILURE MODE 2: PERCENTILE WINDOW SENSITIVITY")
        print("="*70)
        
        targets_df = pd.read_csv(DATA_DIR / "targets_lcc.csv")
        targets = targets_df[targets_df['compound'] == 'Hyperforin']['gene_symbol'].tolist()
        
        print(f"\nHyperforin targets (k={len(targets)})")
        
        for pct_win in [0.01, 0.05, 0.1, 0.2, 0.5]:
            config = ECNPConfig(percentile_window=pct_win)
            pool = self._get_pool(targets, config)
            
            if len(pool) < 5:
                print(f"\npct_win={pct_win}: Pool too small ({len(pool)})")
                continue
            
            pool_m = self.m_vector[pool]
            target_m = self.m_vector[[t for t in targets if t in self.m_vector.index]]
            
            I_T = target_m.sum()
            mu_T = len(target_m) * pool_m.mean()
            
            print(f"\npct_win={pct_win}:")
            print(f"  Pool size: {len(pool)}")
            print(f"  Pool mean: {pool_m.mean():.4f}")
            print(f"  Target mean: {target_m.mean():.4f}")
            print(f"  I(T)={I_T:.4f}, mu_T={mu_T:.4f}")
            print(f"  I(T) - mu = {I_T - mu_T:.4f}")
            
            if pct_win == 0.01:
                print(f"  PROBLEM: Tight window makes pool = target clones")
                print(f"           mu_T approaches I(T), numerator shrinks")
                print(f"           sigma also shrinks, but not as fast")
                print(f"           Result: Z can go negative if pool.mean > target.mean")
    
    def analyze_k_scaling(self):
        """Why does large k give inflated Z?"""
        print("\n" + "="*70)
        print("FAILURE MODE 3: LARGE K INFLATES Z")
        print("="*70)
        
        # Test different k values
        np.random.seed(42)
        
        print("\nZ vs k for random target sets (10 trials each):")
        print(f"{'k':>5} | {'Mean Z':>8} | {'Std Z':>8} | {'Pool Size':>10}")
        print("-" * 40)
        
        for k in [10, 50, 100, 200, 500]:
            z_vals = []
            pool_sizes = []
            
            for _ in range(10):
                targets = list(np.random.choice(self.node_list, k, replace=False))
                r = self._quick_ecnp(targets)
                if r is not None:
                    z_vals.append(r['Z'])
                    pool_sizes.append(r['pool_size'])
            
            if z_vals:
                print(f"{k:>5} | {np.mean(z_vals):>8.2f} | {np.std(z_vals):>8.2f} | {np.mean(pool_sizes):>10.0f}")
        
        print("\n  ANALYSIS:")
        print("    As k grows, random selection increasingly covers entire influence distribution")
        print("    Pool is degree+percentile matched, but union over k targets")
        print("    More targets = larger pool coverage = pool mean regresses to network mean")
        print("    If target set has ANY high-influence nodes, I(T) >> k*pool_mean")
        print("    Result: Z inflates with k")
    
    def _quick_ecnp(self, targets):
        config = ECNPConfig()
        target_idx = [self.node_to_idx[t] for t in targets if t in self.node_to_idx]
        k = len(target_idx)
        if k == 0:
            return None
        
        I_T = sum(self.m_vector.get(t, 0) for t in targets)
        pool = self._get_pool(targets, config)
        if len(pool) < 20:
            return None
        
        pool_m = [self.m_vector.get(p, 0) for p in pool]
        mu_T = k * np.mean(pool_m)
        
        # Redundancy
        M_D = self.M[self.dili_idx, :][:, target_idx]
        norms = np.linalg.norm(M_D, axis=0)
        norms[norms == 0] = 1e-10
        M_D_norm = M_D / norms[np.newaxis, :]
        rho = M_D_norm.T @ M_D_norm
        mean_rho = (np.sum(rho) - np.trace(rho)) / (k * (k - 1)) if k > 1 else 0
        
        sigma_sq = np.var(pool_m, ddof=1) / k + config.lambda_redundancy * mean_rho
        sigma = np.sqrt(sigma_sq) if sigma_sq > 0 else 1e-10
        
        Z = (I_T - mu_T) / sigma
        
        return {'Z': Z, 'pool_size': len(pool)}
    
    def propose_fixes(self):
        """Propose potential fixes for each failure mode."""
        print("\n" + "="*70)
        print("PROPOSED FIXES")
        print("="*70)
        
        print("""
FIX 1: EXTREME INFLUENCE TARGETS
--------------------------------
Problem: Pool cannot match extreme percentiles.
Options:
  a) Cap percentile at 95th (treat 95th-100th as same stratum)
  b) Use log-influence for percentile ranking
  c) Detect and warn when targets are in extreme tail

FIX 2: PERCENTILE WINDOW SENSITIVITY  
-------------------------------------
Problem: 10% window is tuned for Hyperforin/Quercetin, not universal.
Options:
  a) Adaptive window based on target distribution
  b) Minimum pool size constraint: expand window until pool > N
  c) Cross-validation to select window per analysis

FIX 3: K SCALING
----------------
Problem: Large k causes mu underestimation.
Options:
  a) Normalize by k in a different way
  b) Use median instead of mean for pool statistics
  c) Per-target matching (each target gets own pool, aggregate)

RECOMMENDED IMMEDIATE FIX:
--------------------------
Implement pool size floor: if pool < 50, expand percentile window
automatically. This prevents degenerate cases without changing core logic.
""")
    
    def run_analysis(self):
        self.analyze_top_influence_failure()
        self.analyze_percentile_sensitivity()
        self.analyze_k_scaling()
        self.propose_fixes()


def main():
    analyzer = FailureModeAnalyzer()
    analyzer.run_analysis()


if __name__ == "__main__":
    main()

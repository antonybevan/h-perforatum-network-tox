"""
Comprehensive ECNP Stress Test Suite

Real tests, no fluff. 6 categories + trivial limits.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from dataclasses import dataclass
from typing import List, Dict, Tuple
from scipy import sparse

PROJECT_ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = PROJECT_ROOT / "data" / "processed"
ECNP_DATA_DIR = PROJECT_ROOT / "research" / "ecnp-closed-form" / "data"


@dataclass
class Config:
    degree_tolerance: float = 0.2
    percentile_window: float = 0.1
    lambda_redundancy: float = 0.0195


class ComprehensiveStressTester:
    
    def __init__(self):
        self._load_data()
    
    def _load_data(self):
        npz = np.load(ECNP_DATA_DIR / "influence_matrix_900.npz", allow_pickle=True)
        self.M = npz['M']
        self.node_list = npz['node_list'].tolist()
        self.node_to_idx = {g: i for i, g in enumerate(self.node_list)}
        self.n_nodes = len(self.node_list)
        
        dili = pd.read_csv(DATA_DIR / "dili_900_lcc.csv")
        self.dili_genes = dili['gene_name'].tolist()
        self.dili_idx = [self.node_to_idx[g] for g in self.dili_genes if g in self.node_to_idx]
        
        self.m_vector = pd.read_csv(ECNP_DATA_DIR / "dili_influence_vector_900.csv")
        self.m_vector = self.m_vector.set_index('gene')['dili_influence']
        
        edges = pd.read_parquet(DATA_DIR / "network_900_liver_lcc.parquet")
        self.degrees = pd.concat([edges['gene1'], edges['gene2']]).value_counts().to_dict()
        
        sorted_m = self.m_vector.sort_values()
        self.pct_ranks = {g: i / len(sorted_m) for i, g in enumerate(sorted_m.index)}
        self.pct_ranks_array = np.array([self.pct_ranks[g] for g in self.node_list])
    
    def _get_pool(self, targets, deg_tol=0.2, pct_win=0.1):
        matched = set()
        for t in targets:
            if t not in self.m_vector.index:
                continue
            t_deg = self.degrees.get(t, 0)
            t_pct = self.pct_ranks.get(t, 0.5)
            for node in self.m_vector.index:
                if node in targets:
                    continue
                if (t_deg * (1-deg_tol) <= self.degrees.get(node, 0) <= t_deg * (1+deg_tol) and
                    max(0, t_pct - pct_win) <= self.pct_ranks.get(node, 0) <= min(1, t_pct + pct_win)):
                    matched.add(node)
        return list(matched)
    
    def _compute_rho(self, target_idx):
        k = len(target_idx)
        if k < 2:
            return 0.0
        M_D = self.M[self.dili_idx, :][:, target_idx]
        norms = np.linalg.norm(M_D, axis=0)
        norms[norms == 0] = 1e-10
        M_D_norm = M_D / norms[np.newaxis, :]
        rho = M_D_norm.T @ M_D_norm
        return (np.sum(rho) - np.trace(rho)) / (k * (k - 1))
    
    def _ecnp(self, targets, config=None):
        if config is None:
            config = Config()
        target_idx = [self.node_to_idx[t] for t in targets if t in self.node_to_idx]
        k = len(target_idx)
        if k == 0:
            return {'Z': np.nan}
        
        I_T = sum(self.m_vector.get(t, 0) for t in targets)
        pool = self._get_pool(targets, config.degree_tolerance, config.percentile_window)
        if len(pool) < 20:
            return {'Z': np.nan, 'error': 'pool_small'}
        
        pool_m = [self.m_vector.get(p, 0) for p in pool]
        mu_T = k * np.mean(pool_m)
        mean_rho = self._compute_rho(target_idx)
        sigma_sq = np.var(pool_m, ddof=1) / k + config.lambda_redundancy * mean_rho
        sigma = np.sqrt(sigma_sq) if sigma_sq > 0 else 1e-10
        Z = (I_T - mu_T) / sigma
        
        return {'Z': Z, 'k': k, 'mu': mu_T, 'sigma': sigma, 'rho': mean_rho, 'pool': len(pool)}
    
    # =========================================================================
    # TEST 1: ADVERSARIAL TARGET SETS
    # =========================================================================
    
    def test_adversarial_targets(self):
        print("\n" + "="*70)
        print("TEST 1: ADVERSARIAL TARGET SETS")
        print("="*70)
        
        np.random.seed(42)
        
        # 1a: Random targets from same influence percentile as Hyperforin
        # Hyperforin is around 95-99th percentile
        print("\n1a. Random targets matching Hyperforin percentile (but not drug targets):")
        high_pct_nodes = [g for g, p in self.pct_ranks.items() if 0.9 < p < 0.99]
        fake_hyp = list(np.random.choice(high_pct_nodes, 10, replace=False))
        r = self._ecnp(fake_hyp)
        print(f"  Z = {r['Z']:.2f} (expect ~0-2)")
        
        # 1b: High degree but low influence
        print("\n1b. High-degree but low-influence nodes:")
        deg_sorted = sorted(self.degrees.items(), key=lambda x: -x[1])
        high_deg = [g for g, _ in deg_sorted[:200]]
        low_inf_high_deg = [g for g in high_deg if self.pct_ranks.get(g, 0) < 0.3][:10]
        r = self._ecnp(low_inf_high_deg)
        print(f"  Z = {r['Z']:.2f} (expect ~0)")
        
        # 1c: Maximally redundant (pick genes very close in the network)
        print("\n1c. Maximally redundant set (similar pathway genes):")
        # Use genes from a common pathway - pick top of one column of M
        col_idx = 100  # arbitrary node
        col_influences = self.M[self.dili_idx, :].sum(axis=0)
        correlations = []
        for i in range(self.n_nodes):
            if i == col_idx:
                continue
            corr = np.corrcoef(self.M[:, col_idx], self.M[:, i])[0, 1]
            correlations.append((i, corr))
        most_similar = sorted(correlations, key=lambda x: -x[1])[:10]
        redundant_set = [self.node_list[col_idx]] + [self.node_list[i] for i, _ in most_similar[:9]]
        r = self._ecnp(redundant_set)
        print(f"  Z = {r['Z']:.2f}, rho = {r.get('rho', 0):.3f} (expect high rho)")
        
        # 1d: Maximally orthogonal (pick genes with minimal correlation)
        print("\n1d. Maximally orthogonal set:")
        least_similar = sorted(correlations, key=lambda x: x[1])[:10]
        orthogonal_set = [self.node_list[col_idx]] + [self.node_list[i] for i, _ in least_similar[:9]]
        r = self._ecnp(orthogonal_set)
        print(f"  Z = {r['Z']:.2f}, rho = {r.get('rho', 0):.3f} (expect low rho)")
    
    # =========================================================================
    # TEST 2: DISEASE MODULE CORRUPTION
    # =========================================================================
    
    def test_disease_corruption(self):
        print("\n" + "="*70)
        print("TEST 2: DISEASE MODULE CORRUPTION")
        print("="*70)
        
        np.random.seed(42)
        
        targets_df = pd.read_csv(DATA_DIR / "targets_lcc.csv")
        hyp = targets_df[targets_df['compound'] == 'Hyperforin']['gene_symbol'].tolist()
        
        # Original DILI
        r_orig = self._ecnp(hyp)
        print(f"\nOriginal DILI (82 genes): Z = {r_orig['Z']:.2f}")
        
        # Gradual corruption - full recomputation for each
        print("\nGradual dilution of disease module:")
        for corruption in [0, 0.25, 0.5, 0.75, 1.0]:
            n_corrupt = int(len(self.dili_genes) * corruption)
            n_keep = len(self.dili_genes) - n_corrupt
            
            # Keep first n_keep real DILI genes, replace rest with random
            kept = self.dili_genes[:n_keep]
            non_dili = [g for g in self.node_list if g not in self.dili_genes]
            random_genes = list(np.random.choice(non_dili, n_corrupt, replace=False))
            corrupted_module = kept + random_genes
            
            # Recompute m_vector for corrupted module
            corrupt_idx = [self.node_to_idx[g] for g in corrupted_module if g in self.node_to_idx]
            m_corrupt = self.M[corrupt_idx, :].sum(axis=0)
            m_corrupt_series = pd.Series(m_corrupt, index=self.node_list)
            
            # Recompute percentile ranks with corrupted m_vector
            sorted_m = m_corrupt_series.sort_values()
            pct_corrupt = {g: i / len(sorted_m) for i, g in enumerate(sorted_m.index)}
            
            # Build pool with corrupted percentiles
            matched = set()
            for t in hyp:
                if t not in m_corrupt_series.index:
                    continue
                t_deg = self.degrees.get(t, 0)
                t_pct = pct_corrupt.get(t, 0.5)
                for node in m_corrupt_series.index:
                    if node in hyp:
                        continue
                    if (t_deg * 0.8 <= self.degrees.get(node, 0) <= t_deg * 1.2 and
                        max(0, t_pct - 0.1) <= pct_corrupt.get(node, 0) <= min(1, t_pct + 0.1)):
                        matched.add(node)
            
            pool = list(matched)
            if len(pool) < 20:
                print(f"  {corruption*100:.0f}% corrupt: pool too small ({len(pool)})")
                continue
            
            # Compute ECNP with corrupted system
            I_T = sum(m_corrupt_series.get(t, 0) for t in hyp)
            pool_m = [m_corrupt_series.get(p, 0) for p in pool]
            mu_T = len(hyp) * np.mean(pool_m)
            sigma = np.std(pool_m, ddof=1) / np.sqrt(len(hyp))
            Z = (I_T - mu_T) / sigma if sigma > 0 else np.nan
            
            print(f"  {corruption*100:.0f}% corrupt: Z = {Z:.2f} (pool={len(pool)})")
    
    # =========================================================================
    # TEST 3: INFLUENCE RANK SCRAMBLING
    # =========================================================================
    
    def test_rank_scrambling(self):
        print("\n" + "="*70)
        print("TEST 3: INFLUENCE RANK SCRAMBLING")
        print("="*70)
        print("\nTests whether percentile conditioning is doing real work.")
        print("If scrambling ranks doesn't change Z, the conditioning is placebo.")
        
        targets_df = pd.read_csv(DATA_DIR / "targets_lcc.csv")
        hyp = targets_df[targets_df['compound'] == 'Hyperforin']['gene_symbol'].tolist()
        
        # Original
        r_orig = self._ecnp(hyp)
        print(f"\nOriginal (real ranks): Z = {r_orig['Z']:.2f}")
        
        # Run 10 scrambles and collect Z values
        print("\nWith 10 different rank scrambles:")
        z_scrambled = []
        
        for seed in range(10):
            np.random.seed(seed)
            scrambled_pct = self.pct_ranks_array.copy()
            np.random.shuffle(scrambled_pct)
            scrambled_ranks = {self.node_list[i]: p for i, p in enumerate(scrambled_pct)}
            
            # Build pool with scrambled ranks
            matched = set()
            for t in hyp:
                if t not in self.m_vector.index:
                    continue
                t_deg = self.degrees.get(t, 0)
                t_pct = scrambled_ranks.get(t, 0.5)
                for node in self.m_vector.index:
                    if node in hyp:
                        continue
                    node_pct = scrambled_ranks.get(node, 0)
                    if (t_deg * 0.8 <= self.degrees.get(node, 0) <= t_deg * 1.2 and
                        max(0, t_pct - 0.1) <= node_pct <= min(1, t_pct + 0.1)):
                        matched.add(node)
            
            pool = list(matched)
            if len(pool) < 20:
                continue
            
            I_T = sum(self.m_vector.get(t, 0) for t in hyp)
            pool_m = [self.m_vector.get(p, 0) for p in pool]
            mu_T = len(hyp) * np.mean(pool_m)
            sigma = np.std(pool_m, ddof=1) / np.sqrt(len(hyp))
            Z = (I_T - mu_T) / sigma if sigma > 0 else np.nan
            z_scrambled.append(Z)
        
        if z_scrambled:
            print(f"  Z values: {[f'{z:.1f}' for z in z_scrambled]}")
            print(f"  Mean: {np.mean(z_scrambled):.1f}, Std: {np.std(z_scrambled):.1f}")
            print(f"  Range: [{min(z_scrambled):.1f}, {max(z_scrambled):.1f}]")
            print(f"\n  Original Z = {r_orig['Z']:.2f}")
            print(f"  If scrambled Z >> original, percentile matching is CRITICAL")
            print(f"  If scrambled Z ~ original, percentile matching doesn't matter")
    
    # =========================================================================
    # TEST 4: ASYMPTOTIC K BEHAVIOR
    # =========================================================================
    
    def test_k_asymptotic(self):
        print("\n" + "="*70)
        print("TEST 4: ASYMPTOTIC K BEHAVIOR")
        print("="*70)
        
        np.random.seed(42)
        
        # Pick nodes from a fixed percentile band (50-70th)
        mid_pct_nodes = [g for g, p in self.pct_ranks.items() if 0.5 < p < 0.7]
        
        print("\nk vs Z (targets from 50-70th percentile):")
        print(f"{'k':>5} | {'mu':>10} | {'sigma':>10} | {'Z':>8} | {'rho':>6}")
        print("-" * 50)
        
        for k in [5, 10, 20, 50, 100, 200]:
            if k > len(mid_pct_nodes):
                break
            targets = list(np.random.choice(mid_pct_nodes, k, replace=False))
            r = self._ecnp(targets)
            if r.get('Z') is not None and not np.isnan(r.get('Z', np.nan)):
                print(f"{k:>5} | {r['mu']:>10.4f} | {r['sigma']:>10.4f} | {r['Z']:>8.2f} | {r.get('rho', 0):>6.3f}")
            else:
                print(f"{k:>5} | FAILED")
        
        print("\nExpect: mu ~ k, sigma ~ sqrt(k_eff), Z should stabilize")
    
    # =========================================================================
    # TEST 5: TRIVIAL LIMITS
    # =========================================================================
    
    def test_trivial_limits(self):
        print("\n" + "="*70)
        print("TEST 5: TRIVIAL LIMIT SANITY CHECKS")
        print("="*70)
        
        # 5a: Single hub node (high degree, high influence)
        print("\n5a. Single hub node (highest influence):")
        top_inf = self.m_vector.sort_values(ascending=False).head(1).index[0]
        r = self._ecnp([top_inf, top_inf])  # Duplicate to meet k>=2
        print(f"  {top_inf}: influence={self.m_vector[top_inf]:.4f}, percentile={self.pct_ranks[top_inf]:.1%}")
        
        # 5b: Single peripheral node (low degree, low influence)
        print("\n5b. Single peripheral node (lowest influence):")
        bot_inf = self.m_vector.sort_values(ascending=True).head(1).index[0]
        r = self._ecnp([bot_inf, bot_inf])
        print(f"  {bot_inf}: influence={self.m_vector[bot_inf]:.4f}, percentile={self.pct_ranks[bot_inf]:.1%}")
        
        # 5c: Two identical targets (check if rho = 1)
        print("\n5c. Two 'identical' targets (same node twice - dummy check):")
        print("  [Skipped - need distinct nodes for rho calculation]")
        
        # 5d: Two very similar targets
        print("\n5d. Two highly correlated targets:")
        col_idx = 0
        correlations = []
        for i in range(1, min(100, self.n_nodes)):
            corr = np.corrcoef(self.M[:, col_idx], self.M[:, i])[0, 1]
            correlations.append((i, corr))
        most_similar = sorted(correlations, key=lambda x: -x[1])[0]
        t1, t2 = self.node_list[col_idx], self.node_list[most_similar[0]]
        r = self._ecnp([t1, t2])
        print(f"  {t1}, {t2}: rho = {r.get('rho', 0):.3f}, Z = {r.get('Z', np.nan):.2f}")
        
        # 5e: Two orthogonal targets
        print("\n5e. Two uncorrelated targets:")
        least_corr = sorted(correlations, key=lambda x: abs(x[1]))[0]
        t1, t2 = self.node_list[col_idx], self.node_list[least_corr[0]]
        r = self._ecnp([t1, t2])
        print(f"  {t1}, {t2}: rho = {r.get('rho', 0):.3f}, Z = {r.get('Z', np.nan):.2f}")
    
    def run_all(self):
        print("="*70)
        print("COMPREHENSIVE ECNP STRESS TEST SUITE")
        print("="*70)
        
        self.test_adversarial_targets()
        self.test_disease_corruption()
        self.test_rank_scrambling()
        self.test_k_asymptotic()
        self.test_trivial_limits()
        
        print("\n" + "="*70)
        print("STRESS TESTING COMPLETE")
        print("="*70)


if __name__ == "__main__":
    tester = ComprehensiveStressTester()
    tester.run_all()

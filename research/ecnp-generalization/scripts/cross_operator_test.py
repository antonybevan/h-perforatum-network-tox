"""
Cross-Operator Sanity Test

Tests whether ECNP generalizes across different propagation operators.

Tests:
1. Different RWR restart probabilities α (0.05, 0.10, 0.15, 0.20, 0.30)
2. Validate that percentile conditioning is still required

Expectations:
- Absolute Z values change with α
- Null logic (percentile matching) still holds
- Relative ordering of compounds persists
"""

import numpy as np
import pandas as pd
from pathlib import Path
from scipy import sparse
from typing import List, Dict

PROJECT_ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = PROJECT_ROOT / "data" / "processed"


class CrossOperatorTester:
    
    def __init__(self):
        self._load_data()
    
    def _load_data(self):
        # Load network
        self.edges_df = pd.read_parquet(DATA_DIR / "network_900_liver_lcc.parquet")
        all_nodes = set(self.edges_df['gene1'].tolist() + self.edges_df['gene2'].tolist())
        self.node_list = sorted(all_nodes)
        self.node_to_idx = {g: i for i, g in enumerate(self.node_list)}
        self.n_nodes = len(self.node_list)
        
        # Load DILI genes
        dili = pd.read_csv(DATA_DIR / "dili_900_lcc.csv")
        self.dili_genes = dili['gene_name'].tolist()
        self.dili_idx = [self.node_to_idx[g] for g in self.dili_genes if g in self.node_to_idx]
        
        # Load targets
        targets_df = pd.read_csv(DATA_DIR / "targets_lcc.csv")
        self.hyperforin = targets_df[targets_df['compound'] == 'Hyperforin']['gene_symbol'].tolist()
        self.quercetin = targets_df[targets_df['compound'] == 'Quercetin']['gene_symbol'].tolist()
        
        print(f"Loaded: {len(self.edges_df)} edges, {self.n_nodes} nodes, {len(self.dili_genes)} DILI genes")
    
    def build_transition_matrix(self) -> sparse.csr_matrix:
        """Build column-stochastic transition matrix."""
        n = self.n_nodes
        rows, cols, data = [], [], []
        
        for _, row in self.edges_df.iterrows():
            i = self.node_to_idx.get(row['gene1'])
            j = self.node_to_idx.get(row['gene2'])
            if i is not None and j is not None:
                rows.extend([i, j])
                cols.extend([j, i])
                data.extend([1, 1])
        
        A = sparse.csr_matrix((data, (rows, cols)), shape=(n, n))
        col_sums = np.array(A.sum(axis=0)).flatten()
        col_sums[col_sums == 0] = 1
        W = A.multiply(1 / col_sums)
        
        return W.tocsr()
    
    def compute_influence_matrix(self, W: sparse.csr_matrix, alpha: float) -> np.ndarray:
        """Compute M = alpha * (I - (1-alpha)*W)^-1 for given alpha."""
        print(f"  Computing M with alpha={alpha}...", end=" ")
        n = W.shape[0]
        I = sparse.eye(n, format='csr')
        A = I - (1 - alpha) * W
        
        A_dense = A.toarray()
        A_inv = np.linalg.inv(A_dense)
        M = alpha * A_inv
        print("done")
        
        return M
    
    def compute_ecnp(self, M: np.ndarray, targets: List[str], dili_idx: List[int],
                     use_percentile: bool = True) -> Dict:
        """Compute ECNP with or without percentile matching."""
        # Compute m_vector
        m_vector = M[dili_idx, :].sum(axis=0)
        m_series = pd.Series(m_vector, index=self.node_list)
        
        # Get target values
        target_idx = [self.node_to_idx[t] for t in targets if t in self.node_to_idx]
        k = len(target_idx)
        if k == 0:
            return {'Z': np.nan}
        
        I_T = sum(m_series.get(t, 0) for t in targets)
        
        if use_percentile:
            # Percentile-matched pool
            sorted_m = m_series.sort_values()
            pct_ranks = {g: i / len(sorted_m) for i, g in enumerate(sorted_m.index)}
            
            matched = set()
            for t in targets:
                if t not in m_series.index:
                    continue
                t_pct = pct_ranks.get(t, 0.5)
                for node in m_series.index:
                    if node in targets:
                        continue
                    node_pct = pct_ranks.get(node, 0)
                    if max(0, t_pct - 0.1) <= node_pct <= min(1, t_pct + 0.1):
                        matched.add(node)
            
            pool = list(matched)
        else:
            # Random pool (no percentile matching)
            non_targets = [g for g in self.node_list if g not in targets]
            pool = list(np.random.choice(non_targets, min(1000, len(non_targets)), replace=False))
        
        if len(pool) < 20:
            return {'Z': np.nan, 'error': 'pool_small'}
        
        pool_m = [m_series.get(p, 0) for p in pool]
        mu_T = k * np.mean(pool_m)
        sigma = np.std(pool_m, ddof=1) / np.sqrt(k)
        Z = (I_T - mu_T) / sigma if sigma > 0 else np.nan
        
        return {'Z': Z, 'I_T': I_T, 'mu_T': mu_T, 'pool': len(pool)}
    
    def test_alpha_variation(self):
        """Test different RWR restart probabilities."""
        print("\n" + "="*70)
        print("TEST 1: RWR RESTART PROBABILITY (alpha) VARIATION")
        print("="*70)
        
        W = self.build_transition_matrix()
        alphas = [0.05, 0.10, 0.15, 0.20, 0.30]
        
        results = []
        
        for alpha in alphas:
            print(f"\nalpha = {alpha}")
            M = self.compute_influence_matrix(W, alpha)
            
            # Run with percentile matching
            hyp_pct = self.compute_ecnp(M, self.hyperforin, self.dili_idx, use_percentile=True)
            que_pct = self.compute_ecnp(M, self.quercetin, self.dili_idx, use_percentile=True)
            
            # Run without percentile matching (random pool)
            np.random.seed(42)
            hyp_rand = self.compute_ecnp(M, self.hyperforin, self.dili_idx, use_percentile=False)
            que_rand = self.compute_ecnp(M, self.quercetin, self.dili_idx, use_percentile=False)
            
            print(f"  Hyperforin: Z_pct={hyp_pct['Z']:.2f}, Z_rand={hyp_rand['Z']:.2f}")
            print(f"  Quercetin:  Z_pct={que_pct['Z']:.2f}, Z_rand={que_rand['Z']:.2f}")
            
            results.append({
                'alpha': alpha,
                'hyperforin_pct': hyp_pct['Z'],
                'quercetin_pct': que_pct['Z'],
                'hyperforin_rand': hyp_rand['Z'],
                'quercetin_rand': que_rand['Z']
            })
        
        return pd.DataFrame(results)
    
    def test_percentile_necessity(self):
        """Test whether percentile matching is still needed across operators."""
        print("\n" + "="*70)
        print("TEST 2: PERCENTILE MATCHING NECESSITY ACROSS alpha")
        print("="*70)
        
        W = self.build_transition_matrix()
        alphas = [0.10, 0.15, 0.20]
        
        print("\nFor each alpha, compare Z with/without percentile matching:")
        print(f"{'alpha':>6} | {'Compound':>12} | {'Z_pct':>8} | {'Z_rand':>8} | {'Ratio':>6}")
        print("-" * 60)
        
        for alpha in alphas:
            M = self.compute_influence_matrix(W, alpha)
            
            for compound, targets in [("Hyperforin", self.hyperforin), ("Quercetin", self.quercetin)]:
                pct_result = self.compute_ecnp(M, targets, self.dili_idx, use_percentile=True)
                
                np.random.seed(42)
                rand_result = self.compute_ecnp(M, targets, self.dili_idx, use_percentile=False)
                
                if not np.isnan(pct_result['Z']) and not np.isnan(rand_result['Z']):
                    ratio = rand_result['Z'] / pct_result['Z'] if pct_result['Z'] != 0 else np.inf
                    print(f"{alpha:>6.2f} | {compound:>12} | {pct_result['Z']:>8.2f} | {rand_result['Z']:>8.2f} | {ratio:>6.2f}x")
    
    def run_all(self):
        print("="*70)
        print("CROSS-OPERATOR SANITY TESTS")
        print("="*70)
        
        results_df = self.test_alpha_variation()
        self.test_percentile_necessity()
        
        print("\n" + "="*70)
        print("SUMMARY")
        print("="*70)
        
        print("\nZ-scores across alpha (with percentile matching):")
        print(results_df[['alpha', 'hyperforin_pct', 'quercetin_pct']].to_string(index=False))
        
        print("\n" + "="*70)
        print("CROSS-OPERATOR TESTS COMPLETE")
        print("="*70)
        
        return results_df


if __name__ == "__main__":
    tester = CrossOperatorTester()
    results = tester.run_all()

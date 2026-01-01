"""
Network Perturbation Stress Test

Tests ECNP robustness to interactome uncertainty.

Tests:
1. Edge dropping (5%, 10%, 15%)
2. Edge rewiring (degree-preserving)
3. Measure Z stability across perturbations

Expectations:
- Absolute Z changes
- Relative disease alignment persists
- λ remains in same order of magnitude
"""

import numpy as np
import pandas as pd
from pathlib import Path
from scipy import sparse
from typing import List, Dict, Tuple

PROJECT_ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = PROJECT_ROOT / "data" / "processed"
ECNP_DATA_DIR = PROJECT_ROOT / "research" / "ecnp-closed-form" / "data"

ALPHA = 0.15  # RWR restart probability


class NetworkPerturbationTester:
    
    def __init__(self):
        self._load_data()
    
    def _load_data(self):
        # Load original network
        self.edges_df = pd.read_parquet(DATA_DIR / "network_900_liver_lcc.parquet")
        print(f"Loaded network: {len(self.edges_df)} edges")
        
        # Build node list
        all_nodes = set(self.edges_df['gene1'].tolist() + self.edges_df['gene2'].tolist())
        self.node_list = sorted(all_nodes)
        self.node_to_idx = {g: i for i, g in enumerate(self.node_list)}
        self.n_nodes = len(self.node_list)
        print(f"Nodes: {self.n_nodes}")
        
        # Load DILI genes
        dili = pd.read_csv(DATA_DIR / "dili_900_lcc.csv")
        self.dili_genes = dili['gene_name'].tolist()
        self.dili_idx = [self.node_to_idx[g] for g in self.dili_genes if g in self.node_to_idx]
        
        # Load targets
        targets_df = pd.read_csv(DATA_DIR / "targets_lcc.csv")
        self.hyperforin = targets_df[targets_df['compound'] == 'Hyperforin']['gene_symbol'].tolist()
        self.quercetin = targets_df[targets_df['compound'] == 'Quercetin']['gene_symbol'].tolist()
    
    def build_transition_matrix(self, edges: pd.DataFrame) -> sparse.csr_matrix:
        """Build column-stochastic transition matrix from edge list."""
        n = self.n_nodes
        
        # Build adjacency
        rows, cols, data = [], [], []
        for _, row in edges.iterrows():
            i = self.node_to_idx.get(row['gene1'])
            j = self.node_to_idx.get(row['gene2'])
            if i is not None and j is not None:
                rows.extend([i, j])
                cols.extend([j, i])
                data.extend([1, 1])
        
        A = sparse.csr_matrix((data, (rows, cols)), shape=(n, n))
        
        # Normalize to column-stochastic
        col_sums = np.array(A.sum(axis=0)).flatten()
        col_sums[col_sums == 0] = 1
        W = A.multiply(1 / col_sums)
        
        return W.tocsr()
    
    def compute_influence_matrix(self, W: sparse.csr_matrix) -> np.ndarray:
        """Compute M = alpha * (I - (1-alpha)*W)^-1"""
        n = W.shape[0]
        I = sparse.eye(n, format='csr')
        A = I - (1 - ALPHA) * W
        
        # Dense inversion for moderate n
        A_dense = A.toarray()
        A_inv = np.linalg.inv(A_dense)
        M = ALPHA * A_inv
        
        return M
    
    def compute_ecnp(self, M: np.ndarray, targets: List[str], dili_idx: List[int]) -> Dict:
        """Compute ECNP with given M."""
        # Compute m_vector
        m_vector = M[dili_idx, :].sum(axis=0)
        m_series = pd.Series(m_vector, index=self.node_list)
        
        # Percentile ranks
        sorted_m = m_series.sort_values()
        pct_ranks = {g: i / len(sorted_m) for i, g in enumerate(sorted_m.index)}
        
        # Degrees from M (diagonal gives self-influence, approximate degree)
        degrees = {}
        for g in self.node_list:
            idx = self.node_to_idx[g]
            degrees[g] = np.sum(M[idx, :] > 0.001)  # Count non-trivial influences
        
        # Get pool
        target_idx = [self.node_to_idx[t] for t in targets if t in self.node_to_idx]
        k = len(target_idx)
        if k == 0:
            return {'Z': np.nan}
        
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
        if len(pool) < 20:
            return {'Z': np.nan, 'error': 'pool_small'}
        
        I_T = sum(m_series.get(t, 0) for t in targets)
        pool_m = [m_series.get(p, 0) for p in pool]
        mu_T = k * np.mean(pool_m)
        sigma = np.std(pool_m, ddof=1) / np.sqrt(k)
        Z = (I_T - mu_T) / sigma if sigma > 0 else np.nan
        
        return {'Z': Z, 'I_T': I_T, 'mu_T': mu_T, 'pool': len(pool)}
    
    def drop_edges(self, fraction: float) -> pd.DataFrame:
        """Randomly drop a fraction of edges."""
        n_drop = int(len(self.edges_df) * fraction)
        drop_idx = np.random.choice(len(self.edges_df), n_drop, replace=False)
        return self.edges_df.drop(self.edges_df.index[drop_idx])
    
    def rewire_edges(self, n_rewires: int) -> pd.DataFrame:
        """
        Rewire edges while preserving degree distribution.
        Uses edge-swap method.
        """
        edges = self.edges_df.copy()
        edge_list = list(zip(edges['gene1'], edges['gene2']))
        n_edges = len(edge_list)
        
        for _ in range(n_rewires):
            # Pick two random edges
            i, j = np.random.choice(n_edges, 2, replace=False)
            e1, e2 = edge_list[i], edge_list[j]
            
            # Swap endpoints: (a,b), (c,d) -> (a,d), (c,b)
            a, b = e1
            c, d = e2
            
            # Avoid self-loops and duplicates
            if a != d and c != b:
                if (a, d) not in edge_list and (d, a) not in edge_list:
                    if (c, b) not in edge_list and (b, c) not in edge_list:
                        edge_list[i] = (a, d)
                        edge_list[j] = (c, b)
        
        return pd.DataFrame({'gene1': [e[0] for e in edge_list], 
                            'gene2': [e[1] for e in edge_list]})
    
    def test_edge_dropping(self):
        """Test robustness to edge removal."""
        print("\n" + "="*70)
        print("TEST 1: EDGE DROPPING")
        print("="*70)
        
        np.random.seed(42)
        
        # Original
        print("\nComputing original ECNP...")
        W_orig = self.build_transition_matrix(self.edges_df)
        M_orig = self.compute_influence_matrix(W_orig)
        
        z_hyp_orig = self.compute_ecnp(M_orig, self.hyperforin, self.dili_idx)['Z']
        z_que_orig = self.compute_ecnp(M_orig, self.quercetin, self.dili_idx)['Z']
        
        print(f"Original: Hyperforin Z={z_hyp_orig:.2f}, Quercetin Z={z_que_orig:.2f}")
        
        # Test different drop fractions
        print("\nWith edge dropping:")
        for drop_frac in [0.05, 0.10, 0.15]:
            print(f"\n{drop_frac*100:.0f}% edges dropped:")
            
            # Run 3 trials
            z_hyp_trials = []
            z_que_trials = []
            
            for trial in range(3):
                np.random.seed(42 + trial)
                edges_dropped = self.drop_edges(drop_frac)
                W = self.build_transition_matrix(edges_dropped)
                M = self.compute_influence_matrix(W)
                
                z_hyp = self.compute_ecnp(M, self.hyperforin, self.dili_idx)['Z']
                z_que = self.compute_ecnp(M, self.quercetin, self.dili_idx)['Z']
                
                if not np.isnan(z_hyp):
                    z_hyp_trials.append(z_hyp)
                if not np.isnan(z_que):
                    z_que_trials.append(z_que)
            
            if z_hyp_trials:
                hyp_mean = np.mean(z_hyp_trials)
                hyp_change = (hyp_mean - z_hyp_orig) / z_hyp_orig * 100
                print(f"  Hyperforin: Z={hyp_mean:.2f} ({hyp_change:+.1f}% from original)")
            
            if z_que_trials:
                que_mean = np.mean(z_que_trials)
                que_change = (que_mean - z_que_orig) / z_que_orig * 100
                print(f"  Quercetin: Z={que_mean:.2f} ({que_change:+.1f}% from original)")
    
    def test_edge_rewiring(self):
        """Test robustness to edge rewiring (degree-preserving)."""
        print("\n" + "="*70)
        print("TEST 2: EDGE REWIRING (degree-preserving)")
        print("="*70)
        
        np.random.seed(42)
        
        # Original
        W_orig = self.build_transition_matrix(self.edges_df)
        M_orig = self.compute_influence_matrix(W_orig)
        
        z_hyp_orig = self.compute_ecnp(M_orig, self.hyperforin, self.dili_idx)['Z']
        z_que_orig = self.compute_ecnp(M_orig, self.quercetin, self.dili_idx)['Z']
        
        print(f"\nOriginal: Hyperforin Z={z_hyp_orig:.2f}, Quercetin Z={z_que_orig:.2f}")
        
        # Test different rewire counts
        n_edges = len(self.edges_df)
        print("\nWith edge rewiring:")
        
        for rewire_frac in [0.05, 0.10, 0.20]:
            n_rewires = int(n_edges * rewire_frac)
            print(f"\n{rewire_frac*100:.0f}% edges rewired ({n_rewires} swaps):")
            
            np.random.seed(42)
            edges_rewired = self.rewire_edges(n_rewires)
            W = self.build_transition_matrix(edges_rewired)
            M = self.compute_influence_matrix(W)
            
            z_hyp = self.compute_ecnp(M, self.hyperforin, self.dili_idx)['Z']
            z_que = self.compute_ecnp(M, self.quercetin, self.dili_idx)['Z']
            
            if not np.isnan(z_hyp):
                hyp_change = (z_hyp - z_hyp_orig) / z_hyp_orig * 100
                print(f"  Hyperforin: Z={z_hyp:.2f} ({hyp_change:+.1f}% from original)")
            
            if not np.isnan(z_que):
                que_change = (z_que - z_que_orig) / z_que_orig * 100
                print(f"  Quercetin: Z={z_que:.2f} ({que_change:+.1f}% from original)")
    
    def run_all(self):
        print("="*70)
        print("NETWORK PERTURBATION STRESS TESTS")
        print("="*70)
        
        self.test_edge_dropping()
        self.test_edge_rewiring()
        
        print("\n" + "="*70)
        print("NETWORK PERTURBATION TESTS COMPLETE")
        print("="*70)


if __name__ == "__main__":
    tester = NetworkPerturbationTester()
    tester.run_all()

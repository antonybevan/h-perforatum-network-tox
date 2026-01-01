"""
Cross-Disease ECNP Validation

Tests whether the ECNP closed-form algorithm with λ=0.0195 generalizes
across different disease modules.

For each disease:
1. Compute disease-specific m_j vector
2. Run ECNP on Hyperforin and Quercetin
3. Compare to DILI baseline

Key question: Does λ=0.0195 work universally, or is it disease-specific?
"""

import numpy as np
import pandas as pd
from pathlib import Path
from dataclasses import dataclass
from typing import List, Dict
import sys

# Paths
PROJECT_ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = PROJECT_ROOT / "data" / "processed"
ECNP_DATA_DIR = PROJECT_ROOT / "research" / "ecnp-closed-form" / "data"
DISEASE_DIR = PROJECT_ROOT / "research" / "ecnp-generalization" / "data"
RESULTS_DIR = PROJECT_ROOT / "research" / "ecnp-generalization" / "results"


@dataclass
class ECNPConfig:
    degree_tolerance: float = 0.2
    percentile_window: float = 0.1
    lambda_redundancy: float = 0.0195


class CrossDiseaseECNP:
    """ECNP analysis across multiple disease modules."""
    
    def __init__(self):
        self._load_base_data()
    
    def _load_base_data(self):
        """Load network and influence matrix (shared across diseases)."""
        # Influence matrix
        npz = np.load(ECNP_DATA_DIR / "influence_matrix_900.npz", allow_pickle=True)
        self.M = npz['M']
        self.node_list = npz['node_list'].tolist()
        self.node_to_idx = {g: i for i, g in enumerate(self.node_list)}
        
        # Network degrees
        edges = pd.read_parquet(DATA_DIR / "network_900_liver_lcc.parquet")
        self.degrees = pd.concat([edges['gene1'], edges['gene2']]).value_counts().to_dict()
        
        # Targets
        targets_df = pd.read_csv(DATA_DIR / "targets_lcc.csv")
        self.hyperforin = targets_df[targets_df['compound'] == 'Hyperforin']['gene_symbol'].tolist()
        self.quercetin = targets_df[targets_df['compound'] == 'Quercetin']['gene_symbol'].tolist()
    
    def compute_disease_influence(self, disease_genes: List[str]) -> pd.Series:
        """Compute per-node influence vector for a disease module."""
        disease_idx = [self.node_to_idx[g] for g in disease_genes if g in self.node_to_idx]
        
        # m_j = sum_{i in D} M[i, j]
        m_vector = self.M[disease_idx, :].sum(axis=0)
        
        return pd.Series(m_vector, index=self.node_list)
    
    def get_percentile_ranks(self, m_vector: pd.Series) -> Dict[str, float]:
        """Compute percentile ranks for a given m_vector."""
        sorted_m = m_vector.sort_values()
        return {g: i / len(sorted_m) for i, g in enumerate(sorted_m.index)}
    
    def get_stratum_pool(self, targets: List[str], m_vector: pd.Series, 
                          config: ECNPConfig) -> List[str]:
        """Get percentile-rank matched pool."""
        pct_ranks = self.get_percentile_ranks(m_vector)
        matched = set()
        
        for t in targets:
            if t not in m_vector.index:
                continue
            
            t_deg = self.degrees.get(t, 0)
            t_pct = pct_ranks.get(t, 0.5)
            
            deg_min = t_deg * (1 - config.degree_tolerance)
            deg_max = t_deg * (1 + config.degree_tolerance)
            pct_min = max(0, t_pct - config.percentile_window)
            pct_max = min(1, t_pct + config.percentile_window)
            
            for node in m_vector.index:
                if node in targets:
                    continue
                node_deg = self.degrees.get(node, 0)
                node_pct = pct_ranks.get(node, 0)
                
                if (deg_min <= node_deg <= deg_max and
                    pct_min <= node_pct <= pct_max):
                    matched.add(node)
        
        return list(matched)
    
    def compute_redundancy(self, targets: List[str], disease_genes: List[str]) -> float:
        """Compute mean pairwise redundancy on disease projection."""
        target_idx = [self.node_to_idx[t] for t in targets if t in self.node_to_idx]
        disease_idx = [self.node_to_idx[g] for g in disease_genes if g in self.node_to_idx]
        k = len(target_idx)
        
        if k < 2:
            return 0.0
        
        M_D = self.M[disease_idx, :][:, target_idx]
        norms = np.linalg.norm(M_D, axis=0)
        norms[norms == 0] = 1e-10
        M_D_norm = M_D / norms[np.newaxis, :]
        rho = M_D_norm.T @ M_D_norm
        
        return (np.sum(rho) - np.trace(rho)) / (k * (k - 1))
    
    def ecnp(self, targets: List[str], m_vector: pd.Series, 
             disease_genes: List[str], config: ECNPConfig) -> Dict:
        """Run ECNP with a specific disease module."""
        target_m = [m_vector.get(t, 0) for t in targets if t in m_vector.index]
        k = len(target_m)
        
        if k == 0:
            return {'Z': np.nan, 'error': 'No targets'}
        
        I_T = sum(target_m)
        
        pool = self.get_stratum_pool(targets, m_vector, config)
        pool_m = [m_vector.get(p, 0) for p in pool if p in m_vector.index]
        
        if len(pool_m) < 20:
            return {'Z': np.nan, 'error': f'Pool too small ({len(pool_m)})'}
        
        pool_mean = np.mean(pool_m)
        pool_var = np.var(pool_m, ddof=1)
        
        mu_T = k * pool_mean
        mean_rho = self.compute_redundancy(targets, disease_genes)
        
        sigma_T_sq = pool_var / k + config.lambda_redundancy * mean_rho
        sigma_T = np.sqrt(sigma_T_sq)
        
        Z = (I_T - mu_T) / sigma_T
        
        return {
            'Z': Z, 'I_T': I_T, 'mu_T': mu_T, 'sigma_T': sigma_T,
            'k': k, 'pool_size': len(pool_m), 'mean_rho': mean_rho
        }
    
    def run_cross_disease(self, config: ECNPConfig = ECNPConfig()):
        """Run ECNP across all disease modules."""
        results = []
        
        # Load disease modules
        diseases = {}
        for f in DISEASE_DIR.glob("*_lcc.csv"):
            name = f.stem.replace("_lcc", "")
            genes = pd.read_csv(f)['gene_name'].tolist()
            diseases[name] = genes
        
        print("="*70)
        print("CROSS-DISEASE ECNP VALIDATION")
        print(f"Lambda = {config.lambda_redundancy}")
        print("="*70)
        
        for disease, genes in diseases.items():
            print(f"\n{disease.upper()} ({len(genes)} genes)")
            print("-" * 50)
            
            # Compute disease-specific m_j
            m_vector = self.compute_disease_influence(genes)
            
            # Run on both compounds
            for compound, targets in [("Hyperforin", self.hyperforin), 
                                       ("Quercetin", self.quercetin)]:
                result = self.ecnp(targets, m_vector, genes, config)
                
                if not np.isnan(result['Z']):
                    print(f"  {compound}: Z={result['Z']:.2f}, "
                          f"pool={result['pool_size']}, rho={result['mean_rho']:.3f}")
                else:
                    print(f"  {compound}: FAILED ({result.get('error', 'unknown')})")
                
                results.append({
                    'disease': disease,
                    'compound': compound,
                    'n_disease_genes': len(genes),
                    **result
                })
        
        return pd.DataFrame(results)


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    
    analyzer = CrossDiseaseECNP()
    
    # Run with default lambda
    results = analyzer.run_cross_disease()
    
    # Save
    results.to_csv(RESULTS_DIR / "cross_disease_ecnp.csv", index=False)
    print(f"\nResults saved to {RESULTS_DIR / 'cross_disease_ecnp.csv'}")
    
    # Summary statistics
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    
    valid = results[~results['Z'].isna()]
    
    print(f"\nHyperforin Z-scores by disease:")
    hyp = valid[valid['compound'] == 'Hyperforin'].set_index('disease')['Z']
    for d, z in hyp.items():
        print(f"  {d}: {z:.2f}")
    print(f"  Range: [{hyp.min():.2f}, {hyp.max():.2f}]")
    
    print(f"\nQuercetin Z-scores by disease:")
    que = valid[valid['compound'] == 'Quercetin'].set_index('disease')['Z']
    for d, z in que.items():
        print(f"  {d}: {z:.2f}")
    print(f"  Range: [{que.min():.2f}, {que.max():.2f}]")


if __name__ == "__main__":
    main()

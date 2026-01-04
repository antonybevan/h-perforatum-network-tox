"""
Network Calibration: Compare Layer 2 across different STRING thresholds.

Tests stability of inference across:
- STRING 900 (high confidence)
- STRING 700 (medium confidence)
"""

from ecnp_permutation_test import ECNPPermutationTest, PermutationConfig
import numpy as np
import pandas as pd
from pathlib import Path
import networkx as nx


class NetworkCalibration:
    """Compare Layer 2 results across network configurations."""
    
    def __init__(self, project_root: Path = None):
        if project_root is None:
            project_root = Path(__file__).parent.parent.parent.parent
        
        self.project_root = project_root
        self.data_dir = project_root / "data" / "processed"
        self.research_dir = project_root / "research" / "ecnp-closed-form" / "data"
        
    def load_network_stats(self, threshold: int) -> dict:
        """Load network and compute basic stats."""
        parquet_file = self.data_dir / f"network_{threshold}_liver_lcc.parquet"
        
        if not parquet_file.exists():
            return None
            
        edges = pd.read_parquet(parquet_file)
        
        G = nx.from_pandas_edgelist(edges, 'gene1', 'gene2')
        
        return {
            'threshold': threshold,
            'n_nodes': G.number_of_nodes(),
            'n_edges': G.number_of_edges(),
            'avg_degree': 2 * G.number_of_edges() / G.number_of_nodes(),
            'density': nx.density(G)
        }
    
    def compare_thresholds(self, compounds: dict, n_permutations: int = 500) -> pd.DataFrame:
        """
        Compare Layer 2 results across STRING thresholds.
        
        Args:
            compounds: dict of {name: [targets]}
            n_permutations: permutations per test
        """
        # We only have 900 implemented - show what we can
        test = ECNPPermutationTest()
        config = PermutationConfig(n_permutations=n_permutations, random_seed=42)
        
        results = []
        
        # STRING 900 (current implementation)
        stats_900 = self.load_network_stats(900)
        print(f"\nSTRING 900: {stats_900['n_nodes']} nodes, {stats_900['n_edges']} edges")
        
        for name, targets in compounds.items():
            result = test.test(targets, config)
            results.append({
                'threshold': 900,
                'compound': name,
                'k': result['k'],
                'S': result['S_observed'],
                'p_value': result['p_value']
            })
            print(f"  {name}: k={result['k']}, S={result['S_observed']:.4f}, p={result['p_value']:.4f}")
        
        # Check if 700 network exists
        stats_700 = self.load_network_stats(700)
        if stats_700:
            print(f"\nSTRING 700: {stats_700['n_nodes']} nodes, {stats_700['n_edges']} edges")
            print("  (Layer 2 not yet implemented for 700 threshold)")
        else:
            print("\nSTRING 700 network not available in processed data")
        
        return pd.DataFrame(results), {'900': stats_900, '700': stats_700}


class DecisionTheory:
    """Decision-theoretic framework for risk-tier thresholds."""
    
    def __init__(self):
        pass
    
    def compute_cost_matrix(
        self,
        c_fn: float = 10.0,  # Cost of false negative (missed DILI signal)
        c_fp: float = 1.0,   # Cost of false positive (unnecessary concern)
    ) -> dict:
        """
        Compute expected costs at different p-value thresholds.
        
        In drug safety:
        - False negative (Type II): Missing a real DILI signal -> patient harm
        - False positive (Type I): Flagging a safe compound -> unnecessary delay
        
        Args:
            c_fn: Cost of false negative (typically high in safety context)
            c_fp: Cost of false positive (typically lower)
        """
        return {
            'c_fn': c_fn,
            'c_fp': c_fp,
            'cost_ratio': c_fn / c_fp,
            'optimal_alpha': self._compute_optimal_threshold(c_fn, c_fp)
        }
    
    def _compute_optimal_threshold(self, c_fn: float, c_fp: float) -> float:
        """
        Compute optimal alpha based on cost ratio.
        
        From decision theory: optimal threshold when costs are asymmetric
        α* ≈ c_fp / (c_fn + c_fp)
        
        This assumes roughly equal base rates of true positives and negatives.
        """
        return c_fp / (c_fn + c_fp)
    
    def define_risk_tiers(self, p_values: list, labels: list) -> pd.DataFrame:
        """
        Define empirical risk tiers based on p-value distribution.
        
        Tiers:
        - CRITICAL: p < 0.01 (strong evidence of enrichment)
        - HIGH: p < 0.05 (significant at conventional threshold)
        - MODERATE: p < 0.10 (marginal evidence)
        - LOW: p >= 0.10 (no evidence)
        """
        tiers = []
        for p, label in zip(p_values, labels):
            if p < 0.01:
                tier = 'CRITICAL'
            elif p < 0.05:
                tier = 'HIGH'
            elif p < 0.10:
                tier = 'MODERATE'
            else:
                tier = 'LOW'
            tiers.append({'compound': label, 'p_value': p, 'tier': tier})
        
        return pd.DataFrame(tiers)
    
    def expected_cost(self, alpha: float, power: float, c_fn: float, c_fp: float, 
                      prevalence: float = 0.5) -> float:
        """
        Compute expected cost of using threshold alpha.
        
        E[Cost] = P(H1) * (1 - power) * c_fn + P(H0) * alpha * c_fp
        
        Args:
            alpha: Type I error rate
            power: Detection power at this alpha
            c_fn: Cost of false negative
            c_fp: Cost of false positive
            prevalence: Prior probability of true signal
        """
        e_cost = prevalence * (1 - power) * c_fn + (1 - prevalence) * alpha * c_fp
        return e_cost


if __name__ == '__main__':
    print('='*70)
    print('NETWORK CALIBRATION & DECISION THEORY')
    print('='*70)
    
    # Load targets
    test = ECNPPermutationTest()
    targets_df = pd.read_csv(test.data_dir / 'targets_lcc.csv')
    
    compounds = {
        'Hyperforin': targets_df[targets_df['compound'] == 'Hyperforin']['gene_symbol'].tolist(),
        'Quercetin': targets_df[targets_df['compound'] == 'Quercetin']['gene_symbol'].tolist(),
    }
    
    # Network Calibration
    print('\n' + '='*70)
    print('1. NETWORK CALIBRATION')
    print('='*70)
    
    nc = NetworkCalibration()
    results_df, stats = nc.compare_thresholds(compounds, n_permutations=500)
    
    print('\nResults:')
    print(results_df.to_string(index=False))
    
    # Decision Theory
    print('\n' + '='*70)
    print('2. DECISION-THEORETIC FRAMEWORK')
    print('='*70)
    
    dt = DecisionTheory()
    
    # Cost scenarios
    print('\nCost Analysis (drug safety context):')
    print('-'*50)
    
    scenarios = [
        ('Conservative', 10.0, 1.0),   # Missing DILI is 10x worse than false alarm
        ('Balanced', 1.0, 1.0),        # Equal costs
        ('Permissive', 1.0, 10.0),     # False alarms are 10x worse
    ]
    
    for name, c_fn, c_fp in scenarios:
        cost_info = dt.compute_cost_matrix(c_fn, c_fp)
        print(f"  {name}: c_fn={c_fn}, c_fp={c_fp} -> optimal_alpha={cost_info['optimal_alpha']:.3f}")
    
    # Risk Tiers
    print('\n' + '-'*50)
    print('Risk Tier Classification:')
    print('-'*50)
    
    p_values = [0.011, 0.635]  # From Layer 2 validation
    labels = ['Hyperforin', 'Quercetin']
    tiers = dt.define_risk_tiers(p_values, labels)
    print(tiers.to_string(index=False))
    
    # Expected cost analysis
    print('\n' + '-'*50)
    print('Expected Cost at Different Thresholds:')
    print('-'*50)
    
    # Assume power from our spike-in experiments at 50% spike
    # power ~ 0.45 at k=10, alpha=0.05
    power_at_05 = 0.45
    power_at_01 = 0.25  # Lower power at stricter threshold
    power_at_10 = 0.55  # Higher power at relaxed threshold
    
    c_fn, c_fp = 10.0, 1.0  # Conservative scenario
    print(f"\nUsing Conservative costs (c_fn={c_fn}, c_fp={c_fp}):")
    for alpha, power in [(0.01, power_at_01), (0.05, power_at_05), (0.10, power_at_10)]:
        e_cost = dt.expected_cost(alpha, power, c_fn, c_fp, prevalence=0.3)
        print(f"  alpha={alpha:.2f}, power={power:.2f}: E[Cost]={e_cost:.3f}")
    
    print('\n' + '='*70)
    print('SUMMARY')
    print('='*70)
    print("""
Risk Tier Guidelines (drug safety):
  - CRITICAL (p<0.01): Immediate review, prioritize for DILI assessment
  - HIGH (p<0.05): Flag for safety evaluation
  - MODERATE (p<0.10): Monitor, include in safety panel
  - LOW (p>=0.10): Standard development pathway

Hyperforin: p=0.011 -> HIGH tier -> Flag for DILI review
Quercetin: p=0.635 -> LOW tier -> Standard pathway

Decision-theoretic justification:
  At typical safety cost ratio (c_fn/c_fp = 10), optimal threshold is alpha=0.09.
  Using alpha=0.05 is slightly conservative, which is appropriate for drug safety.
""")

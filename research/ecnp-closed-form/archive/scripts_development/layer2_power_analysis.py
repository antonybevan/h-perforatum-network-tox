"""
Layer 2 Power Analysis: Effect size vs detection probability.

Spike-in experiment: inject known signal by replacing random targets
with high-influence nodes, measure detection rate.
"""

from ecnp_permutation_test import ECNPPermutationTest, PermutationConfig
import numpy as np
import pandas as pd
from pathlib import Path


class PowerAnalysis:
    """Power analysis via spike-in experiments."""
    
    def __init__(self):
        self.test = ECNPPermutationTest()
        self.m_vector = self.test.m_vector
        self.nodes = self.test.node_list
        
        # Pre-compute influence percentiles
        vals = self.test.m_array
        self.percentiles = {
            node: np.mean(vals <= self.m_vector.get(node, 0))
            for node in self.nodes
        }
        
        # Define high-influence nodes (top 5%)
        self.high_influence_nodes = [
            n for n in self.nodes if self.percentiles[n] >= 0.95
        ]
        
        # Define low-influence nodes (bottom 50%)
        self.low_influence_nodes = [
            n for n in self.nodes if self.percentiles[n] <= 0.50
        ]
    
    def run_spike_in(
        self,
        k: int = 10,
        spike_fractions: list = None,
        n_trials: int = 100,
        n_permutations: int = 1000,
        alpha: float = 0.05,
        seed: int = 42
    ) -> pd.DataFrame:
        """
        Spike-in experiment: measure power at different effect sizes.
        
        Args:
            k: Number of targets per synthetic compound
            spike_fractions: Fraction of targets replaced with high-influence nodes
            n_trials: Number of random trials per spike fraction
            n_permutations: Permutations per test
            alpha: Significance level
            seed: Random seed
            
        Returns:
            DataFrame with columns: spike_fraction, power, mean_p, std_p
        """
        if spike_fractions is None:
            spike_fractions = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        
        rng = np.random.default_rng(seed)
        config = PermutationConfig(n_permutations=n_permutations, random_seed=seed)
        
        results = []
        
        for frac in spike_fractions:
            n_high = int(k * frac)
            n_low = k - n_high
            
            p_values = []
            
            for trial in range(n_trials):
                try:
                    # Sample targets: n_low from low-influence, n_high from high-influence
                    low_sample = rng.choice(self.low_influence_nodes, size=n_low, replace=False) if n_low > 0 else []
                    high_sample = rng.choice(self.high_influence_nodes, size=n_high, replace=False) if n_high > 0 else []
                    targets = list(low_sample) + list(high_sample)
                    
                    # Run Layer 2 test
                    result = self.test.test(targets, config)
                    p_values.append(result['p_value'])
                except Exception as e:
                    # Skip failed trials (strata issues)
                    continue
            
            if len(p_values) == 0:
                print(f"  spike={frac:.1f} (n_high={n_high}): FAILED - no valid trials")
                continue
                
            p_values = np.array(p_values)
            power = np.mean(p_values < alpha)
            
            results.append({
                'spike_fraction': frac,
                'n_high': n_high,
                'n_low': n_low,
                'power': power,
                'mean_p': np.mean(p_values),
                'std_p': np.std(p_values),
                'median_p': np.median(p_values),
                'n_valid': len(p_values)
            })
            
            print(f"  spike={frac:.1f} (n_high={n_high}): power={power:.3f}, mean_p={np.mean(p_values):.3f} ({len(p_values)} trials)")
        
        return pd.DataFrame(results)
    
    def run_effect_size_curve(
        self,
        k_values: list = None,
        n_trials: int = 50,
        n_permutations: int = 500,
        seed: int = 42
    ) -> pd.DataFrame:
        """
        Effect size curve: power vs k at fixed spike fraction.
        
        Fixes spike_fraction=0.5, varies k to see how target count affects power.
        """
        if k_values is None:
            k_values = [5, 10, 15, 20, 30, 50]
        
        rng = np.random.default_rng(seed)
        config = PermutationConfig(n_permutations=n_permutations, random_seed=seed)
        
        results = []
        
        for k in k_values:
            n_high = k // 2
            n_low = k - n_high
            
            p_values = []
            
            for trial in range(n_trials):
                try:
                    low_sample = rng.choice(self.low_influence_nodes, size=n_low, replace=False) if n_low > 0 else []
                    high_sample = rng.choice(self.high_influence_nodes, size=n_high, replace=False) if n_high > 0 else []
                    targets = list(low_sample) + list(high_sample)
                    
                    result = self.test.test(targets, config)
                    p_values.append(result['p_value'])
                except Exception:
                    continue
            
            if len(p_values) == 0:
                print(f"  k={k}: FAILED")
                continue
                
            p_values = np.array(p_values)
            power = np.mean(p_values < 0.05)
            
            results.append({
                'k': k,
                'n_high': n_high,
                'power': power,
                'mean_p': np.mean(p_values),
                'n_valid': len(p_values)
            })
            
            print(f"  k={k} (n_high={n_high}): power={power:.3f} ({len(p_values)} trials)")
        
        return pd.DataFrame(results)


if __name__ == '__main__':
    print('='*70)
    print('LAYER 2 POWER ANALYSIS')
    print('='*70)
    
    pa = PowerAnalysis()
    
    print(f"\nHigh-influence pool: {len(pa.high_influence_nodes)} nodes (top 5%)")
    print(f"Low-influence pool: {len(pa.low_influence_nodes)} nodes (bottom 50%)")
    
    # Spike-in experiment
    print('\n' + '-'*70)
    print('SPIKE-IN EXPERIMENT (k=10, varying spike fraction)')
    print('-'*70)
    
    spike_results = pa.run_spike_in(
        k=10,
        spike_fractions=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
        n_trials=50,
        n_permutations=200,
        seed=42
    )
    
    print('\n' + '-'*70)
    print('EFFECT SIZE CURVE (spike=50%, varying k)')
    print('-'*70)
    
    effect_results = pa.run_effect_size_curve(
        k_values=[5, 10, 20, 30],
        n_trials=30,
        n_permutations=200,
        seed=42
    )
    
    print('\n' + '='*70)
    print('SUMMARY')
    print('='*70)
    print('\nSpike-in Results:')
    print(spike_results[['spike_fraction', 'n_high', 'power', 'mean_p']].to_string(index=False))
    
    print('\nEffect Size Results:')
    print(effect_results[['k', 'n_high', 'power', 'mean_p']].to_string(index=False))
    
    # Save results
    output_dir = Path(__file__).parent.parent / 'results'
    output_dir.mkdir(exist_ok=True)
    spike_results.to_csv(output_dir / 'power_spike_in.csv', index=False)
    effect_results.to_csv(output_dir / 'power_effect_size.csv', index=False)
    print(f"\nResults saved to {output_dir}")

"""
Validate spectral variance against empirical resampling.

We need to verify that sigma_exact matches the empirical variance
from stratum-matched null sampling.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import numpy as np
import pandas as pd
from ecnp_optimized import ECNPOptimized, ECNPConfig


def empirical_variance_test():
    """Compare exact variance to empirical."""
    
    print("EMPIRICAL VARIANCE VALIDATION")
    print("=" * 60)
    
    np.random.seed(42)
    ecnp = ECNPOptimized()
    config = ECNPConfig()
    
    # Load targets
    df = pd.read_csv(ecnp.data_dir / "targets_lcc.csv")
    hyperforin = df[df['compound'] == 'Hyperforin']['gene_symbol'].tolist()
    quercetin = df[df['compound'] == 'Quercetin']['gene_symbol'].tolist()
    
    for name, targets in [("Hyperforin", hyperforin), ("Quercetin", quercetin)]:
        print(f"\n{name}:")
        
        # Get target indices
        target_indices = np.array([ecnp.node_to_idx[t] for t in targets 
                                   if t in ecnp.node_to_idx])
        k = len(target_indices)
        
        # Get target properties
        target_info = []
        for t_idx in target_indices:
            target_info.append({
                'idx': t_idx,
                'deg': ecnp.degrees_array[t_idx],
                'pct': ecnp.percentiles_array[t_idx]
            })
        
        # Resample stratum-matched null many times
        n_resamples = 200
        null_Is = []
        
        for _ in range(n_resamples):
            matched_idx = []
            matched_nodes = set()
            
            for t in target_info:
                deg_min = t['deg'] * (1 - config.degree_tolerance)
                deg_max = t['deg'] * (1 + config.degree_tolerance)
                pct_min = max(0, t['pct'] - config.percentile_window)
                pct_max = min(1, t['pct'] + config.percentile_window)
                
                candidates = []
                for j in range(ecnp.n_nodes):
                    node = ecnp.node_list[j]
                    if j in target_indices or node in matched_nodes:
                        continue
                    if deg_min <= ecnp.degrees_array[j] <= deg_max:
                        if pct_min <= ecnp.percentiles_array[j] <= pct_max:
                            candidates.append(j)
                
                if candidates:
                    chosen = np.random.choice(candidates)
                    matched_idx.append(chosen)
                    matched_nodes.add(ecnp.node_list[chosen])
            
            if len(matched_idx) == k:
                I_null = np.sum(ecnp.m_array[matched_idx])
                null_Is.append(I_null)
        
        null_Is = np.array(null_Is)
        
        # Empirical statistics
        emp_mean = np.mean(null_Is)
        emp_sigma = np.std(null_Is, ddof=1)
        
        # Get algorithm result
        result = ecnp.compute(targets)
        
        # Compute Z with empirical sigma
        I_T = result['I_T']
        Z_empirical = (I_T - emp_mean) / emp_sigma
        
        print(f"  k = {k}")
        print(f"  n_resamples = {len(null_Is)}")
        print()
        print(f"  Algorithm mu_T:     {result['mu_T']:.4f}")
        print(f"  Empirical mean:     {emp_mean:.4f}")
        print()
        print(f"  Closed-form sigma:  {result['sigma_T']:.4f}")
        print(f"  Empirical sigma:    {emp_sigma:.4f}")
        print(f"  Ratio (emp/CF):     {emp_sigma / result['sigma_T']:.2f}x")
        print()
        print(f"  I_T:                {I_T:.4f}")
        print(f"  Z_raw (CF):         {result['Z_raw']:.2f}")
        print(f"  Z_adj (VIF):        {result['Z']:.2f}")
        print(f"  Z_empirical (gold): {Z_empirical:.2f}")
        print()
        
        # What should the correct VIF be?
        correct_VIF = (emp_sigma / result['sigma_T']) ** 2
        print(f"  Correct VIF = (emp/CF)^2 = {correct_VIF:.2f}")
        print(f"  Current VIF = {result['VIF']:.2f}")
    
    print("\n" + "=" * 60)
    print("CONCLUSION:")
    print("  Z_empirical is the ground truth.")
    print("  Correct VIF should make Z_adj = Z_empirical")


if __name__ == "__main__":
    empirical_variance_test()

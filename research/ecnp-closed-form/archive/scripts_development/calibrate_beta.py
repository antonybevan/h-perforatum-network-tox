"""
Calibrate shrinkage parameter beta for VIF.

Correct formula: VIF = 1 + beta * (k-1) * rho_target

We need to find beta such that:
  sigma_CF * sqrt(VIF) = sigma_empirical

Using multiple calibration points (different compounds, k values).
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import numpy as np
import pandas as pd
from ecnp_optimized import ECNPOptimized, ECNPConfig


def calibrate_beta():
    """Calibrate beta to match empirical variance."""
    
    print("CALIBRATING SHRINKAGE PARAMETER BETA")
    print("=" * 60)
    
    np.random.seed(42)
    ecnp = ECNPOptimized()
    config = ECNPConfig()
    
    # Load targets
    df = pd.read_csv(ecnp.data_dir / "targets_lcc.csv")
    
    calibration_data = []
    
    # Test multiple compounds
    compounds = ['Hyperforin', 'Quercetin']
    
    for compound in compounds:
        targets = df[df['compound'] == compound]['gene_symbol'].tolist()
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
        
        # Empirical resampling (200 samples)
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
        
        emp_sigma = np.std(null_Is, ddof=1)
        
        # Get closed-form result
        result = ecnp.compute(targets)
        cf_sigma = result['sigma_T']
        
        # Target correlation (rho)
        rho = result['mean_rho']
        
        # VIF needed to match empirical
        VIF_needed = (emp_sigma / cf_sigma) ** 2
        
        # beta = (VIF - 1) / ((k-1) * rho)
        if (k - 1) * rho > 0:
            beta = (VIF_needed - 1) / ((k - 1) * rho)
        else:
            beta = 0
        
        calibration_data.append({
            'compound': compound,
            'k': k,
            'rho': rho,
            'cf_sigma': cf_sigma,
            'emp_sigma': emp_sigma,
            'VIF_needed': VIF_needed,
            'beta': beta
        })
        
        print(f"\n{compound}:")
        print(f"  k = {k}, rho = {rho:.4f}")
        print(f"  sigma_CF = {cf_sigma:.4f}")
        print(f"  sigma_emp = {emp_sigma:.4f}")
        print(f"  VIF_needed = {VIF_needed:.2f}")
        print(f"  => beta = {beta:.2f}")
    
    # Also test random k-sized samples
    print("\n\nTesting random target sets for calibration...")
    
    for k_test in [5, 10, 20, 30]:
        # Sample random nodes
        for trial in range(3):
            rand_idx = np.random.choice(ecnp.n_nodes, k_test, replace=False)
            rand_targets = [ecnp.node_list[i] for i in rand_idx]
            
            result = ecnp.compute(rand_targets)
            if result['status'].value != 'success':
                continue
            
            cf_sigma = result['sigma_T']
            rho = result['mean_rho']
            
            # Quick empirical estimate (50 samples)
            target_info = []
            for t_idx in rand_idx:
                target_info.append({
                    'idx': t_idx,
                    'deg': ecnp.degrees_array[t_idx],
                    'pct': ecnp.percentiles_array[t_idx]
                })
            
            null_Is = []
            for _ in range(50):
                matched_idx = []
                matched_nodes = set()
                
                for t in target_info:
                    deg_min = t['deg'] * (1 - config.degree_tolerance)
                    deg_max = t['deg'] * (1 + config.degree_tolerance)
                    pct_min = max(0, t['pct'] - config.percentile_window)
                    pct_max = min(1, t['pct'] + config.percentile_window)
                    
                    candidates = []
                    for j in range(ecnp.n_nodes):
                        if j in rand_idx or ecnp.node_list[j] in matched_nodes:
                            continue
                        if deg_min <= ecnp.degrees_array[j] <= deg_max:
                            if pct_min <= ecnp.percentiles_array[j] <= pct_max:
                                candidates.append(j)
                    
                    if candidates:
                        chosen = np.random.choice(candidates)
                        matched_idx.append(chosen)
                        matched_nodes.add(ecnp.node_list[chosen])
                
                if len(matched_idx) == k_test:
                    I_null = np.sum(ecnp.m_array[matched_idx])
                    null_Is.append(I_null)
            
            if len(null_Is) >= 20:
                emp_sigma = np.std(null_Is, ddof=1)
                VIF_needed = (emp_sigma / cf_sigma) ** 2
                if (k_test - 1) * rho > 0:
                    beta = (VIF_needed - 1) / ((k_test - 1) * rho)
                    calibration_data.append({
                        'compound': f'random_k{k_test}_{trial}',
                        'k': k_test,
                        'rho': rho,
                        'cf_sigma': cf_sigma,
                        'emp_sigma': emp_sigma,
                        'VIF_needed': VIF_needed,
                        'beta': beta
                    })
    
    # Summary
    print("\n" + "=" * 60)
    print("CALIBRATION SUMMARY")
    print("=" * 60)
    
    betas = [d['beta'] for d in calibration_data if d['beta'] > 0 and d['beta'] < 100]
    
    print(f"\nAll calibrated betas:")
    for d in calibration_data:
        print(f"  {d['compound']}: k={d['k']}, rho={d['rho']:.3f}, beta={d['beta']:.2f}")
    
    if betas:
        mean_beta = np.mean(betas)
        median_beta = np.median(betas)
        print(f"\nMean beta: {mean_beta:.2f}")
        print(f"Median beta: {median_beta:.2f}")
        
        print(f"\nRECOMMENDATION:")
        print(f"  Use VIF = 1 + {median_beta:.1f} * (k-1) * rho")
    
    return calibration_data


if __name__ == "__main__":
    calibrate_beta()

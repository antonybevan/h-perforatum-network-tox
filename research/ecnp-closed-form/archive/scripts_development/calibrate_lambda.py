"""
Lambda Recalibration for Variance Match

Objective: Find λ* such that closed-form variance matches empirical variance.

Method:
1. Generate stratum-matched null samples across k, percentiles
2. Compute empirical variance of I(T) for each sample
3. Compute closed-form variance using current λ
4. Solve for λ* such that E[σ²_CF(λ*)] ≈ E[σ²_empirical]

This is ONE λ*, per network, per disease module.
No dependence on compound identity.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import numpy as np
import pandas as pd
from scipy import optimize
from ecnp_optimized import ECNPOptimized, ECNPConfig, ECNPStatus


def generate_stratum_samples(ecnp, k_values=[5, 10, 20, 30], n_per_k=50):
    """Generate stratum-matched null samples across k values."""
    
    np.random.seed(42)
    samples = []
    
    # Get percentile distribution to sample from
    # Sample strata proportionally to real compounds
    targets_df = pd.read_csv(ecnp.data_dir / "targets_lcc.csv")
    
    for k in k_values:
        print(f"  Generating samples for k={k}...")
        
        for trial in range(n_per_k):
            # Sample k nodes with realistic percentile distribution
            # Use uniform sampling across percentile strata
            sampled_targets = []
            
            # Sample from different percentile bands
            for _ in range(k):
                # Pick a random percentile target
                pct_target = np.random.uniform(0.1, 0.99)
                
                # Find nodes in that stratum
                candidates = []
                for i in range(ecnp.n_nodes):
                    if ecnp.node_list[i] in sampled_targets:
                        continue
                    if abs(ecnp.percentiles_array[i] - pct_target) < 0.1:
                        candidates.append(i)
                
                if candidates:
                    idx = np.random.choice(candidates)
                    sampled_targets.append(ecnp.node_list[idx])
            
            if len(sampled_targets) >= 2:
                samples.append({
                    'k': len(sampled_targets),
                    'targets': sampled_targets
                })
    
    return samples


def compute_empirical_variances(ecnp, samples, n_resamples=50):
    """
    For each sample, compute empirical variance via stratum resampling.
    """
    print("\nComputing empirical variances...")
    
    results = []
    
    for i, sample in enumerate(samples):
        if i % 20 == 0:
            print(f"  Processing sample {i+1}/{len(samples)}...")
        
        targets = sample['targets']
        k = sample['k']
        
        # Get target indices and properties
        target_info = []
        for t in targets:
            if t in ecnp.node_to_idx:
                idx = ecnp.node_to_idx[t]
                target_info.append({
                    'node': t,
                    'idx': idx,
                    'deg': ecnp.degrees_array[idx],
                    'pct': ecnp.percentiles_array[idx]
                })
        
        if len(target_info) < 2:
            continue
        
        # Resample stratum-matched nulls
        null_Is = []
        
        for _ in range(n_resamples):
            matched_idx = []
            matched_nodes = set()
            
            for t_info in target_info:
                deg_min = t_info['deg'] * 0.8
                deg_max = t_info['deg'] * 1.2
                pct_min = max(0, t_info['pct'] - 0.1)
                pct_max = min(1, t_info['pct'] + 0.1)
                
                candidates = []
                for j in range(ecnp.n_nodes):
                    node = ecnp.node_list[j]
                    if node in matched_nodes or node in targets:
                        continue
                    if deg_min <= ecnp.degrees_array[j] <= deg_max:
                        if pct_min <= ecnp.percentiles_array[j] <= pct_max:
                            candidates.append(j)
                
                if candidates:
                    chosen = np.random.choice(candidates)
                    matched_idx.append(chosen)
                    matched_nodes.add(ecnp.node_list[chosen])
            
            if matched_idx:
                I_null = np.sum(ecnp.m_array[matched_idx])
                null_Is.append(I_null)
        
        if len(null_Is) >= 10:
            emp_var = np.var(null_Is, ddof=1)
            
            # Compute closed-form quantities
            pool_m = []
            for t_info in target_info:
                deg_min = t_info['deg'] * 0.8
                deg_max = t_info['deg'] * 1.2
                pct_min = max(0, t_info['pct'] - 0.1)
                pct_max = min(1, t_info['pct'] + 0.1)
                
                for j in range(ecnp.n_nodes):
                    if ecnp.node_list[j] in targets:
                        continue
                    if deg_min <= ecnp.degrees_array[j] <= deg_max:
                        if pct_min <= ecnp.percentiles_array[j] <= pct_max:
                            pool_m.append(ecnp.m_array[j])
            
            if pool_m:
                pool_var = np.var(pool_m, ddof=1)
                
                # Compute mean rho for this sample
                target_idx = np.array([t['idx'] for t in target_info])
                M_targets = ecnp.M_dili[:, target_idx]
                norms = np.linalg.norm(M_targets, axis=0, keepdims=True)
                norms = np.where(norms == 0, 1e-10, norms)
                M_norm = M_targets / norms
                rho_matrix = M_norm.T @ M_norm
                
                if len(target_idx) > 1:
                    mask = np.ones_like(rho_matrix, dtype=bool)
                    np.fill_diagonal(mask, False)
                    mean_rho = np.mean(rho_matrix[mask])
                else:
                    mean_rho = 0
                
                results.append({
                    'k': len(target_info),
                    'emp_var': emp_var,
                    'pool_var': pool_var,
                    'mean_rho': mean_rho
                })
    
    return results


def solve_for_lambda(results):
    """
    Solve for λ* such that E[σ²_CF(λ*)] ≈ E[σ²_empirical]
    
    σ²_CF = pool_var/k + λ * ρ
    
    We want: E[pool_var/k + λ * ρ] = E[emp_var]
    
    Rearranging: λ = (E[emp_var] - E[pool_var/k]) / E[ρ]
    """
    
    # Compute expected values
    emp_vars = np.array([r['emp_var'] for r in results])
    pool_term = np.array([r['pool_var'] / r['k'] for r in results])
    rhos = np.array([r['mean_rho'] for r in results])
    
    E_emp_var = np.mean(emp_vars)
    E_pool_term = np.mean(pool_term)
    E_rho = np.mean(rhos)
    
    # Solve for λ*
    lambda_star = (E_emp_var - E_pool_term) / E_rho if E_rho > 0 else 0.02
    
    print(f"\nCalibration Results:")
    print(f"  E[emp_var]:    {E_emp_var:.6f}")
    print(f"  E[pool_var/k]: {E_pool_term:.6f}")
    print(f"  E[rho]:        {E_rho:.4f}")
    print(f"  λ* solved:     {lambda_star:.4f}")
    
    # Verify: compute variance inflation factor with new λ
    cf_vars_new = pool_term + lambda_star * rhos
    kappa = np.mean(emp_vars / cf_vars_new)
    
    print(f"\n  With λ* = {lambda_star:.4f}:")
    print(f"    E[emp_var / CF_var]: {kappa:.2f} (target: 1.0)")
    
    return lambda_star, kappa


def verify_calibration(ecnp, lambda_star):
    """Verify the new λ on validation compounds."""
    
    print("\n" + "=" * 60)
    print("VERIFICATION ON VALIDATION COMPOUNDS")
    print("=" * 60)
    
    # Create new config with calibrated λ
    config = ECNPConfig()
    config.lambda_redundancy = lambda_star
    
    targets_df = pd.read_csv(ecnp.data_dir / "targets_lcc.csv")
    hyperforin = targets_df[targets_df['compound'] == 'Hyperforin']['gene_symbol'].tolist()
    quercetin = targets_df[targets_df['compound'] == 'Quercetin']['gene_symbol'].tolist()
    
    # Original config
    orig_config = ECNPConfig()
    
    for name, targets in [("Hyperforin", hyperforin), ("Quercetin", quercetin)]:
        r_orig = ecnp.compute(targets, orig_config)
        r_new = ecnp.compute(targets, config)
        
        print(f"\n{name}:")
        print(f"  λ_old = {orig_config.lambda_redundancy:.4f}: Z = {r_orig['Z']:.2f}, σ = {r_orig['sigma_T']:.4f}")
        print(f"  λ_new = {lambda_star:.4f}: Z = {r_new['Z']:.2f}, σ = {r_new['sigma_T']:.4f}")
        
        # Compare to MC reference
        mc_ref = 10.27 if name == "Hyperforin" else 4.42
        error_new = abs(r_new['Z'] - mc_ref) / mc_ref * 100
        print(f"  MC ref: {mc_ref:.2f}, Error with λ*: {error_new:.1f}%")


def main():
    print("=" * 60)
    print("LAMBDA RECALIBRATION FOR VARIANCE MATCH")
    print("=" * 60)
    
    ecnp = ECNPOptimized()
    
    print("\n1. GENERATING STRATUM SAMPLES")
    print("-" * 40)
    samples = generate_stratum_samples(ecnp, k_values=[5, 10, 15, 20, 30], n_per_k=30)
    print(f"   Generated {len(samples)} samples")
    
    print("\n2. COMPUTING EMPIRICAL VARIANCES")
    print("-" * 40)
    results = compute_empirical_variances(ecnp, samples, n_resamples=30)
    print(f"   Computed variances for {len(results)} samples")
    
    print("\n3. SOLVING FOR lambda*")
    print("-" * 40)
    lambda_star, kappa = solve_for_lambda(results)
    
    verify_calibration(ecnp, lambda_star)
    
    print("\n" + "=" * 60)
    print("FINAL RECOMMENDATION")
    print("=" * 60)
    print(f"""
    Calibrated lambda* = {lambda_star:.4f}
    
    Update ECNPConfig with:
        lambda_redundancy: float = {lambda_star:.4f}
    
    This lambda* is calibrated to match empirical null variance.
    It is network-level, disease-module-level constant.
    No per-compound tuning.
    """)
    
    return lambda_star


if __name__ == "__main__":
    main()

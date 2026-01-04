"""
Theoretical Investigation: Why VIF ≈ 10?

The closed-form variance formula is:
  σ²_CF = pool_var / k + λ * mean_rho

But empirical variance is ~10x higher. Why?

Let's derive the EXACT variance from first principles.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import numpy as np
import pandas as pd
from ecnp_optimized import ECNPOptimized, ECNPConfig


def investigate_vif():
    """Investigate the theoretical basis for VIF ≈ 10."""
    
    print("=" * 70)
    print("THEORETICAL INVESTIGATION: WHY VIF ~ 10?")
    print("=" * 70)
    
    ecnp = ECNPOptimized()
    config = ECNPConfig()
    
    # Load Hyperforin targets
    df = pd.read_csv(ecnp.data_dir / "targets_lcc.csv")
    targets = df[df['compound'] == 'Hyperforin']['gene_symbol'].tolist()
    
    target_indices = np.array([ecnp.node_to_idx[t] for t in targets 
                               if t in ecnp.node_to_idx])
    k = len(target_indices)
    
    print(f"\nHyperforin: k = {k} targets")
    
    # Get pool indices
    pool_mask = np.zeros(ecnp.n_nodes, dtype=bool)
    for t_idx in target_indices:
        t_deg = ecnp.degrees_array[t_idx]
        t_pct = ecnp.percentiles_array[t_idx]
        
        deg_min = t_deg * (1 - config.degree_tolerance)
        deg_max = t_deg * (1 + config.degree_tolerance)
        pct_min = max(0.0, t_pct - config.percentile_window)
        pct_max = min(1.0, t_pct + config.percentile_window)
        
        deg_ok = (ecnp.degrees_array >= deg_min) & (ecnp.degrees_array <= deg_max)
        pct_ok = (ecnp.percentiles_array >= pct_min) & (ecnp.percentiles_array <= pct_max)
        
        pool_mask |= (deg_ok & pct_ok)
    
    pool_indices = np.where(pool_mask)[0]
    pool_indices = pool_indices[~np.isin(pool_indices, target_indices)]
    N = len(pool_indices)
    
    print(f"Pool: N = {N} nodes")
    
    # Get m values
    pool_m = ecnp.m_array[pool_indices]
    target_m = ecnp.m_array[target_indices]
    
    # Pool statistics
    pool_mean = np.mean(pool_m)
    pool_var = np.var(pool_m, ddof=1)
    pool_sigma = np.sqrt(pool_var)
    
    print(f"\n--- POOL STATISTICS ---")
    print(f"Pool mean:  {pool_mean:.4f}")
    print(f"Pool sigma: {pool_sigma:.4f}")
    print(f"Pool var:   {pool_var:.6f}")
    
    # Current closed-form variance
    # σ²_CF = pool_var / k + λ * mean_rho
    
    # Compute mean_rho for targets
    M_targets = ecnp.M_dili[:, target_indices]
    norms = np.linalg.norm(M_targets, axis=0, keepdims=True)
    norms = np.where(norms == 0, 1e-10, norms)
    M_norm = M_targets / norms
    rho_matrix = M_norm.T @ M_norm
    mask = ~np.eye(k, dtype=bool)
    mean_rho = np.mean(rho_matrix[mask])
    
    lambda_val = config.lambda_redundancy
    sigma_sq_CF = pool_var / k + lambda_val * mean_rho
    sigma_CF = np.sqrt(sigma_sq_CF)
    
    print(f"\n--- CLOSED-FORM VARIANCE ---")
    print(f"pool_var / k:  {pool_var/k:.6f}")
    print(f"lambda * rho:  {lambda_val * mean_rho:.6f}")
    print(f"sigma²_CF:     {sigma_sq_CF:.6f}")
    print(f"sigma_CF:      {sigma_CF:.4f}")
    
    # EXACT variance of sum of k samples from pool
    # 
    # For the sum of k CORRELATED random variables:
    # Var[X_1 + ... + X_k] = sum_i Var[X_i] + 2 * sum_{i<j} Cov[X_i, X_j]
    #                      = k * sigma² + k*(k-1) * sigma² * rho_avg
    #                      = k * sigma² * [1 + (k-1) * rho_avg]
    # 
    # This is the EXACT formula for equi-correlated samples!
    
    # Compute pool rho_avg
    M_pool = ecnp.M_dili[:, pool_indices]
    norms_pool = np.linalg.norm(M_pool, axis=0, keepdims=True)
    norms_pool = np.where(norms_pool == 0, 1e-10, norms_pool)
    M_pool_norm = M_pool / norms_pool
    
    # Sample rho_avg (full matrix is N x N, too large)
    np.random.seed(42)
    sample_size = min(500, N)
    sample_idx = np.random.choice(N, sample_size, replace=False)
    M_sample = M_pool_norm[:, sample_idx]
    rho_sample = M_sample.T @ M_sample
    mask_sample = ~np.eye(sample_size, dtype=bool)
    rho_pool_avg = np.mean(rho_sample[mask_sample])
    
    print(f"\n--- EXACT VARIANCE (equi-correlated formula) ---")
    print(f"Pool rho_avg:  {rho_pool_avg:.4f}")
    
    # Exact variance for sum of k equi-correlated samples
    sigma_sq_exact = k * pool_var * (1 + (k - 1) * rho_pool_avg)
    sigma_exact = np.sqrt(sigma_sq_exact)
    
    print(f"sigma²_exact = k * pool_var * [1 + (k-1)*rho]")
    print(f"            = {k} * {pool_var:.6f} * [1 + {k-1}*{rho_pool_avg:.4f}]")
    print(f"            = {sigma_sq_exact:.6f}")
    print(f"sigma_exact:   {sigma_exact:.4f}")
    
    # Compare to empirical
    print(f"\n--- COMPARISON ---")
    print(f"sigma_CF:      {sigma_CF:.4f}")
    print(f"sigma_exact:   {sigma_exact:.4f}")
    print(f"sigma_emp:     ~0.176 (from 200 resamples)")
    print(f"Ratio exact/CF: {sigma_exact/sigma_CF:.2f}x")
    
    # THE KEY INSIGHT
    print(f"\n" + "=" * 70)
    print("THE KEY INSIGHT")
    print("=" * 70)
    print("""
The closed-form formula uses:
  σ²_CF = pool_var / k + λ * mean_rho

But the CORRECT formula for sum of k equi-correlated samples is:
  σ²_exact = k * pool_var * [1 + (k-1) * rho_avg]

The error is in how variance scales with k:
  - WRONG: pool_var / k  (variance DECREASES as k increases)
  - RIGHT: k * pool_var  (variance INCREASES as k increases)

The current formula assumes we're estimating the MEAN of k samples.
But we're actually computing the SUM of k samples.

For the SUM:
  - Var[sum of k independent] = k * sigma²
  - Var[sum of k correlated] = k * sigma² * [1 + (k-1)*rho]
""")
    
    # Correct VIF
    VIF_ratio = sigma_sq_exact / sigma_sq_CF
    print(f"\nTheoretical VIF = σ²_exact / σ²_CF = {VIF_ratio:.1f}")
    
    # Derive VIF formula
    # VIF = k * pool_var * [1 + (k-1)*rho] / (pool_var/k + λ*rho)
    #     ≈ k² * [1 + (k-1)*rho] / [1 + k*λ*rho/pool_var]
    
    print(f"\n--- THEORETICAL VIF FORMULA ---")
    print(f"VIF = k² * [1 + (k-1)*rho] / [1 + k*λ*rho/pool_var]")
    print(f"    = {k}² * [1 + {k-1}*{rho_pool_avg:.4f}] / [1 + {k}*{lambda_val}*{mean_rho:.4f}/{pool_var:.6f}]")
    
    numerator = k**2 * (1 + (k-1) * rho_pool_avg)
    denominator = 1 + k * lambda_val * mean_rho / pool_var
    VIF_formula = numerator / denominator
    print(f"    = {numerator:.2f} / {denominator:.2f}")
    print(f"    = {VIF_formula:.1f}")
    
    return VIF_formula


if __name__ == "__main__":
    investigate_vif()

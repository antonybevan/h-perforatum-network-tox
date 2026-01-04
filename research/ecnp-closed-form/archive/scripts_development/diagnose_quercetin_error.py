"""Investigate Quercetin error and find optimal parameters."""
import numpy as np
import pandas as pd
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from ecnp_optimized import ECNPOptimized, ECNPConfig

def main():
    ecnp = ECNPOptimized()
    targets_df = pd.read_csv(ecnp.data_dir / "targets_lcc.csv")
    quercetin = targets_df[targets_df['compound'] == 'Quercetin']['gene_symbol'].tolist()
    hyperforin = targets_df[targets_df['compound'] == 'Hyperforin']['gene_symbol'].tolist()
    
    ref_que = 4.42
    ref_hyp = 10.27
    
    print("=" * 70)
    print("QUERCETIN ERROR ANALYSIS")
    print("=" * 70)
    
    # 1. Lambda sensitivity
    print("\n1. LAMBDA SENSITIVITY:")
    print(f"{'Lambda':<10} | {'Hyp Z':<8} | {'Hyp Err':<8} | {'Que Z':<8} | {'Que Err':<8}")
    print("-" * 55)
    
    for lam in [0.005, 0.01, 0.015, 0.0195, 0.025, 0.03, 0.04, 0.05, 0.06, 0.08]:
        config = ECNPConfig(lambda_redundancy=lam)
        r_h = ecnp.compute(hyperforin, config)
        r_q = ecnp.compute(quercetin, config)
        hyp_err = abs(r_h['Z'] - ref_hyp) / ref_hyp * 100
        que_err = abs(r_q['Z'] - ref_que) / ref_que * 100
        mark_h = "*" if hyp_err < 5 else ""
        mark_q = "*" if que_err < 5 else ""
        print(f"{lam:<10.4f} | {r_h['Z']:<8.2f} | {hyp_err:<7.1f}%{mark_h} | {r_q['Z']:<8.2f} | {que_err:<7.1f}%{mark_q}")
    
    # 2. Percentile window sensitivity
    print("\n2. PERCENTILE WINDOW SENSITIVITY (lambda=0.0195):")
    print(f"{'Window':<10} | {'Hyp Z':<8} | {'Hyp Err':<8} | {'Que Z':<8} | {'Que Err':<8}")
    print("-" * 55)
    
    for pct in [0.05, 0.08, 0.10, 0.12, 0.15, 0.20]:
        config = ECNPConfig(percentile_window=pct)
        r_h = ecnp.compute(hyperforin, config)
        r_q = ecnp.compute(quercetin, config)
        hyp_err = abs(r_h['Z'] - ref_hyp) / ref_hyp * 100
        que_err = abs(r_q['Z'] - ref_que) / ref_que * 100
        print(f"{pct:<10.2f} | {r_h['Z']:<8.2f} | {hyp_err:<7.1f}% | {r_q['Z']:<8.2f} | {que_err:<7.1f}%")
    
    # 3. Find optimal lambda for Quercetin
    print("\n3. OPTIMAL LAMBDA SEARCH FOR QUERCETIN:")
    best_lam = None
    best_combined_err = float('inf')
    
    for lam in np.linspace(0.01, 0.10, 50):
        config = ECNPConfig(lambda_redundancy=lam)
        r_h = ecnp.compute(hyperforin, config)
        r_q = ecnp.compute(quercetin, config)
        hyp_err = abs(r_h['Z'] - ref_hyp) / ref_hyp * 100
        que_err = abs(r_q['Z'] - ref_que) / ref_que * 100
        combined = max(hyp_err, que_err)  # Minimize max error
        
        if combined < best_combined_err:
            best_combined_err = combined
            best_lam = lam
            best_hyp = (r_h['Z'], hyp_err)
            best_que = (r_q['Z'], que_err)
    
    print(f"Best lambda (minimizing max error): {best_lam:.4f}")
    print(f"  Hyperforin: Z={best_hyp[0]:.2f}, error={best_hyp[1]:.1f}%")
    print(f"  Quercetin: Z={best_que[0]:.2f}, error={best_que[1]:.1f}%")
    
    # 4. Diagnose Quercetin specifically
    print("\n4. QUERCETIN DIAGNOSTICS:")
    config = ECNPConfig()
    r = ecnp.compute(quercetin, config)
    
    target_idx = np.array([ecnp.node_to_idx[t] for t in quercetin if t in ecnp.node_to_idx])
    target_pcts = ecnp.percentiles_array[target_idx]
    target_degs = ecnp.degrees_array[target_idx]
    
    print(f"  k = {r['k']}")
    print(f"  Pool size = {r['pool_size']}")
    print(f"  Target percentiles: min={target_pcts.min():.2%}, max={target_pcts.max():.2%}, mean={target_pcts.mean():.2%}")
    print(f"  Target degrees: min={target_degs.min()}, max={target_degs.max()}, mean={target_degs.mean():.1f}")
    print(f"  I(T) = {r['I_T']:.4f}")
    print(f"  mu_T = {r['mu_T']:.4f}")
    print(f"  sigma_T = {r['sigma_T']:.4f}")
    print(f"  mean_rho = {r['mean_rho']:.4f}")
    print(f"  Z = {r['Z']:.2f}")
    print(f"  Gap = I(T) - mu_T = {r['I_T'] - r['mu_T']:.4f}")
    
    # 5. What Z would we need?
    print("\n5. REVERSE ENGINEERING:")
    target_Z = ref_que
    needed_sigma = (r['I_T'] - r['mu_T']) / target_Z
    print(f"  To get Z={target_Z}: need sigma={needed_sigma:.4f} (current: {r['sigma_T']:.4f})")
    print(f"  Ratio: {needed_sigma / r['sigma_T']:.2f}x")
    
    # What lambda gives that sigma?
    pool_var = r['sigma_T']**2 * r['k'] - 0.0195 * r['mean_rho']
    needed_sigma_sq = needed_sigma**2
    needed_lam = (needed_sigma_sq - pool_var / r['k']) / r['mean_rho']
    print(f"  Needed lambda for that sigma: {needed_lam:.4f}")

if __name__ == "__main__":
    main()

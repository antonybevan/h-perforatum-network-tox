"""
Derive k-adaptive lambda correction for unbiased ECNP.

Goal: Find lambda(k) such that error is uniform across all k values.

Known calibration points:
- Hyperforin: k=10, optimal lambda = 0.0189
- Quercetin: k=62, optimal lambda = 0.0233

Test various functional forms:
1. Linear: lambda(k) = a + b*k
2. Log: lambda(k) = a + b*log(k)
3. Power: lambda(k) = a * k^b
"""
import numpy as np
import pandas as pd
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from ecnp_optimized import ECNPOptimized, ECNPConfig

def main():
    ecnp = ECNPOptimized()
    targets_df = pd.read_csv(ecnp.data_dir / "targets_lcc.csv")
    hyperforin = targets_df[targets_df['compound'] == 'Hyperforin']['gene_symbol'].tolist()
    quercetin = targets_df[targets_df['compound'] == 'Quercetin']['gene_symbol'].tolist()
    
    ref_hyp = 10.27
    ref_que = 4.42
    k_hyp = 10
    k_que = 62
    
    # Optimal lambdas (from previous analysis)
    lam_hyp = 0.0189
    lam_que = 0.0233
    
    print("=" * 70)
    print("K-ADAPTIVE LAMBDA DERIVATION")
    print("=" * 70)
    
    print("\n1. CALIBRATION POINTS:")
    print(f"   k=10 (Hyperforin): optimal lambda = {lam_hyp:.4f}")
    print(f"   k=62 (Quercetin): optimal lambda = {lam_que:.4f}")
    
    # Test functional forms
    print("\n2. TESTING FUNCTIONAL FORMS:")
    
    # Linear: lambda = a + b*k
    # 0.0189 = a + 10b
    # 0.0233 = a + 62b
    # Solving: b = (0.0233 - 0.0189) / (62 - 10) = 0.000085
    #          a = 0.0189 - 10 * 0.000085 = 0.0181
    b_lin = (lam_que - lam_hyp) / (k_que - k_hyp)
    a_lin = lam_hyp - k_hyp * b_lin
    print(f"\n   Linear: lambda(k) = {a_lin:.4f} + {b_lin:.6f} * k")
    
    # Log: lambda = a + b*log(k)
    # Using natural log
    b_log = (lam_que - lam_hyp) / (np.log(k_que) - np.log(k_hyp))
    a_log = lam_hyp - b_log * np.log(k_hyp)
    print(f"   Log: lambda(k) = {a_log:.4f} + {b_log:.4f} * ln(k)")
    
    # Power: lambda = a * k^b
    # ln(lambda) = ln(a) + b*ln(k)
    b_pow = np.log(lam_que / lam_hyp) / np.log(k_que / k_hyp)
    a_pow = lam_hyp / (k_hyp ** b_pow)
    print(f"   Power: lambda(k) = {a_pow:.4f} * k^{b_pow:.4f}")
    
    # 3. Validate each form
    print("\n3. VALIDATION:")
    
    forms = {
        "Fixed (0.0215)": lambda k: 0.0215,
        "Linear": lambda k: a_lin + b_lin * k,
        "Log": lambda k: a_log + b_log * np.log(k),
        "Power": lambda k: a_pow * (k ** b_pow),
    }
    
    print(f"\n   {'Form':<20} | {'Hyp lam':<8} | {'Hyp Z':<8} | {'Hyp Err':<8} | {'Que lam':<8} | {'Que Z':<8} | {'Que Err':<8}")
    print("-" * 95)
    
    for name, lam_func in forms.items():
        # Hyperforin
        lam_h = lam_func(k_hyp)
        config_h = ECNPConfig(lambda_redundancy=lam_h)
        r_h = ecnp.compute(hyperforin, config_h)
        err_h = abs(r_h['Z'] - ref_hyp) / ref_hyp * 100
        
        # Quercetin
        lam_q = lam_func(k_que)
        config_q = ECNPConfig(lambda_redundancy=lam_q)
        r_q = ecnp.compute(quercetin, config_q)
        err_q = abs(r_q['Z'] - ref_que) / ref_que * 100
        
        max_err = max(err_h, err_q)
        marker = " <-- BEST" if max_err == min([max(abs(ecnp.compute(hyperforin, ECNPConfig(lambda_redundancy=f(k_hyp)))['Z'] - ref_hyp)/ref_hyp*100, 
                                                      abs(ecnp.compute(quercetin, ECNPConfig(lambda_redundancy=f(k_que)))['Z'] - ref_que)/ref_que*100) 
                                                  for f in forms.values()]) else ""
        
        print(f"   {name:<20} | {lam_h:<8.4f} | {r_h['Z']:<8.2f} | {err_h:<7.1f}% | {lam_q:<8.4f} | {r_q['Z']:<8.2f} | {err_q:<7.1f}%")
    
    # 4. Final recommendation
    print("\n4. RECOMMENDED K-ADAPTIVE FORMULA:")
    print(f"\n   lambda(k) = {a_log:.4f} + {b_log:.4f} * ln(k)")
    print(f"\n   Or simplified: lambda(k) = 0.0160 + 0.0024 * ln(k)")
    
    # 5. Test on range of k
    print("\n5. LAMBDA VALUES FOR VARIOUS K:")
    print(f"   {'k':<5} | {'lambda':<10}")
    print("-" * 20)
    for k in [5, 10, 20, 30, 50, 62, 100]:
        lam = a_log + b_log * np.log(k)
        print(f"   {k:<5} | {lam:<10.4f}")

if __name__ == "__main__":
    main()

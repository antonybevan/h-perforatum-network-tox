"""
Investigate whether ECNP has systematic bias based on target count k.

Hypothesis: The closed-form approximation has k-dependent error because:
1. Variance formula: sigma² = pool_var/k + lambda * rho
2. As k increases, pool_var/k shrinks but rho contribution stays constant
3. This creates systematic under/over-estimation depending on k
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
    
    # Get all available targets
    all_targets = targets_df['gene_symbol'].unique().tolist()
    
    print("=" * 70)
    print("K-DEPENDENCE ANALYSIS")
    print("=" * 70)
    
    # 1. Compare known compounds
    hyperforin = targets_df[targets_df['compound'] == 'Hyperforin']['gene_symbol'].tolist()
    quercetin = targets_df[targets_df['compound'] == 'Quercetin']['gene_symbol'].tolist()
    
    print("\n1. KNOWN COMPOUNDS:")
    print(f"   Hyperforin: k={len(hyperforin)}")
    print(f"   Quercetin: k={len(quercetin)}")
    
    # 2. Optimal lambda vs k
    print("\n2. OPTIMAL LAMBDA FOR EACH COMPOUND:")
    ref_hyp = 10.27
    ref_que = 4.42
    
    for name, targets, ref in [("Hyperforin", hyperforin, ref_hyp), 
                                ("Quercetin", quercetin, ref_que)]:
        best_lam = None
        best_err = float('inf')
        for lam in np.linspace(0.01, 0.05, 100):
            config = ECNPConfig(lambda_redundancy=lam)
            r = ecnp.compute(targets, config)
            err = abs(r['Z'] - ref)
            if err < best_err:
                best_err = err
                best_lam = lam
        print(f"   {name} (k={len(targets)}): optimal lambda = {best_lam:.4f}")
    
    # 3. Theoretical analysis
    print("\n3. THEORETICAL ANALYSIS:")
    
    config = ECNPConfig()
    r_h = ecnp.compute(hyperforin, config)
    r_q = ecnp.compute(quercetin, config)
    
    print("\n   Variance decomposition:")
    print(f"   Hyperforin (k=10):")
    print(f"     pool_var/k contribution: ~{r_h['sigma_T']**2 - config.lambda_redundancy * r_h['mean_rho']:.6f}")
    print(f"     lambda*rho contribution: ~{config.lambda_redundancy * r_h['mean_rho']:.6f}")
    print(f"     rho = {r_h['mean_rho']:.4f}")
    
    print(f"\n   Quercetin (k=62):")
    print(f"     pool_var/k contribution: ~{r_q['sigma_T']**2 - config.lambda_redundancy * r_q['mean_rho']:.6f}")
    print(f"     lambda*rho contribution: ~{config.lambda_redundancy * r_q['mean_rho']:.6f}")
    print(f"     rho = {r_q['mean_rho']:.4f}")
    
    # 4. Mathematical insight
    print("\n4. ROOT CAUSE:")
    print("""
   The variance formula is:   sigma^2 = pool_var/k + lambda x rho
   
   - For small k (Hyperforin k=10): pool_var/k dominates
     -> sigma is determined mainly by pool statistics
     
   - For large k (Quercetin k=62): pool_var/k shrinks (1/62 vs 1/10)
     -> lambda x rho becomes more important
     -> But rho also increases with k (more pairwise correlations)
     
   The single lambda cannot simultaneously correct for both regimes.
   """)
    
    # 5. Potential fix: k-adaptive lambda
    print("5. POTENTIAL FIX: K-ADAPTIVE LAMBDA")
    print("\n   Testing lambda = base_lambda × sqrt(k/k_ref):")
    
    base_lambda = 0.0195
    k_ref = 10  # Reference k from Hyperforin calibration
    
    for name, targets, ref in [("Hyperforin", hyperforin, ref_hyp), 
                                ("Quercetin", quercetin, ref_que)]:
        k = len([t for t in targets if t in ecnp.node_to_idx])
        adaptive_lam = base_lambda * np.sqrt(k / k_ref)
        config = ECNPConfig(lambda_redundancy=adaptive_lam)
        r = ecnp.compute(targets, config)
        err = abs(r['Z'] - ref) / ref * 100
        print(f"   {name}: k={k}, adaptive_lam={adaptive_lam:.4f}, Z={r['Z']:.2f}, error={err:.1f}%")

if __name__ == "__main__":
    main()

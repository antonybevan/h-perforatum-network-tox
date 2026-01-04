"""
Biological Realism Check for ECNP Algorithm

Core question: Does the closed-form algorithm preserve the biological hierarchy 
established in the manuscript?

Expected: Hyperforin > Quercetin by ~2.3x in Z-score
Expected: Hyperforin targets are high-leverage hubs (PXR, CYPs)
Expected: PTNI ratio ~17-22x
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from ecnp_optimized import ECNPOptimized, ECNPConfig
import pandas as pd
import numpy as np

def main():
    ecnp = ECNPOptimized()
    targets_df = pd.read_csv(ecnp.data_dir / "targets_lcc.csv")
    hyperforin = targets_df[targets_df['compound'] == 'Hyperforin']['gene_symbol'].tolist()
    quercetin = targets_df[targets_df['compound'] == 'Quercetin']['gene_symbol'].tolist()

    print("=" * 70)
    print("BIOLOGICAL REALISM CHECK")
    print("=" * 70)

    print("\n1. BASIC METRICS:")
    r_h = ecnp.compute(hyperforin)
    r_q = ecnp.compute(quercetin)

    print(f"   Hyperforin: k={r_h['k']}, Z={r_h['Z']:.2f}")
    print(f"   Quercetin: k={r_q['k']}, Z={r_q['Z']:.2f}")
    z_ratio = r_h['Z'] / r_q['Z']
    print(f"   Ranking: Hyperforin > Quercetin by {z_ratio:.1f}x")

    print("\n2. PER-TARGET INFLUENCE (PTNI):")
    ptni_h = r_h['I_T'] / r_h['k']
    ptni_q = r_q['I_T'] / r_q['k']
    print(f"   Hyperforin PTNI: {ptni_h:.4f}")
    print(f"   Quercetin PTNI: {ptni_q:.4f}")
    ptni_ratio = ptni_h / ptni_q
    print(f"   Fold difference: {ptni_ratio:.1f}x")

    print("\n3. HYPERFORIN TARGET ANALYSIS:")
    print("   (These should be regulatory hubs: PXR, CYPs, etc.)")
    for t in hyperforin:
        if t in ecnp.node_to_idx:
            idx = ecnp.node_to_idx[t]
            pct = ecnp.percentiles_array[idx]
            deg = ecnp.degrees_array[idx]
            m = ecnp.m_array[idx]
            print(f"   {t:12s}: degree={deg:3.0f}, influence_pct={pct:5.1%}, m={m:.4f}")

    print("\n4. TOP 10 QUERCETIN TARGETS BY INFLUENCE:")
    que_data = []
    for t in quercetin:
        if t in ecnp.node_to_idx:
            idx = ecnp.node_to_idx[t]
            que_data.append((t, ecnp.m_array[idx], ecnp.degrees_array[idx], ecnp.percentiles_array[idx]))
    que_sorted = sorted(que_data, key=lambda x: -x[1])[:10]
    for t, m, deg, pct in que_sorted:
        print(f"   {t:12s}: degree={deg:3.0f}, influence_pct={pct:5.1%}, m={m:.4f}")

    print("\n5. BIOLOGICAL HIERARCHY VALIDATION:")
    print("-" * 50)
    
    checks = [
        ("Hyperforin Z > Quercetin Z", r_h['Z'] > r_q['Z']),
        ("Z ratio in expected range (1.5-3.0x)", 1.5 <= z_ratio <= 3.0),
        ("PTNI ratio > 10x", ptni_ratio > 10),
        ("Hyperforin targets in high percentiles", 
         np.mean([ecnp.percentiles_array[ecnp.node_to_idx[t]] for t in hyperforin if t in ecnp.node_to_idx]) > 0.7),
    ]
    
    all_pass = True
    for name, result in checks:
        status = "PASS" if result else "FAIL"
        if not result:
            all_pass = False
        print(f"   [{status}] {name}")
    
    print("\n6. MANUSCRIPT COMPARISON:")
    print("-" * 50)
    print(f"   Manuscript RWR Z-scores: Hyp=+10.3, Que=+4.4, ratio=2.3x")
    print(f"   ECNP closed-form:        Hyp={r_h['Z']:+.1f}, Que={r_q['Z']:+.1f}, ratio={z_ratio:.1f}x")
    
    mc_hyp = 10.27
    mc_que = 4.42
    hyp_err = abs(r_h['Z'] - mc_hyp) / mc_hyp * 100
    que_err = abs(r_q['Z'] - mc_que) / mc_que * 100
    print(f"   Error vs MC: Hyp={hyp_err:.1f}%, Que={que_err:.1f}%")
    
    print("\n" + "=" * 70)
    if all_pass and hyp_err < 5 and que_err < 5:
        print("BIOLOGICAL REALISM: PRESERVED")
        print("The closed-form algorithm correctly identifies Hyperforin as the")
        print("higher-leverage perturbation despite having fewer targets.")
    else:
        print("BIOLOGICAL REALISM: ISSUES DETECTED")
        print("Review the failed checks above.")
    print("=" * 70)

if __name__ == "__main__":
    main()

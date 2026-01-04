"""
CRITICAL SANITY CHECK: Does ECNP Add Value Over Simple Counting?

This script tests the null hypothesis:
H0: ECNP score ≈ direct DILI hit count (RWR is decorative)
H1: ECNP captures additional signal beyond direct hits

Tests:
1. Correlation between ECNP Z-score and direct hit count
2. ECNP score for targets with ZERO direct hits (should still rank if RWR matters)
3. Ablation: Remove direct hits, check if ranking preserved

If ECNP ≈ counting, the RWR machinery is meaningless.

Author: Research Session 2026-01-02
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path

# Add paths
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[4]  # validation -> src -> ecnp-closed-form -> research -> h-perforatum-network-tox
sys.path.insert(0, str(SCRIPT_DIR.parent / "core"))

from ecnp_optimized import ECNPOptimized

PROJECT_ROOT = Path(__file__).resolve().parents[4]
DATA_DIR = PROJECT_ROOT / "data" / "processed"


def load_dili_genes():
    """Load DILI gene set."""
    dili_df = pd.read_csv(DATA_DIR / "dili_900_lcc.csv")
    return set(dili_df['gene_name'])


def analyze_compound(ecnp, targets, dili_genes, compound_name):
    """Analyze a single compound."""
    # ECNP score
    result = ecnp.compute(targets)
    
    # Direct hit analysis
    targets_in_network = [t for t in targets if t in ecnp.node_to_idx]
    direct_hits = [t for t in targets_in_network if t in dili_genes]
    
    # Neighbor analysis (1-hop from DILI)
    dili_indices = {ecnp.node_to_idx[g] for g in dili_genes if g in ecnp.node_to_idx}
    
    # Per-target influence breakdown
    target_influences = []
    for t in targets_in_network:
        idx = ecnp.node_to_idx[t]
        m = ecnp.m_array[idx]
        is_dili = t in dili_genes
        target_influences.append({
            'target': t,
            'influence': m,
            'is_dili': is_dili,
            'percentile': ecnp.percentiles_array[idx]
        })
    
    return {
        'compound': compound_name,
        'n_targets': len(targets),
        'n_targets_in_network': len(targets_in_network),
        'n_direct_hits': len(direct_hits),
        'direct_hit_fraction': len(direct_hits) / len(targets_in_network) if targets_in_network else 0,
        'I_T': result['I_T'],
        'Z': result['Z'],
        'target_details': target_influences
    }


def main():
    print("=" * 80)
    print("ECNP SANITY CHECK: Does RWR Add Value Over Simple Counting?")
    print("=" * 80)
    
    # Load data
    ecnp = ECNPOptimized(project_root=PROJECT_ROOT)
    dili_genes = load_dili_genes()
    targets_df = pd.read_csv(DATA_DIR / "targets_lcc.csv")
    
    # Get all compounds
    compounds = targets_df['compound'].unique()
    
    print(f"\nLoaded {len(compounds)} compounds")
    print(f"DILI genes in network: {len(dili_genes)}")
    
    # Analyze each compound
    results = []
    for compound in compounds:
        targets = targets_df[targets_df['compound'] == compound]['gene_symbol'].tolist()
        if len(targets) >= 2:  # Minimum for ECNP
            analysis = analyze_compound(ecnp, targets, dili_genes, compound)
            results.append(analysis)
    
    results_df = pd.DataFrame([{k: v for k, v in r.items() if k != 'target_details'} for r in results])
    
    print("\n" + "-" * 80)
    print("TEST 1: Correlation between ECNP Z-score and Direct Hit Count")
    print("-" * 80)
    
    correlation = results_df['Z'].corr(results_df['n_direct_hits'])
    correlation_frac = results_df['Z'].corr(results_df['direct_hit_fraction'])
    
    print(f"\n  Correlation(Z, n_direct_hits) = {correlation:.4f}")
    print(f"  Correlation(Z, direct_hit_fraction) = {correlation_frac:.4f}")
    
    if correlation > 0.90:
        print("\n  ⚠️  HIGH CORRELATION: ECNP may just be counting direct hits!")
    elif correlation > 0.70:
        print("\n  ⚠️  MODERATE CORRELATION: Direct hits explain significant variance")
    else:
        print("\n  ✅ LOW CORRELATION: ECNP captures more than direct hits")
    
    print("\n" + "-" * 80)
    print("TEST 2: Compounds with ZERO Direct Hits")
    print("-" * 80)
    
    zero_hit = results_df[results_df['n_direct_hits'] == 0]
    has_hits = results_df[results_df['n_direct_hits'] > 0]
    
    print(f"\n  Compounds with 0 direct DILI hits: {len(zero_hit)}")
    print(f"  Compounds with 1+ direct DILI hits: {len(has_hits)}")
    
    if len(zero_hit) > 0:
        print(f"\n  Zero-hit compounds still have Z-scores:")
        for _, row in zero_hit.iterrows():
            print(f"    {row['compound']}: Z = {row['Z']:.2f}, k = {row['n_targets']}")
        
        zero_hit_z_mean = zero_hit['Z'].mean()
        has_hit_z_mean = has_hits['Z'].mean()
        
        print(f"\n  Mean Z (0 hits): {zero_hit_z_mean:.2f}")
        print(f"  Mean Z (1+ hits): {has_hit_z_mean:.2f}")
        
        if zero_hit_z_mean < has_hit_z_mean:
            print("  ✅ RWR does differentiate: compounds with direct hits score higher")
        else:
            print("  ⚠️  No differentiation: direct hits don't predict Z-score")
    
    print("\n" + "-" * 80)
    print("TEST 3: Hyperforin/Quercetin Deep Dive")
    print("-" * 80)
    
    for compound in ['Hyperforin', 'Quercetin']:
        analysis = next((r for r in results if r['compound'] == compound), None)
        if analysis:
            print(f"\n  {compound}:")
            print(f"    Targets: {analysis['n_targets_in_network']}")
            print(f"    Direct DILI hits: {analysis['n_direct_hits']} ({analysis['direct_hit_fraction']:.1%})")
            print(f"    ECNP Z-score: {analysis['Z']:.2f}")
            print(f"    Total influence I(T): {analysis['I_T']:.4f}")
            
            # Breakdown by direct vs non-direct
            details = analysis['target_details']
            direct = [d for d in details if d['is_dili']]
            non_direct = [d for d in details if not d['is_dili']]
            
            direct_influence = sum(d['influence'] for d in direct)
            non_direct_influence = sum(d['influence'] for d in non_direct)
            total_influence = direct_influence + non_direct_influence
            
            print(f"\n    Influence from DIRECT DILI targets: {direct_influence:.4f} ({direct_influence/total_influence*100:.1f}%)")
            print(f"    Influence from NON-DILI targets: {non_direct_influence:.4f} ({non_direct_influence/total_influence*100:.1f}%)")
            
            print(f"\n    Top 5 targets by influence:")
            sorted_details = sorted(details, key=lambda x: -x['influence'])[:5]
            for d in sorted_details:
                dili_marker = "★ DILI" if d['is_dili'] else ""
                print(f"      {d['target']:12s}: m = {d['influence']:.4f} (pct: {d['percentile']:.1%}) {dili_marker}")
    
    print("\n" + "-" * 80)
    print("TEST 4: Remove Direct Hits - Does RWR Still Work?")
    print("-" * 80)
    
    # Ablation: score compounds without their direct DILI hits
    for compound in ['Hyperforin', 'Quercetin']:
        targets = targets_df[targets_df['compound'] == compound]['gene_symbol'].tolist()
        targets_in_network = [t for t in targets if t in ecnp.node_to_idx]
        
        # Original score
        original_result = ecnp.compute(targets)
        
        # Remove direct hits
        non_dili_targets = [t for t in targets_in_network if t not in dili_genes]
        
        if len(non_dili_targets) >= 2:
            ablated_result = ecnp.compute(non_dili_targets)
            
            print(f"\n  {compound}:")
            print(f"    Original: k={len(targets_in_network)}, Z={original_result['Z']:.2f}")
            print(f"    Without direct hits: k={len(non_dili_targets)}, Z={ablated_result['Z']:.2f}")
            
            if ablated_result['Z'] > 2:
                print(f"    ✅ RWR still detects signal without direct hits")
            elif ablated_result['Z'] > 0:
                print(f"    ⚠️  Weak signal remains")
            else:
                print(f"    ❌ No signal without direct hits")
        else:
            print(f"\n  {compound}: Not enough non-DILI targets for ablation")
    
    print("\n" + "-" * 80)
    print("TEST 5: Theoretical Grounding Check")
    print("-" * 80)
    
    print("""
  Q1: What does RWR influence matrix actually measure?
      A: M[i,j] = probability of reaching node i starting from j with random walk
         This captures network topology, not just adjacency
  
  Q2: Why might ECNP work beyond counting?
      A: A target that is NOT a DILI gene but is a HUB connected to many DILI genes
         will have high m_j even though direct_hit = 0
  
  Q3: What would prove RWR adds value?
      A: Compounds with zero direct hits but high Z-scores (network effects)
         Or: Removing direct hits doesn't collapse the ranking
    """)
    
    print("\n" + "=" * 80)
    print("SUMMARY: ECNP Validity Assessment")
    print("=" * 80)
    
    # Final assessment
    validity_checks = []
    
    # Check 1: Correlation
    if correlation < 0.70:
        validity_checks.append(("Correlation test", "PASS", "RWR captures more than counting"))
    elif correlation < 0.90:
        validity_checks.append(("Correlation test", "WARN", f"Moderate correlation ({correlation:.2f})"))
    else:
        validity_checks.append(("Correlation test", "FAIL", f"High correlation ({correlation:.2f}) - may be just counting"))
    
    # Check 2: Zero-hit differentiation
    if len(zero_hit) > 0 and len(has_hits) > 0:
        if has_hits['Z'].mean() > zero_hit['Z'].mean():
            validity_checks.append(("Zero-hit differentiation", "PASS", "Direct hits predict higher Z"))
        else:
            validity_checks.append(("Zero-hit differentiation", "FAIL", "No differentiation"))
    else:
        validity_checks.append(("Zero-hit differentiation", "N/A", "Insufficient data"))
    
    # Check 3: Hyperforin non-DILI influence
    hyp_analysis = next((r for r in results if r['compound'] == 'Hyperforin'), None)
    if hyp_analysis:
        details = hyp_analysis['target_details']
        non_dili = [d for d in details if not d['is_dili']]
        non_dili_influence = sum(d['influence'] for d in non_dili)
        total = sum(d['influence'] for d in details)
        non_dili_frac = non_dili_influence / total if total > 0 else 0
        
        if non_dili_frac > 0.20:
            validity_checks.append(("Non-DILI contribution", "PASS", f"{non_dili_frac:.0%} from non-DILI targets"))
        else:
            validity_checks.append(("Non-DILI contribution", "WARN", f"Only {non_dili_frac:.0%} from non-DILI"))
    
    print("\n")
    for name, status, msg in validity_checks:
        icon = "✅" if status == "PASS" else ("⚠️" if status == "WARN" else "❌")
        print(f"  {icon} {name}: {status} - {msg}")
    
    n_pass = sum(1 for _, s, _ in validity_checks if s == "PASS")
    n_total = len(validity_checks)
    
    print(f"\n  Overall: {n_pass}/{n_total} checks passed")
    
    if n_pass == n_total:
        print("\n  🟢 ECNP APPEARS VALID: RWR adds value beyond simple counting")
    elif n_pass >= n_total // 2:
        print("\n  🟡 ECNP PARTIALLY VALID: Some signal from RWR, but direct hits dominate")
    else:
        print("\n  🔴 ECNP QUESTIONABLE: May be equivalent to counting direct hits")
    
    print("\n" + "=" * 80)


if __name__ == "__main__":
    main()

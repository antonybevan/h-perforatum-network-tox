"""
Diagnostic: Why 706-Drug Model Underperforms 202-Drug Model
============================================================

Investigating the AUC gap: 0.62 vs 0.87
"""
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score
import warnings
warnings.filterwarnings('ignore')

ROOT = Path(r'v:\new\h-perforatum-network-tox')

print("="*70)
print("DIAGNOSTIC: 706 vs 202 DRUG MODEL PERFORMANCE GAP")
print("="*70)

# =============================================================================
# 1. LOAD BOTH DATASETS
# =============================================================================
print("\n--- Loading Datasets ---")

# 706-drug dataset
df_706 = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'results' / 'dilirank_706_with_ecnp.csv')
print(f"\n706-Drug Dataset:")
print(f"  Total: {len(df_706)}")
print(f"  Valid ECNP: {df_706['ecnp_z'].notna().sum()}")

# Original 202-drug dataset
try:
    df_202 = pd.read_csv(ROOT / 'data' / 'processed' / 'dili_with_ecnp.csv')
    print(f"\n202-Drug Dataset:")
    print(f"  Total: {len(df_202)}")
    has_202 = True
except:
    # Try alternative locations
    try:
        df_202 = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'data' / 'curated' / 'labeled_compounds_improved.csv')
        print(f"\n202-Drug Dataset (from curated):")
        print(f"  Total: {len(df_202)}")
        has_202 = True
    except:
        print("\n202-Drug dataset not found")
        has_202 = False

# =============================================================================
# 2. CLASS BALANCE COMPARISON
# =============================================================================
print("\n" + "="*70)
print("1. CLASS BALANCE")
print("="*70)

df_706_valid = df_706[df_706['ecnp_z'].notna()]
dili_rate_706 = df_706_valid['is_dili'].mean()

print(f"\n706-Drug Dataset (valid ECNP):")
print(f"  DILI+: {df_706_valid['is_dili'].sum()} ({dili_rate_706*100:.1f}%)")
print(f"  DILI-: {(df_706_valid['is_dili']==0).sum()} ({(1-dili_rate_706)*100:.1f}%)")

if has_202 and ('is_dili' in df_202.columns or 'label' in df_202.columns):
    label_col = 'is_dili' if 'is_dili' in df_202.columns else 'label'
    dili_rate_202 = df_202[label_col].mean()
    print(f"\n202-Drug Dataset:")
    print(f"  DILI+: {df_202[label_col].sum()} ({dili_rate_202*100:.1f}%)")
    print(f"  DILI-: {(df_202[label_col]==0).sum()} ({(1-dili_rate_202)*100:.1f}%)")

# =============================================================================
# 3. ECNP COVERAGE & QUALITY
# =============================================================================
print("\n" + "="*70)
print("2. ECNP COVERAGE & QUALITY")
print("="*70)

# Target mapping stats for 706
print(f"\n706-Drug Dataset:")
df_706_valid = df_706[df_706['ecnp_z'].notna()]
print(f"  Mean targets per drug: {df_706_valid['n_targets_mapped'].mean():.1f}")
print(f"  Mean pool size: {df_706_valid['pool_size'].mean():.1f}")
print(f"  Drugs with pool_size > 0: {(df_706_valid['pool_size'] > 0).sum()}")

# ECNP Z distribution
print(f"\n  ECNP Z-score distribution:")
print(f"    Mean: {df_706_valid['ecnp_z'].mean():.3f}")
print(f"    Std:  {df_706_valid['ecnp_z'].std():.3f}")
print(f"    Min:  {df_706_valid['ecnp_z'].min():.3f}")
print(f"    Max:  {df_706_valid['ecnp_z'].max():.3f}")

# =============================================================================
# 4. ECNP DISCRIMINATIVE POWER BY SUBGROUP
# =============================================================================
print("\n" + "="*70)
print("3. ECNP DISCRIMINATIVE POWER BY SUBGROUP")
print("="*70)

# By target count
print("\nBy number of mapped targets:")
for k_min, k_max, label in [(1,2,'k=1-2'), (3,5,'k=3-5'), (6,10,'k=6-10'), (11,100,'k>10')]:
    subset = df_706_valid[(df_706_valid['n_targets_mapped'] >= k_min) & (df_706_valid['n_targets_mapped'] <= k_max)]
    if len(subset) >= 10:
        n_pos = subset['is_dili'].sum()
        n_neg = (subset['is_dili']==0).sum()
        if n_pos >= 3 and n_neg >= 3:
            auc = roc_auc_score(subset['is_dili'], subset['ecnp_z'])
            print(f"  {label}: n={len(subset)}, {n_pos}:{n_neg}, AUC = {auc:.3f}")
        else:
            print(f"  {label}: n={len(subset)}, {n_pos}:{n_neg}, (imbalanced)")

# By pool size
print("\nBy pool size (network coverage):")
for p_min, p_max, label in [(1,100,'pool=1-100'), (101,300,'pool=101-300'), (301,1000,'pool=301-1000'), (1001,5000,'pool>1000')]:
    subset = df_706_valid[(df_706_valid['pool_size'] >= p_min) & (df_706_valid['pool_size'] <= p_max)]
    if len(subset) >= 10:
        n_pos = subset['is_dili'].sum()
        n_neg = (subset['is_dili']==0).sum()
        if n_pos >= 3 and n_neg >= 3:
            auc = roc_auc_score(subset['is_dili'], subset['ecnp_z'])
            print(f"  {label}: n={len(subset)}, {n_pos}:{n_neg}, AUC = {auc:.3f}")
        else:
            print(f"  {label}: n={len(subset)}, {n_pos}:{n_neg}, (imbalanced)")

# =============================================================================
# 5. OVERLAP ANALYSIS - Are 202 drugs IN the 706?
# =============================================================================
print("\n" + "="*70)
print("4. OVERLAP ANALYSIS")
print("="*70)

if has_202:
    # Check for drugbank_id overlap
    if 'drugbank_id' in df_202.columns and 'drugbank_id' in df_706.columns:
        overlap_ids = set(df_202['drugbank_id'].dropna()) & set(df_706['drugbank_id'].dropna())
        print(f"\nDrugBank ID overlap: {len(overlap_ids)} drugs")
        
        if len(overlap_ids) > 0:
            # Performance on overlapping subset
            df_overlap = df_706[df_706['drugbank_id'].isin(overlap_ids) & df_706['ecnp_z'].notna()]
            if len(df_overlap) >= 10:
                n_pos = df_overlap['is_dili'].sum()
                n_neg = (df_overlap['is_dili']==0).sum()
                if n_pos >= 3 and n_neg >= 3:
                    auc_overlap = roc_auc_score(df_overlap['is_dili'], df_overlap['ecnp_z'])
                    print(f"ECNP AUC on overlapping drugs: {auc_overlap:.3f}")
                    print(f"  (n={len(df_overlap)}, {n_pos}:{n_neg})")
    else:
        print("\nNo drugbank_id column for overlap analysis")

# =============================================================================
# 6. DIAGNOSIS: DILI+ vs DILI- ECNP DISTRIBUTIONS
# =============================================================================
print("\n" + "="*70)
print("5. DILI+ vs DILI- SEPARATION")
print("="*70)

dili_pos = df_706_valid[df_706_valid['is_dili']==1]['ecnp_z']
dili_neg = df_706_valid[df_706_valid['is_dili']==0]['ecnp_z']

print(f"\nDILI+ drugs (n={len(dili_pos)}):")
print(f"  Mean ECNP Z: {dili_pos.mean():.3f}")
print(f"  Median:      {dili_pos.median():.3f}")

print(f"\nDILI- drugs (n={len(dili_neg)}):")
print(f"  Mean ECNP Z: {dili_neg.mean():.3f}")
print(f"  Median:      {dili_neg.median():.3f}")

diff = dili_pos.mean() - dili_neg.mean()
print(f"\nSeparation (DILI+ - DILI-): {diff:+.3f}")

if diff > 0.1:
    print("  → ECNP is separating classes (DILI+ > DILI-)")
elif diff < -0.1:
    print("  → ECNP is INVERSELY separating (DILI- > DILI+) - PROBLEM!")
else:
    print("  → Minimal separation - ECNP not discriminative")

# =============================================================================
# 7. SUMMARY & HYPOTHESIS
# =============================================================================
print("\n" + "="*70)
print("SUMMARY: LIKELY CAUSES OF PERFORMANCE GAP")
print("="*70)

print("""
1. CLASS IMBALANCE
   - 706-drug dataset: ~67% DILI+ (highly imbalanced)
   - Harder classification task with majority class dominating

2. ECNP COVERAGE QUALITY
   - Only 48% network overlap (vs likely higher in curated 202)
   - Gene symbol mapping introduces noise

3. ECNP DISCRIMINATIVE POWER
   - Check if DILI+ has higher ECNP than DILI- (should be positive)
   - Current separation is minimal

4. DATASET COMPOSITION
   - 706 drugs from DILIrank (diverse sources)
   - 202 drugs were likely curated for ECNP suitability

NEXT STEPS:
- Filter to drugs with higher network coverage (pool_size > 100)
- Stratify by mechanism if available
- Compare only overlapping drugs between datasets
""")

"""
Deep Accuracy Gap Analysis
===========================

Systematic analysis of where our model fails:
1. Error distribution by compound properties
2. Confusion matrix breakdown
3. Specific drug failure analysis
4. Regime-specific error patterns
5. Feature blind spots
"""
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score, confusion_matrix, classification_report
import warnings
warnings.filterwarnings('ignore')

ROOT = Path(r'v:\new\h-perforatum-network-tox')

# Load best model results
df = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'results' / 'ecfp_model_results.csv')
print(f"Compounds: {len(df)}")
print(f"DILI+: {df['is_dili'].sum()}, DILI-: {(df['is_dili']==0).sum()}")

# Use the combined prediction
if 'y_pred_combined' in df.columns:
    df['y_pred'] = df['y_pred_combined']
elif 'y_pred_ecfp' in df.columns:
    df['y_pred'] = df['y_pred_ecfp']
else:
    print("No prediction column found!")
    exit()

# Binary predictions at 0.5 threshold
df['y_pred_binary'] = (df['y_pred'] > 0.5).astype(int)
df['correct'] = df['y_pred_binary'] == df['is_dili']
df['error'] = ~df['correct']

# Error types
df['fp'] = (df['y_pred_binary'] == 1) & (df['is_dili'] == 0)  # False positive
df['fn'] = (df['y_pred_binary'] == 0) & (df['is_dili'] == 1)  # False negative
df['tp'] = (df['y_pred_binary'] == 1) & (df['is_dili'] == 1)  # True positive
df['tn'] = (df['y_pred_binary'] == 0) & (df['is_dili'] == 0)  # True negative

# =============================================================================
# 1. OVERALL ERROR STATISTICS
# =============================================================================

print("\n" + "="*70)
print("1. OVERALL ERROR STATISTICS")
print("="*70)

print(f"\nConfusion Matrix:")
print(f"  TP (correct DILI+): {df['tp'].sum()}")
print(f"  TN (correct DILI-): {df['tn'].sum()}")
print(f"  FP (false alarm):   {df['fp'].sum()}")
print(f"  FN (missed DILI):   {df['fn'].sum()}")

accuracy = df['correct'].mean()
precision = df['tp'].sum() / (df['tp'].sum() + df['fp'].sum())
recall = df['tp'].sum() / (df['tp'].sum() + df['fn'].sum())
specificity = df['tn'].sum() / (df['tn'].sum() + df['fp'].sum())

print(f"\nMetrics:")
print(f"  Accuracy:    {accuracy:.3f}")
print(f"  Precision:   {precision:.3f}")
print(f"  Recall:      {recall:.3f}")
print(f"  Specificity: {specificity:.3f}")

# =============================================================================
# 2. ERROR BY COMPOUND PROPERTIES
# =============================================================================

print("\n" + "="*70)
print("2. ERROR DISTRIBUTION BY PROPERTIES")
print("="*70)

# By LogP
if 'logp' in df.columns:
    print("\nError rate by LogP:")
    for low, high, label in [(0, 2, 'Low (0-2)'), (2, 4, 'Medium (2-4)'), (4, 10, 'High (>4)')]:
        subset = df[(df['logp'] >= low) & (df['logp'] < high)]
        if len(subset) > 0:
            err = subset['error'].mean() * 100
            print(f"  {label}: {subset['error'].sum()}/{len(subset)} ({err:.1f}%)")

# By target count
if 'k' in df.columns or 'n_targets' in df.columns:
    k_col = 'k' if 'k' in df.columns else 'n_targets'
    print(f"\nError rate by target count ({k_col}):")
    for low, high, label in [(0, 4, '1-3'), (4, 8, '4-7'), (8, 50, '8+')]:
        subset = df[(df[k_col] >= low) & (df[k_col] < high)]
        if len(subset) > 0:
            err = subset['error'].mean() * 100
            print(f"  {label}: {subset['error'].sum()}/{len(subset)} ({err:.1f}%)")

# By ECNP Z-score
if 'ecnp_z' in df.columns:
    print("\nError rate by ECNP Z-score:")
    for low, high, label in [(-5, 0, 'Negative'), (0, 1, 'Low (0-1)'), (1, 5, 'High (>1)')]:
        subset = df[(df['ecnp_z'] >= low) & (df['ecnp_z'] < high)]
        if len(subset) > 0:
            err = subset['error'].mean() * 100
            print(f"  {label}: {subset['error'].sum()}/{len(subset)} ({err:.1f}%)")

# =============================================================================
# 3. SPECIFIC FAILURE CASES
# =============================================================================

print("\n" + "="*70)
print("3. SPECIFIC FAILURE CASES")
print("="*70)

# False negatives (Missed DILI - most dangerous)
fn_drugs = df[df['fn']].copy()
print(f"\nFALSE NEGATIVES (Missed DILI drugs): {len(fn_drugs)}")
fn_drugs = fn_drugs.sort_values('y_pred', ascending=False)
for _, row in fn_drugs.head(15).iterrows():
    logp = row.get('logp', 'N/A')
    ecnp = row.get('ecnp_z', 'N/A')
    k = row.get('k', row.get('n_targets', 'N/A'))
    print(f"  {row['drug_name']:25s} pred={row['y_pred']:.2f} logp={logp:.1f} ecnp_z={ecnp:.2f} k={k}")

# False positives (False alarms)
fp_drugs = df[df['fp']].copy()
print(f"\nFALSE POSITIVES (False alarms): {len(fp_drugs)}")
fp_drugs = fp_drugs.sort_values('y_pred', ascending=False)
for _, row in fp_drugs.head(15).iterrows():
    logp = row.get('logp', 'N/A')
    ecnp = row.get('ecnp_z', 'N/A')
    k = row.get('k', row.get('n_targets', 'N/A'))
    print(f"  {row['drug_name']:25s} pred={row['y_pred']:.2f} logp={logp:.1f} ecnp_z={ecnp:.2f} k={k}")

# =============================================================================
# 4. REGIME-SPECIFIC ERROR ANALYSIS
# =============================================================================

print("\n" + "="*70)
print("4. REGIME-SPECIFIC ERROR PATTERNS")
print("="*70)

# Define regimes
df['regime'] = 'unknown'
if 'ecnp_eligible' in df.columns:
    df.loc[df['ecnp_eligible'] == 1, 'regime'] = 'network'

if 'logp' in df.columns:
    df.loc[(df['logp'] >= 3) & (df['regime'] == 'unknown'), 'regime'] = 'high_logp'

for regime in df['regime'].unique():
    subset = df[df['regime'] == regime]
    if len(subset) > 0:
        print(f"\n{regime.upper()} (n={len(subset)}):")
        print(f"  Error rate: {subset['error'].mean()*100:.1f}%")
        print(f"  FN: {subset['fn'].sum()}, FP: {subset['fp'].sum()}")
        
        # DILI rate in this regime
        dili_rate = subset['is_dili'].mean() * 100
        print(f"  DILI rate: {dili_rate:.0f}%")

# =============================================================================
# 5. FEATURE BLIND SPOTS
# =============================================================================

print("\n" + "="*70)
print("5. FEATURE BLIND SPOTS")
print("="*70)

# Compare mean features for errors vs correct
feature_cols = ['logp', 'mw', 'ecnp_z', 'k', 'n_targets', 'hbd', 'hba', 'tpsa']
feature_cols = [f for f in feature_cols if f in df.columns]

print("\nMean feature values:")
print(f"{'Feature':<15} {'Correct':<12} {'Error':<12} {'Diff':<10}")
print("-" * 50)

for feat in feature_cols:
    correct_mean = df[df['correct']][feat].mean()
    error_mean = df[df['error']][feat].mean()
    diff = error_mean - correct_mean
    print(f"{feat:<15} {correct_mean:<12.2f} {error_mean:<12.2f} {diff:+.2f}")

# Check for patterns in errors
print("\n\nPOTENTIAL BLIND SPOTS:")

# Low target compounds
if 'k' in df.columns or 'n_targets' in df.columns:
    k_col = 'k' if 'k' in df.columns else 'n_targets'
    low_k = df[df[k_col] <= 3]
    if len(low_k) > 0:
        err_rate = low_k['error'].mean() * 100
        print(f"  Low target (k<=3): {err_rate:.0f}% error rate (n={len(low_k)})")

# Compounds with no network signal
if 'ecnp_eligible' in df.columns:
    non_elig = df[df['ecnp_eligible'] == 0]
    if len(non_elig) > 0:
        err_rate = non_elig['error'].mean() * 100
        print(f"  Non-ECNP-eligible: {err_rate:.0f}% error rate (n={len(non_elig)})")

# =============================================================================
# 6. SUMMARY
# =============================================================================

print("\n" + "="*70)
print("ACCURACY GAP SUMMARY")
print("="*70)

print(f"""
OVERALL PERFORMANCE:
  AUC: 0.884
  Accuracy: {accuracy:.1%}
  Errors: {df['error'].sum()}/{len(df)}

ERROR BREAKDOWN:
  False negatives (missed DILI): {df['fn'].sum()} - CRITICAL
  False positives (false alarms): {df['fp'].sum()}

KEY GAPS IDENTIFIED:
""")

# Identify key gaps
gaps = []
if 'logp' in df.columns:
    low_logp_err = df[df['logp'] < 2]['error'].mean() * 100
    if low_logp_err > 25:
        gaps.append(f"  - Low LogP compounds ({low_logp_err:.0f}% error)")

if 'ecnp_eligible' in df.columns:
    non_elig_err = df[df['ecnp_eligible'] == 0]['error'].mean() * 100
    if non_elig_err > 25:
        gaps.append(f"  - Non-ECNP-eligible compounds ({non_elig_err:.0f}% error)")

fn_count = df['fn'].sum()
if fn_count > 10:
    gaps.append(f"  - {fn_count} DILI drugs missed (false negatives)")

if gaps:
    for g in gaps:
        print(g)
else:
    print("  - No major gaps identified")

# Save analysis
output = ROOT / 'research' / 'ecnp-generalization' / 'results' / 'accuracy_gap_analysis.txt'
print(f"\nAnalysis complete.")

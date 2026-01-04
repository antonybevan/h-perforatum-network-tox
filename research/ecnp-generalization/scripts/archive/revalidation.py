"""
ML Pipeline Revalidation & Gap Analysis
========================================

Deeper validation to check for:
1. Overfitting / data leakage
2. Class imbalance effects
3. Holdout consistency
4. Feature correlation issues
5. Regime-specific failure modes
"""
import sys
sys.path.insert(0, r'v:\new\h-perforatum-network-tox\research\ecnp-closed-form\src\core')

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score, precision_recall_curve, average_precision_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict, train_test_split
from sklearn.preprocessing import StandardScaler
from scipy.stats import spearmanr
import warnings
warnings.filterwarnings('ignore')

try:
    from xgboost import XGBClassifier
    XGBOOST_AVAILABLE = True
except:
    XGBOOST_AVAILABLE = False
    from sklearn.linear_model import LogisticRegression

ROOT = Path(r'v:\new\h-perforatum-network-tox')

# Load previous results
df = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'results' / 'ml_pipeline_results.csv')
print(f"Loaded: {len(df)} compounds")
print(f"DILI+: {df['is_dili'].sum()}, DILI-: {(df['is_dili']==0).sum()}")

# =============================================================================
# 1. CHECK FOR DATA LEAKAGE
# =============================================================================

print("\n" + "="*60)
print("1. DATA LEAKAGE CHECK")
print("="*60)

# Feature correlation with target
feature_cols = ['logp', 'mw', 'hbd', 'hba', 'tpsa', 'k', 'n_dili_hits', 
                'n_network_hits', 'ecnp_z', 'dili_fraction', 'network_fraction']
                
for col in feature_cols:
    if col in df.columns:
        corr, p = spearmanr(df[col], df['is_dili'])
        if abs(corr) > 0.3:
            print(f"  WARNING: {col} has high correlation with target: r={corr:.3f}")
        elif abs(corr) > 0.15:
            print(f"  {col}: r={corr:.3f} (moderate)")

# Check if n_dili_hits is leaking
print(f"\nn_dili_hits correlation with is_dili: {spearmanr(df['n_dili_hits'], df['is_dili'])[0]:.3f}")
print("  (This is expected - but verify it's from network, not labels)")

# =============================================================================
# 2. TRAIN-TEST SPLIT VALIDATION
# =============================================================================

print("\n" + "="*60)
print("2. HOLDOUT VALIDATION (20% test)")
print("="*60)

X_cols = ['logp', 'mw', 'hbd', 'hba', 'tpsa', 'k', 'has_network_targets', 
          'has_transporter_targets', 'high_logp', 'ecnp_z', 'n_dili_hits', 
          'n_network_hits', 'dili_fraction']
          
X = df[X_cols].values
y = df['is_dili'].values

# Multiple random splits
test_aucs = []
for seed in [42, 123, 456, 789, 101112]:
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, 
                                                         random_state=seed, stratify=y)
    scaler = StandardScaler()
    X_train_s = scaler.fit_transform(X_train)
    X_test_s = scaler.transform(X_test)
    
    if XGBOOST_AVAILABLE:
        model = XGBClassifier(n_estimators=100, max_depth=3, random_state=42, eval_metric='logloss')
    else:
        model = LogisticRegression(max_iter=1000, random_state=42)
    
    model.fit(X_train_s, y_train)
    y_pred = model.predict_proba(X_test_s)[:, 1]
    auc = roc_auc_score(y_test, y_pred)
    test_aucs.append(auc)
    print(f"  Split {seed}: Test AUC = {auc:.3f}")

print(f"\nMean test AUC: {np.mean(test_aucs):.3f} (+/- {np.std(test_aucs):.3f})")

if np.std(test_aucs) > 0.05:
    print("  WARNING: High variance across splits - model may be unstable")

# =============================================================================
# 3. CLASS IMBALANCE ANALYSIS
# =============================================================================

print("\n" + "="*60)
print("3. CLASS IMBALANCE EFFECTS")
print("="*60)

# AUPRC (better for imbalanced data)
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

if XGBOOST_AVAILABLE:
    model = XGBClassifier(n_estimators=100, max_depth=3, random_state=42, eval_metric='logloss')
else:
    model = LogisticRegression(max_iter=1000, random_state=42)

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
y_pred = cross_val_predict(model, X_scaled, y, cv=cv, method='predict_proba')[:, 1]

auc = roc_auc_score(y, y_pred)
auprc = average_precision_score(y, y_pred)
baseline_pr = y.sum() / len(y)

print(f"ROC-AUC: {auc:.3f}")
print(f"AUPRC: {auprc:.3f} (baseline = {baseline_pr:.3f})")
print(f"AUPRC lift: {auprc / baseline_pr:.2f}x over random")

# =============================================================================
# 4. FAILURE MODE ANALYSIS
# =============================================================================

print("\n" + "="*60)
print("4. FAILURE MODE ANALYSIS")
print("="*60)

df['y_pred'] = y_pred
df['error'] = (df['y_pred'] > 0.5).astype(int) != df['is_dili']
df['fp'] = ((df['y_pred'] > 0.5) & (df['is_dili'] == 0)).astype(int)
df['fn'] = ((df['y_pred'] <= 0.5) & (df['is_dili'] == 1)).astype(int)

print(f"\nFalse positives: {df['fp'].sum()}")
print(f"False negatives: {df['fn'].sum()}")

# False negatives (missed DILI)
fn_drugs = df[df['fn'] == 1].nlargest(10, 'y_pred')
print(f"\nTop false negatives (missed DILI drugs):")
for _, row in fn_drugs.iterrows():
    print(f"  {row['drug_name']}: pred={row['y_pred']:.2f}, Z={row['ecnp_z']:.2f}, n_targets={row['k']}")

# False positives (false alarms)
fp_drugs = df[df['fp'] == 1].nlargest(10, 'y_pred')
print(f"\nTop false positives:")
for _, row in fp_drugs.iterrows():
    print(f"  {row['drug_name']}: pred={row['y_pred']:.2f}, Z={row['ecnp_z']:.2f}, n_targets={row['k']}")

# =============================================================================
# 5. REGIME-SPECIFIC GAPS
# =============================================================================

print("\n" + "="*60)
print("5. REGIME-SPECIFIC GAP ANALYSIS")
print("="*60)

# By ECNP eligibility
for eligible in [0, 1]:
    subset = df[df['ecnp_eligible'] == eligible]
    errors = subset['error'].sum()
    total = len(subset)
    error_rate = errors / total * 100
    
    label = "ECNP-eligible" if eligible else "Non-eligible"
    print(f"\n{label} (n={total}):")
    print(f"  Errors: {errors} ({error_rate:.1f}%)")
    print(f"  FN: {subset['fn'].sum()}, FP: {subset['fp'].sum()}")

# By target count
for k_range in [(3, 4), (5, 7), (8, 15), (16, 100)]:
    subset = df[(df['k'] >= k_range[0]) & (df['k'] <= k_range[1])]
    if len(subset) >= 10:
        y_t = subset['is_dili'].values
        y_p = subset['y_pred'].values
        n_pos = y_t.sum()
        n_neg = len(y_t) - n_pos
        if n_pos >= 3 and n_neg >= 3:
            auc = roc_auc_score(y_t, y_p)
            print(f"\nk={k_range[0]}-{k_range[1]} (n={len(subset)}): AUC = {auc:.3f}")

# =============================================================================
# 6. SUMMARY
# =============================================================================

print("\n" + "="*60)
print("REVALIDATION SUMMARY")
print("="*60)

print(f"""
Cross-validation AUC:  {auc:.3f}
Mean holdout AUC:      {np.mean(test_aucs):.3f} (+/- {np.std(test_aucs):.3f})
AUPRC:                 {auprc:.3f} ({auprc/baseline_pr:.1f}x baseline)

IDENTIFIED GAPS:
""")

gaps = []

if np.std(test_aucs) > 0.05:
    gaps.append("- High variance across splits (unstable model)")
    
if auprc < 0.8:
    gaps.append("- AUPRC below 0.8 suggests precision issues")
    
fn_rate = df['fn'].sum() / df['is_dili'].sum() * 100
if fn_rate > 20:
    gaps.append(f"- False negative rate {fn_rate:.0f}% - missing DILI drugs")

ecnp_non_eligible_errors = df[df['ecnp_eligible'] == 0]['error'].sum() / len(df[df['ecnp_eligible'] == 0]) * 100
if ecnp_non_eligible_errors > 30:
    gaps.append(f"- High error rate ({ecnp_non_eligible_errors:.0f}%) in non-ECNP-eligible compounds")

if gaps:
    for gap in gaps:
        print(gap)
else:
    print("- No major gaps identified")

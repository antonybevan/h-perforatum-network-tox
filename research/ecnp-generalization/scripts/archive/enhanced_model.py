"""
Enhanced DILI Model with Rule-of-Two and CYP Metabolism
=========================================================

Adds mechanistically significant features:
1. FDA Rule-of-Two (RO2): dose ≥100mg AND LogP ≥3
2. CYP metabolism proxy: CYP target + high LogP
3. Improved regime classification

Goal: Improve prediction on dose-dependent intrinsic DILI
"""
import sys
sys.path.insert(0, r'v:\new\h-perforatum-network-tox\research\ecnp-closed-form\src\core')

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

try:
    from xgboost import XGBClassifier
    XGBOOST_AVAILABLE = True
except:
    XGBOOST_AVAILABLE = False
    from sklearn.linear_model import LogisticRegression

ROOT = Path(r'v:\new\h-perforatum-network-tox')

# =============================================================================
# LOAD DATA
# =============================================================================

print("Loading data...")
df = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'results' / 'ml_pipeline_with_alerts.csv')

# Add dose data if available
pk = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'data' / 'curated' / 'pk_data.csv')
pk['drug_name_lower'] = pk['drug_name'].str.lower().str.strip()
pk_unique = pk.groupby('drug_name_lower').first().reset_index()

df['drug_name_lower'] = df['drug_name'].str.lower().str.strip()
df = df.merge(pk_unique[['drug_name_lower', 'dose_mg', 'high_dose']], 
               on='drug_name_lower', how='left', suffixes=('', '_pk'))

# Use high_dose from PK if available, otherwise from existing
if 'high_dose' not in df.columns or df['high_dose'].isna().all():
    if 'high_dose_pk' in df.columns:
        df['high_dose'] = df['high_dose_pk']
df['high_dose'] = df['high_dose'].fillna(0)
df['dose_mg'] = df['dose_mg'].fillna(df['dose_mg'].median())

print(f"Compounds: {len(df)}")

# =============================================================================
# CREATE MECHANISTICALLY SIGNIFICANT FEATURES
# =============================================================================

print("\nCreating mechanistic features...")

# 1. FDA Rule-of-Two (RO2)
# High DILI risk = dose >= 100mg AND LogP >= 3
df['high_logp'] = (df['logp'] >= 3).astype(int)
df['rule_of_two'] = ((df['high_dose'] >= 1) & (df['high_logp'] >= 1)).astype(int)

print(f"Rule-of-Two satisfied: {df['rule_of_two'].sum()}/{len(df)}")

# 2. CYP metabolism proxy
# CYP targets + high LogP = bioactivation risk
CYP_GENES = {'CYP1A1', 'CYP1A2', 'CYP2A6', 'CYP2B6', 'CYP2C8', 'CYP2C9', 
             'CYP2C19', 'CYP2D6', 'CYP2E1', 'CYP3A4', 'CYP3A5'}

import ast
def parse_targets(t):
    try:
        return ast.literal_eval(t) if isinstance(t, str) else []
    except:
        return []

if 'targets' in df.columns:
    df['target_list'] = df['targets'].apply(parse_targets)
    df['n_cyp_targets'] = df['target_list'].apply(lambda t: len(set(t) & CYP_GENES))
else:
    df['n_cyp_targets'] = 0

df['has_cyp_targets'] = (df['n_cyp_targets'] > 0).astype(int)
df['bioactivation_risk'] = ((df['has_cyp_targets'] == 1) & (df['high_logp'] == 1)).astype(int)

print(f"CYP targets present: {df['has_cyp_targets'].sum()}/{len(df)}")
print(f"Bioactivation risk: {df['bioactivation_risk'].sum()}/{len(df)}")

# 3. Improved regime classification
# Intrinsic: high dose + high LogP + CYP
# Network: DILI gene targets + ECNP eligible
# Unknown: neither

df['intrinsic_regime'] = (
    (df['rule_of_two'] == 1) | 
    (df['bioactivation_risk'] == 1)
).astype(int)

df['network_regime'] = df.get('ecnp_eligible', 0)

df['regime'] = 'unknown'
df.loc[df['intrinsic_regime'] == 1, 'regime'] = 'intrinsic'
df.loc[(df['network_regime'] == 1) & (df['intrinsic_regime'] == 0), 'regime'] = 'network'

print(f"\nRegime distribution:")
print(df['regime'].value_counts())

# =============================================================================
# ANALYZE NEW FEATURES
# =============================================================================

print("\n" + "="*60)
print("NEW FEATURE ANALYSIS")
print("="*60)

# DILI rate by new features
for feature in ['rule_of_two', 'bioactivation_risk', 'intrinsic_regime']:
    for val in [0, 1]:
        subset = df[df[feature] == val]
        if len(subset) > 0:
            dili_rate = subset['is_dili'].mean() * 100
            label = "Yes" if val == 1 else "No"
            print(f"{feature} = {label}: {subset['is_dili'].sum()}/{len(subset)} ({dili_rate:.0f}%)")
    print()

# =============================================================================
# BUILD ENHANCED MODEL
# =============================================================================

print("="*60)
print("ENHANCED MODEL WITH RO2 + CYP")
print("="*60)

# All features including new ones
all_features = [
    # Chemistry
    'logp', 'mw', 'hbd', 'hba', 'tpsa',
    # Targets
    'k', 'has_network_targets', 'n_dili_hits', 'dili_fraction',
    # ECNP
    'ecnp_z',
    # Alerts
    'n_alerts',
    # NEW: RO2 and CYP
    'rule_of_two', 'bioactivation_risk', 'intrinsic_regime',
    'has_cyp_targets', 'high_logp'
]

# Filter valid columns
feature_cols = [f for f in all_features if f in df.columns]

valid_mask = df[feature_cols].notna().all(axis=1)
df_valid = df[valid_mask].copy()
print(f"\nValid samples: {len(df_valid)}")

X = df_valid[feature_cols].values
y = df_valid['is_dili'].values

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

if XGBOOST_AVAILABLE:
    model = XGBClassifier(n_estimators=100, max_depth=3, random_state=42, eval_metric='logloss')
else:
    model = LogisticRegression(max_iter=1000, random_state=42)

y_pred = cross_val_predict(model, X_scaled, y, cv=cv, method='predict_proba')[:, 1]
auc_enhanced = roc_auc_score(y, y_pred)

print(f"\nEnhanced model AUC: {auc_enhanced:.3f}")
print(f"Previous model AUC: 0.802")
print(f"Improvement: {auc_enhanced - 0.802:+.3f}")

# =============================================================================
# FEATURE IMPORTANCE
# =============================================================================

print("\n" + "="*60)
print("FEATURE IMPORTANCE")
print("="*60)

model.fit(X_scaled, y)
if XGBOOST_AVAILABLE:
    importances = model.feature_importances_
else:
    importances = np.abs(model.coef_[0])

for name, imp in sorted(zip(feature_cols, importances), key=lambda x: x[1], reverse=True):
    print(f"  {name:20s}: {imp:.3f}")

# =============================================================================
# REGIME-SPECIFIC PERFORMANCE
# =============================================================================

print("\n" + "="*60)
print("REGIME-SPECIFIC AUC")
print("="*60)

df_valid['y_pred_enhanced'] = y_pred

for regime in ['intrinsic', 'network', 'unknown']:
    subset = df_valid[df_valid['regime'] == regime]
    if len(subset) >= 15:
        y_t = subset['is_dili'].values
        y_p = subset['y_pred_enhanced'].values
        n_pos = y_t.sum()
        n_neg = len(y_t) - n_pos
        if n_pos >= 3 and n_neg >= 3:
            auc = roc_auc_score(y_t, y_p)
            error = ((y_p > 0.5).astype(int) != y_t).sum()
            error_rate = error / len(y_t) * 100
            print(f"\n{regime.upper()} regime (n={len(subset)}, {n_pos}:{n_neg}):")
            print(f"  AUC = {auc:.3f}")
            print(f"  Error rate = {error_rate:.1f}%")
        else:
            print(f"\n{regime.upper()} regime: insufficient class balance")
    else:
        print(f"\n{regime.upper()} regime: n={len(subset)} (too few)")

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "="*60)
print("SUMMARY")
print("="*60)

print(f"""
New features added:
  - Rule-of-Two (dose ≥100mg + LogP ≥3): {df['rule_of_two'].sum()} compounds
  - Bioactivation risk (CYP + high LogP): {df['bioactivation_risk'].sum()} compounds
  - Intrinsic regime: {(df['regime']=='intrinsic').sum()} compounds

Model performance:
  - Previous AUC: 0.802
  - Enhanced AUC: {auc_enhanced:.3f}
  - Change: {auc_enhanced - 0.802:+.3f}
""")

# Save results
output = ROOT / 'research' / 'ecnp-generalization' / 'results' / 'enhanced_model_results.csv'
df_valid.to_csv(output, index=False)
print(f"Saved: {output}")

"""
Dose-Based Stratification for Non-Eligible Regime
==================================================

Uses PK data (half-life, dose) to stratify non-eligible compounds
and check if dose-aware modeling reduces error rate.
"""
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
pk = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'data' / 'curated' / 'pk_data.csv')

print(f"Compounds: {len(df)}")
print(f"PK data rows: {len(pk)}")

# Match PK by drug name
df['drug_name_lower'] = df['drug_name'].str.lower().str.strip()
pk['drug_name_lower'] = pk['drug_name'].str.lower().str.strip()

# Get unique PK per drug (take first if duplicates)
pk_unique = pk.groupby('drug_name_lower').first().reset_index()
print(f"Unique PK drugs: {len(pk_unique)}")

# Merge
merged = df.merge(pk_unique[['drug_name_lower', 'half_life_hours', 'half_life_bin', 
                              'dose_mg', 'high_dose', 'route']], 
                   on='drug_name_lower', how='left')
print(f"Merged rows: {len(merged)}")

# Check PK coverage
print(f"\nPK coverage:")
print(f"  Half-life: {merged['half_life_bin'].notna().sum()}/{len(merged)}")
print(f"  Dose: {merged['dose_mg'].notna().sum()}/{len(merged)}")
print(f"  High dose: {merged['high_dose'].notna().sum()}/{len(merged)}")

# =============================================================================
# ANALYZE NON-ELIGIBLE BY DOSE
# =============================================================================

print("\n" + "="*60)
print("NON-ELIGIBLE REGIME: DOSE STRATIFICATION")
print("="*60)

non_elig = merged[merged['ecnp_eligible'] == 0].copy()
print(f"Non-eligible compounds: {len(non_elig)}")

# Fill missing dose with median
non_elig['dose_mg'] = non_elig['dose_mg'].fillna(non_elig['dose_mg'].median())
non_elig['high_dose'] = non_elig['high_dose'].fillna(0)
non_elig['half_life_bin'] = non_elig['half_life_bin'].fillna('unknown')

# Stratify by dose
print(f"\nHigh dose distribution:")
print(non_elig['high_dose'].value_counts())

print(f"\nHalf-life distribution:")
print(non_elig['half_life_bin'].value_counts())

# DILI by dose
high_dose = non_elig[non_elig['high_dose'] == 1]
low_dose = non_elig[non_elig['high_dose'] == 0]

print(f"\nDILI by dose:")
if len(high_dose) > 0:
    print(f"  High dose: {high_dose['is_dili'].sum()}/{len(high_dose)} ({high_dose['is_dili'].mean()*100:.0f}%)")
if len(low_dose) > 0:
    print(f"  Low dose: {low_dose['is_dili'].sum()}/{len(low_dose)} ({low_dose['is_dili'].mean()*100:.0f}%)")

# =============================================================================
# BUILD DOSE-AWARE MODEL
# =============================================================================

print("\n" + "="*60)
print("DOSE-AWARE MODEL FOR NON-ELIGIBLE")
print("="*60)

# Features with dose
dose_features = [
    'logp', 'mw', 'hbd', 'hba', 'tpsa', 'k',
    'n_alerts', 'dose_mg', 'high_dose'
]

# Check feature availability
print(f"Feature availability:")
for f in dose_features:
    if f in non_elig.columns:
        available = non_elig[f].notna().sum()
        print(f"  {f}: {available}/{len(non_elig)}")

# Filter to valid rows
valid_mask = non_elig[dose_features].notna().all(axis=1)
non_elig_valid = non_elig[valid_mask].copy()
print(f"\nValid for dose-aware model: {len(non_elig_valid)}")

if len(non_elig_valid) >= 30:
    X = non_elig_valid[dose_features].values
    y = non_elig_valid['is_dili'].values
    
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    
    if XGBOOST_AVAILABLE:
        model = XGBClassifier(n_estimators=100, max_depth=3, random_state=42, eval_metric='logloss')
    else:
        model = LogisticRegression(max_iter=1000, random_state=42)
    
    y_pred = cross_val_predict(model, X_scaled, y, cv=cv, method='predict_proba')[:, 1]
    auc = roc_auc_score(y, y_pred)
    
    print(f"\nDose-aware model AUC: {auc:.3f}")
    print(f"Previous non-eligible AUC: 0.739")
    print(f"Improvement: {auc - 0.739:+.3f}")
    
    # Error rate
    errors = ((y_pred > 0.5).astype(int) != y).sum()
    error_rate = errors / len(y) * 100
    print(f"\nError rate: {error_rate:.1f}% (was 30.7%)")
    
    # Feature importance
    model.fit(X_scaled, y)
    if XGBOOST_AVAILABLE:
        importances = model.feature_importances_
    else:
        importances = np.abs(model.coef_[0])
    
    print(f"\nFeature importance:")
    for name, imp in sorted(zip(dose_features, importances), key=lambda x: x[1], reverse=True):
        print(f"  {name:12s}: {imp:.3f}")

else:
    print(f"Not enough valid samples for dose-aware model: {len(non_elig_valid)}")

# =============================================================================
# FULL PIPELINE WITH DOSE
# =============================================================================

print("\n" + "="*60)
print("FULL PIPELINE WITH DOSE FEATURES")
print("="*60)

# Add dose to full dataset
merged['dose_mg'] = merged['dose_mg'].fillna(merged['dose_mg'].median())
merged['high_dose'] = merged['high_dose'].fillna(0)

full_features = [
    'logp', 'mw', 'hbd', 'hba', 'tpsa', 'k',
    'has_network_targets', 'ecnp_z', 'n_dili_hits', 'dili_fraction',
    'n_alerts', 'dose_mg', 'high_dose'
]

valid_full = merged[full_features].notna().all(axis=1)
merged_valid = merged[valid_full].copy()
print(f"Valid for full model: {len(merged_valid)}")

X_full = merged_valid[full_features].values
y_full = merged_valid['is_dili'].values

scaler_full = StandardScaler()
X_full_scaled = scaler_full.fit_transform(X_full)

if XGBOOST_AVAILABLE:
    model_full = XGBClassifier(n_estimators=100, max_depth=3, random_state=42, eval_metric='logloss')
else:
    model_full = LogisticRegression(max_iter=1000, random_state=42)

y_pred_full = cross_val_predict(model_full, X_full_scaled, y_full, cv=cv, method='predict_proba')[:, 1]
auc_full = roc_auc_score(y_full, y_pred_full)

print(f"\nFull pipeline with dose AUC: {auc_full:.3f}")
print(f"Previous (no dose): 0.802")
print(f"Improvement: {auc_full - 0.802:+.3f}")

# Regime-specific
merged_valid['y_pred_dose'] = y_pred_full

for regime, mask in [('ECNP-eligible', merged_valid['ecnp_eligible'] == 1),
                      ('Non-eligible', merged_valid['ecnp_eligible'] == 0)]:
    subset = merged_valid[mask]
    if len(subset) >= 10:
        y_t = subset['is_dili'].values
        y_p = subset['y_pred_dose'].values
        if y_t.sum() >= 3 and (len(y_t) - y_t.sum()) >= 3:
            auc_r = roc_auc_score(y_t, y_p)
            errors = ((y_p > 0.5).astype(int) != y_t).sum()
            error_rate = errors / len(y_t) * 100
            print(f"\n{regime}:")
            print(f"  AUC = {auc_r:.3f}")
            print(f"  Error rate = {error_rate:.1f}%")

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "="*60)
print("DOSE STRATIFICATION SUMMARY")
print("="*60)

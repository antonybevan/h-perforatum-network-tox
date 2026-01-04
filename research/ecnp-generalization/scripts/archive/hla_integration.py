"""
HLA-DILI Integration
====================

Adds HLA genetic risk flag to compounds and checks impact on false negatives.
"""
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.ensemble import RandomForestClassifier
import warnings
warnings.filterwarnings('ignore')

ROOT = Path(r'v:\new\h-perforatum-network-tox')

# Load data
print("Loading data...")
df = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'results' / 'ecfp_model_results.csv')
hla = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'data' / 'curated' / 'hla_dili_associations.csv')

print(f"Compounds: {len(df)}")
print(f"HLA associations: {len(hla)}")

# Create lookup of drugs with HLA risk
hla_drugs = set(hla['drug_name'].str.lower().str.strip())
print(f"Unique HLA-risk drugs: {len(hla_drugs)}")

# Match drugs in our dataset
df['drug_name_lower'] = df['drug_name'].str.lower().str.strip()
df['hla_risk'] = df['drug_name_lower'].isin(hla_drugs).astype(int)

print(f"\nDrugs matched with HLA risk: {df['hla_risk'].sum()}/{len(df)}")

# Show HLA-risk drugs in our dataset
hla_matched = df[df['hla_risk'] == 1][['drug_name', 'is_dili', 'y_pred_combined']].copy()
if len(hla_matched) > 0:
    print(f"\nHLA-risk drugs found:")
    for _, row in hla_matched.iterrows():
        status = "DILI+" if row['is_dili'] == 1 else "DILI-"
        pred = row.get('y_pred_combined', row.get('y_pred_ecfp', 'N/A'))
        print(f"  {row['drug_name']}: {status}, pred={pred:.2f}")

# =============================================================================
# CHECK FALSE NEGATIVES WITH HLA RISK
# =============================================================================

print("\n" + "="*60)
print("FALSE NEGATIVE ANALYSIS WITH HLA RISK")
print("="*60)

# Identify false negatives from previous analysis
if 'y_pred_combined' in df.columns:
    df['y_pred'] = df['y_pred_combined']
elif 'y_pred_ecfp' in df.columns:
    df['y_pred'] = df['y_pred_ecfp']

df['y_pred_binary'] = (df['y_pred'] > 0.5).astype(int)
df['fn'] = (df['y_pred_binary'] == 0) & (df['is_dili'] == 1)

fn_drugs = df[df['fn']]
print(f"\nFalse negatives: {len(fn_drugs)}")

# Check HLA overlap with false negatives
fn_with_hla = fn_drugs[fn_drugs['hla_risk'] == 1]
print(f"FN with HLA risk: {len(fn_with_hla)}")

if len(fn_with_hla) > 0:
    print("\nFalse negatives that have HLA risk (can now flag):")
    for _, row in fn_with_hla.iterrows():
        print(f"  {row['drug_name']}: pred={row['y_pred']:.2f}")

# =============================================================================
# BUILD ENHANCED MODEL WITH HLA FLAG
# =============================================================================

print("\n" + "="*60)
print("MODEL WITH HLA RISK FLAG")
print("="*60)

# Add HLA flag to features
feature_cols = ['logp', 'mw', 'hbd', 'hba', 'tpsa', 'k', 'ecnp_z', 
                'has_network_targets', 'dili_fraction', 'hla_risk']
feature_cols = [f for f in feature_cols if f in df.columns]

valid_mask = df[feature_cols].notna().all(axis=1)
df_valid = df[valid_mask].copy()

X = df_valid[feature_cols].values
y = df_valid['is_dili'].values

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
model = RandomForestClassifier(n_estimators=100, max_depth=10, random_state=42, n_jobs=-1)

y_pred_hla = cross_val_predict(model, X, y, cv=cv, method='predict_proba')[:, 1]
auc_with_hla = roc_auc_score(y, y_pred_hla)

print(f"AUC with HLA flag: {auc_with_hla:.3f}")
print(f"Previous AUC: 0.884")
print(f"Change: {auc_with_hla - 0.884:+.3f}")

# Check new false negatives
df_valid['y_pred_hla'] = y_pred_hla
df_valid['fn_new'] = (df_valid['y_pred_hla'] <= 0.5) & (df_valid['is_dili'] == 1)
new_fn_count = df_valid['fn_new'].sum()
print(f"\nFalse negatives with HLA model: {new_fn_count} (was 9)")

# =============================================================================
# HYBRID APPROACH: HLA FLAG + MODEL
# =============================================================================

print("\n" + "="*60)
print("HYBRID APPROACH: MODEL + HLA OVERRIDE")
print("="*60)

# If a drug has HLA risk and model pred < 0.5, override to high risk
df_valid['y_pred_hybrid'] = df_valid['y_pred_hla'].copy()
hla_override = (df_valid['hla_risk'] == 1) & (df_valid['y_pred_hla'] < 0.7)
df_valid.loc[hla_override, 'y_pred_hybrid'] = 0.8  # Boost HLA-risk drugs

auc_hybrid = roc_auc_score(df_valid['is_dili'], df_valid['y_pred_hybrid'])
print(f"AUC with HLA override: {auc_hybrid:.3f}")

# New performance
df_valid['fn_hybrid'] = (df_valid['y_pred_hybrid'] <= 0.5) & (df_valid['is_dili'] == 1)
df_valid['fp_hybrid'] = (df_valid['y_pred_hybrid'] > 0.5) & (df_valid['is_dili'] == 0)
print(f"False negatives: {df_valid['fn_hybrid'].sum()}")
print(f"False positives: {df_valid['fp_hybrid'].sum()}")

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "="*60)
print("SUMMARY: HLA INTEGRATION IMPACT")
print("="*60)

print(f"""
HLA-DILI associations loaded: {len(hla)}
Drugs matched in dataset: {df['hla_risk'].sum()}

Model Performance:
  Previous AUC:       0.884
  With HLA flag:      {auc_with_hla:.3f}
  With HLA override:  {auc_hybrid:.3f}

False Negatives:
  Before: 9
  With HLA model: {new_fn_count}
  With HLA override: {df_valid['fn_hybrid'].sum()}
""")

# Save
output = ROOT / 'research' / 'ecnp-generalization' / 'results' / 'hla_enhanced_results.csv'
df_valid.to_csv(output, index=False)
print(f"Saved: {output}")

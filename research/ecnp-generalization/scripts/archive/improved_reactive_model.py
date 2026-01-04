"""
Improved Reactive/Unknown Regime Model
=======================================

The non-eligible regime has 29% error rate.
This module adds:
1. Structural alerts for reactive metabolites (SMARTS patterns)
2. Improved chemistry features
3. Calibrated risk scoring for non-network compounds
"""
import sys
sys.path.insert(0, r'v:\new\h-perforatum-network-tox\research\ecnp-closed-form\src\core')

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score, precision_recall_curve, average_precision_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("RDKit required for structural alerts")

try:
    from xgboost import XGBClassifier
    XGBOOST_AVAILABLE = True
except:
    XGBOOST_AVAILABLE = False
    from sklearn.linear_model import LogisticRegression

ROOT = Path(r'v:\new\h-perforatum-network-tox')

# =============================================================================
# STRUCTURAL ALERTS FOR REACTIVE METABOLITES
# =============================================================================

# Known hepatotoxicity structural alerts (SMARTS patterns)
# From literature: Stepan et al., Chen et al., FDA guidance
STRUCTURAL_ALERTS = {
    # Quinone precursors
    'hydroquinone': 'c1cc(O)ccc1O',
    'catechol': 'c1ccc(O)c(O)c1',
    'aminophenol': 'c1cc(N)ccc1O',
    
    # Reactive metabolite precursors
    'aniline': 'c1ccc(N)cc1',
    'nitroaromatic': 'c1ccccc1[N+](=O)[O-]',
    'thiophene': 'c1ccsc1',
    'furan': 'c1ccoc1',
    
    # Acyl glucuronide risk
    'carboxylic_acid': 'C(=O)O',
    
    # Covalent binding risk
    'epoxide': 'C1OC1',
    'michael_acceptor': 'C=CC=O',
    'halogenated_aromatic': 'c1ccc(F)cc1',
    'chlorinated_aromatic': 'c1ccc(Cl)cc1',
    
    # Mitochondrial toxicity
    'lipophilic_amine': 'CCCCN',
    'benzimidazole': 'c1ccc2[nH]cnc2c1',
    
    # Idiosyncratic DILI signals
    'hydrazine': 'NN',
    'isoniazid_like': 'c1ccncc1C(=O)NN',
}

def count_structural_alerts(smiles):
    """Count number of structural alerts in molecule."""
    if not RDKIT_AVAILABLE:
        return 0, {}
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 0, {}
    except:
        return 0, {}
    
    alert_counts = {}
    total = 0
    
    for name, smarts in STRUCTURAL_ALERTS.items():
        try:
            pattern = Chem.MolFromSmarts(smarts)
            if pattern:
                matches = mol.GetSubstructMatches(pattern)
                if matches:
                    alert_counts[name] = len(matches)
                    total += len(matches)
        except:
            pass
    
    return total, alert_counts

# =============================================================================
# LOAD DATA
# =============================================================================

print("Loading data...")
df = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'results' / 'ml_pipeline_results.csv')

# Add structural alerts
print("Computing structural alerts...")
df['n_alerts'] = 0
df['alert_names'] = ''

for idx, row in df.iterrows():
    if pd.notna(row.get('smiles')):
        n_alerts, alerts = count_structural_alerts(row['smiles'])
        df.at[idx, 'n_alerts'] = n_alerts
        df.at[idx, 'alert_names'] = ','.join(alerts.keys()) if alerts else ''

print(f"Compounds with alerts: {(df['n_alerts'] > 0).sum()}/{len(df)}")

# =============================================================================
# ANALYZE NON-ELIGIBLE COMPOUNDS
# =============================================================================

print("\n" + "="*60)
print("NON-ELIGIBLE REGIME ANALYSIS")
print("="*60)

non_eligible = df[df['ecnp_eligible'] == 0].copy()
print(f"Non-eligible compounds: {len(non_eligible)}")
print(f"DILI+: {non_eligible['is_dili'].sum()}, DILI-: {(non_eligible['is_dili']==0).sum()}")

# Alert distribution
print(f"\nStructural alerts in non-eligible:")
print(f"  Mean alerts: {non_eligible['n_alerts'].mean():.2f}")
print(f"  With any alert: {(non_eligible['n_alerts'] > 0).sum()}")

# Alert by DILI status
print(f"\n  DILI+ mean alerts: {non_eligible[non_eligible['is_dili']==1]['n_alerts'].mean():.2f}")
print(f"  DILI- mean alerts: {non_eligible[non_eligible['is_dili']==0]['n_alerts'].mean():.2f}")

# =============================================================================
# BUILD IMPROVED REACTIVE MODEL
# =============================================================================

print("\n" + "="*60)
print("BUILDING IMPROVED REACTIVE MODEL")
print("="*60)

# Features for non-eligible regime
reactive_features = [
    'logp', 'mw', 'hbd', 'hba', 'tpsa',  # Chemistry
    'n_alerts',  # Structural alerts
    'k',  # Target count (limited info)
]

# Filter to non-eligible with valid features
valid_mask = non_eligible[reactive_features].notna().all(axis=1)
non_eligible_valid = non_eligible[valid_mask].copy()
print(f"Valid non-eligible: {len(non_eligible_valid)}")

X_ne = non_eligible_valid[reactive_features].values
y_ne = non_eligible_valid['is_dili'].values

# Scale
scaler = StandardScaler()
X_ne_scaled = scaler.fit_transform(X_ne)

# Cross-validate
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

if XGBOOST_AVAILABLE:
    model = XGBClassifier(n_estimators=100, max_depth=3, random_state=42, eval_metric='logloss')
else:
    model = LogisticRegression(max_iter=1000, random_state=42)

y_pred_ne = cross_val_predict(model, X_ne_scaled, y_ne, cv=cv, method='predict_proba')[:, 1]
auc_ne = roc_auc_score(y_ne, y_pred_ne)
print(f"\nNon-eligible regime AUC: {auc_ne:.3f}")

# Compare to previous (chemistry only)
if 'y_pred' in non_eligible_valid.columns:
    auc_prev = roc_auc_score(y_ne, non_eligible_valid['y_pred'].values)
    print(f"Previous model on non-eligible: {auc_prev:.3f}")
    print(f"Improvement: {auc_ne - auc_prev:+.3f}")

# =============================================================================
# FEATURE IMPORTANCE
# =============================================================================

print("\n" + "="*60)
print("FEATURE IMPORTANCE (Reactive Model)")
print("="*60)

model.fit(X_ne_scaled, y_ne)
if XGBOOST_AVAILABLE:
    importances = model.feature_importances_
else:
    importances = np.abs(model.coef_[0])

for name, imp in sorted(zip(reactive_features, importances), key=lambda x: x[1], reverse=True):
    print(f"  {name:15s}: {imp:.3f}")

# =============================================================================
# COMBINED MODEL (FULL PIPELINE UPDATE)
# =============================================================================

print("\n" + "="*60)
print("UPDATED FULL PIPELINE")
print("="*60)

# Add alert features to full dataset
all_features = [
    'logp', 'mw', 'hbd', 'hba', 'tpsa', 'k', 
    'has_network_targets', 'has_transporter_targets', 'high_logp',
    'ecnp_z', 'n_dili_hits', 'n_network_hits', 'dili_fraction',
    'n_alerts'  # NEW: structural alerts
]

X_all = df[all_features].values
y_all = df['is_dili'].values

scaler_all = StandardScaler()
X_all_scaled = scaler_all.fit_transform(X_all)

if XGBOOST_AVAILABLE:
    model_full = XGBClassifier(n_estimators=100, max_depth=3, random_state=42, eval_metric='logloss')
else:
    model_full = LogisticRegression(max_iter=1000, random_state=42)

y_pred_full = cross_val_predict(model_full, X_all_scaled, y_all, cv=cv, method='predict_proba')[:, 1]
auc_full = roc_auc_score(y_all, y_pred_full)

print(f"Previous full model AUC: 0.792")
print(f"Updated with alerts AUC: {auc_full:.3f}")
print(f"Improvement: {auc_full - 0.792:+.3f}")

# Regime-specific with alerts
df['y_pred_updated'] = y_pred_full

ecnp_elig = df[df['ecnp_eligible'] == 1]
non_elig = df[df['ecnp_eligible'] == 0]

if len(ecnp_elig) > 10:
    y_t = ecnp_elig['is_dili'].values
    y_p = ecnp_elig['y_pred_updated'].values
    if y_t.sum() >= 3 and (len(y_t) - y_t.sum()) >= 3:
        print(f"\nECNP-eligible AUC: {roc_auc_score(y_t, y_p):.3f}")

if len(non_elig) > 10:
    y_t = non_elig['is_dili'].values
    y_p = non_elig['y_pred_updated'].values
    if y_t.sum() >= 3 and (len(y_t) - y_t.sum()) >= 3:
        print(f"Non-eligible AUC: {roc_auc_score(y_t, y_p):.3f}")

# Error rate check
non_elig['error_updated'] = (non_elig['y_pred_updated'] > 0.5).astype(int) != non_elig['is_dili']
error_rate_new = non_elig['error_updated'].sum() / len(non_elig) * 100
print(f"\nNon-eligible error rate: {error_rate_new:.1f}% (was 29.2%)")

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "="*60)
print("SUMMARY: REACTIVE MODEL IMPROVEMENT")
print("="*60)

print(f"""
Structural alerts added: {len(STRUCTURAL_ALERTS)} patterns
Compounds with alerts: {(df['n_alerts'] > 0).sum()}/{len(df)}

Non-eligible regime:
  Previous AUC: ~0.78
  Updated AUC:  {roc_auc_score(non_elig['is_dili'], non_elig['y_pred_updated']):.3f}
  Error rate:   {error_rate_new:.1f}% (was 29.2%)

Full pipeline:
  Previous AUC: 0.792
  Updated AUC:  {auc_full:.3f}
""")

# Save updated results
output = ROOT / 'research' / 'ecnp-generalization' / 'results' / 'ml_pipeline_with_alerts.csv'
df.to_csv(output, index=False)
print(f"Saved: {output}")

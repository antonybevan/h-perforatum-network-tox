"""
ECFP4 Fingerprint Model Test
============================

Tests whether ECFP4 molecular fingerprints can achieve the high AUC (0.95+)
claimed in literature. Compares to our current feature-based approach.
"""
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
import warnings
warnings.filterwarnings('ignore')

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("ERROR: RDKit required for ECFP fingerprints")
    exit()

try:
    from xgboost import XGBClassifier
    XGBOOST_AVAILABLE = True
except:
    XGBOOST_AVAILABLE = False

ROOT = Path(r'v:\new\h-perforatum-network-tox')

# =============================================================================
# LOAD DATA
# =============================================================================

print("Loading data...")
df = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'results' / 'ml_pipeline_with_alerts.csv')
print(f"Compounds: {len(df)}")

# =============================================================================
# GENERATE ECFP4 FINGERPRINTS
# =============================================================================

print("\nGenerating ECFP4 fingerprints...")

def smiles_to_ecfp(smiles, radius=2, nBits=1024):
    """Convert SMILES to ECFP4 fingerprint (radius=2, 1024 bits)."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits)
            return np.array(fp)
    except:
        pass
    return None

# Generate fingerprints
fingerprints = []
valid_idx = []

for idx, row in df.iterrows():
    if pd.notna(row.get('smiles')):
        fp = smiles_to_ecfp(row['smiles'])
        if fp is not None:
            fingerprints.append(fp)
            valid_idx.append(idx)

print(f"Valid fingerprints: {len(fingerprints)}/{len(df)}")

# Create fingerprint matrix
X_fp = np.array(fingerprints)
df_valid = df.loc[valid_idx].copy()
y = df_valid['is_dili'].values

print(f"Fingerprint shape: {X_fp.shape}")
print(f"Bits set (mean): {X_fp.sum(axis=1).mean():.1f}")

# =============================================================================
# MODEL 1: ECFP4 ONLY
# =============================================================================

print("\n" + "="*60)
print("MODEL 1: ECFP4 FINGERPRINTS ONLY")
print("="*60)

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# Random Forest (commonly used for fingerprints)
rf_model = RandomForestClassifier(n_estimators=100, max_depth=10, random_state=42, n_jobs=-1)
y_pred_rf = cross_val_predict(rf_model, X_fp, y, cv=cv, method='predict_proba')[:, 1]
auc_rf = roc_auc_score(y, y_pred_rf)
print(f"Random Forest AUC: {auc_rf:.3f}")

# XGBoost
if XGBOOST_AVAILABLE:
    xgb_model = XGBClassifier(n_estimators=100, max_depth=5, random_state=42, eval_metric='logloss')
    y_pred_xgb = cross_val_predict(xgb_model, X_fp, y, cv=cv, method='predict_proba')[:, 1]
    auc_xgb = roc_auc_score(y, y_pred_xgb)
    print(f"XGBoost AUC: {auc_xgb:.3f}")

# =============================================================================
# MODEL 2: ECFP4 + TABULAR FEATURES
# =============================================================================

print("\n" + "="*60)
print("MODEL 2: ECFP4 + TABULAR FEATURES")
print("="*60)

# Combine fingerprints with our existing features
tabular_features = ['logp', 'mw', 'hbd', 'hba', 'tpsa', 'k', 'ecnp_z', 
                    'has_network_targets', 'dili_fraction', 'rule_of_two']

X_tab = df_valid[[f for f in tabular_features if f in df_valid.columns]].values

# Standardize tabular features
scaler = StandardScaler()
X_tab_scaled = scaler.fit_transform(X_tab)

# Combine
X_combined = np.hstack([X_fp, X_tab_scaled])
print(f"Combined features: {X_combined.shape[1]} ({X_fp.shape[1]} FP + {X_tab_scaled.shape[1]} tab)")

# Random Forest
y_pred_combined_rf = cross_val_predict(rf_model, X_combined, y, cv=cv, method='predict_proba')[:, 1]
auc_combined_rf = roc_auc_score(y, y_pred_combined_rf)
print(f"Random Forest AUC: {auc_combined_rf:.3f}")

# XGBoost
if XGBOOST_AVAILABLE:
    y_pred_combined_xgb = cross_val_predict(xgb_model, X_combined, y, cv=cv, method='predict_proba')[:, 1]
    auc_combined_xgb = roc_auc_score(y, y_pred_combined_xgb)
    print(f"XGBoost AUC: {auc_combined_xgb:.3f}")

# =============================================================================
# MODEL 3: TABULAR ONLY (BASELINE)
# =============================================================================

print("\n" + "="*60)
print("MODEL 3: TABULAR FEATURES ONLY (BASELINE)")
print("="*60)

y_pred_tab_rf = cross_val_predict(rf_model, X_tab_scaled, y, cv=cv, method='predict_proba')[:, 1]
auc_tab_rf = roc_auc_score(y, y_pred_tab_rf)
print(f"Random Forest AUC: {auc_tab_rf:.3f}")

if XGBOOST_AVAILABLE:
    y_pred_tab_xgb = cross_val_predict(xgb_model, X_tab_scaled, y, cv=cv, method='predict_proba')[:, 1]
    auc_tab_xgb = roc_auc_score(y, y_pred_tab_xgb)
    print(f"XGBoost AUC: {auc_tab_xgb:.3f}")

# =============================================================================
# COMPARE TO LITERATURE
# =============================================================================

print("\n" + "="*60)
print("COMPARISON TO LITERATURE CLAIMS")
print("="*60)

print(f"""
Literature claims:
  CNN + ECFP4 (1,597 cpds):  AUC = 0.96
  Deep Learning (475 cpds): AUC = 0.955
  
Our results (202 cpds):
  ECFP4 only (RF):          AUC = {auc_rf:.3f}
  ECFP4 only (XGB):         AUC = {auc_xgb if XGBOOST_AVAILABLE else 'N/A':.3f}
  Combined (RF):            AUC = {auc_combined_rf:.3f}
  Combined (XGB):           AUC = {auc_combined_xgb if XGBOOST_AVAILABLE else 'N/A':.3f}
  Tabular only (XGB):       AUC = {auc_tab_xgb if XGBOOST_AVAILABLE else auc_tab_rf:.3f}
  
Previous best:              AUC = 0.800
""")

# Best result
best_auc = max(auc_rf, auc_combined_rf, auc_tab_rf)
if XGBOOST_AVAILABLE:
    best_auc = max(best_auc, auc_xgb, auc_combined_xgb, auc_tab_xgb)

print(f"Best AUC achieved: {best_auc:.3f}")
if best_auc >= 0.90:
    print("  -> Approaching literature claims!")
elif best_auc >= 0.85:
    print("  -> Good improvement, but not at literature level")
else:
    print("  -> Limited improvement - data size may be the bottleneck")

# =============================================================================
# SAVE RESULTS
# =============================================================================

df_valid['y_pred_ecfp'] = y_pred_rf
df_valid['y_pred_combined'] = y_pred_combined_rf
output = ROOT / 'research' / 'ecnp-generalization' / 'results' / 'ecfp_model_results.csv'
df_valid.to_csv(output, index=False)
print(f"\nSaved: {output}")

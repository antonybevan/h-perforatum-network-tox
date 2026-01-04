"""
Train Full Model on 706-Drug ECNP Dataset
==========================================

Uses same pipeline as original 202-drug model:
- Chemistry features (logP, MW, TPSA, etc.)
- Target topology features
- ECNP Z-score
- Compare to original 202-drug performance (AUC ~0.87)
"""
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
import warnings
warnings.filterwarnings('ignore')

try:
    from xgboost import XGBClassifier
    XGBOOST_AVAILABLE = True
except ImportError:
    XGBOOST_AVAILABLE = False

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("WARNING: RDKit not available - chemistry features will be missing")

ROOT = Path(r'v:\new\h-perforatum-network-tox')

print("="*70)
print("706-DRUG ECNP MODEL - FULL PIPELINE")
print("="*70)

# =============================================================================
# 1. LOAD DATA
# =============================================================================
print("\n--- Loading 706-Drug Dataset ---")
df = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'results' / 'dilirank_706_with_ecnp.csv')
print(f"Total drugs: {len(df)}")
print(f"With valid ECNP: {df['ecnp_z'].notna().sum()}")

# =============================================================================
# 2. COMPUTE CHEMISTRY FEATURES
# =============================================================================
print("\n--- Computing Chemistry Features ---")

if RDKIT_AVAILABLE:
    def get_chem(smiles):
        try:
            mol = Chem.MolFromSmiles(smiles) if pd.notna(smiles) else None
            if mol:
                return {
                    'logp': Descriptors.MolLogP(mol),
                    'mw': Descriptors.MolWt(mol),
                    'hbd': Lipinski.NumHDonors(mol),
                    'hba': Lipinski.NumHAcceptors(mol),
                    'tpsa': Descriptors.TPSA(mol),
                    'rotatable': Descriptors.NumRotatableBonds(mol)
                }
        except:
            pass
        return {k: np.nan for k in ['logp', 'mw', 'hbd', 'hba', 'tpsa', 'rotatable']}
    
    chem = df['smiles'].apply(get_chem).apply(pd.Series)
    df = pd.concat([df, chem], axis=1)
    print(f"Chemistry features computed for {chem['logp'].notna().sum()} drugs")
else:
    df['logp'] = np.nan
    df['mw'] = np.nan

# =============================================================================
# 3. FEATURE ENGINEERING
# =============================================================================
print("\n--- Engineering Features ---")

# Target features
df['k'] = df['n_targets_mapped']
df['high_logp'] = (df['logp'] > 3).astype(int) if 'logp' in df.columns else 0

# =============================================================================
# 4. ABLATION STUDY
# =============================================================================
print("\n" + "="*70)
print("ABLATION STUDY")
print("="*70)

# Filter to valid samples
feature_cols_ml = ['logp', 'mw', 'hbd', 'hba', 'tpsa', 'k']
feature_cols_full = feature_cols_ml + ['ecnp_z']

df_valid = df[df[feature_cols_full].notna().all(axis=1)].copy()
print(f"\nValid samples for full model: {len(df_valid)}")
print(f"  DILI+: {df_valid['is_dili'].sum()}")
print(f"  DILI-: {(df_valid['is_dili']==0).sum()}")

X_ml = df_valid[feature_cols_ml].values
X_full = df_valid[feature_cols_full].values  
y = df_valid['is_dili'].values

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
results = {}

# Model 1: ML only (no ECNP)
print("\n1. ML WITHOUT ECNP (chemistry + targets)")
X_ml_scaled = StandardScaler().fit_transform(X_ml)

if XGBOOST_AVAILABLE:
    model_ml = XGBClassifier(n_estimators=100, max_depth=3, random_state=42, eval_metric='logloss', verbosity=0)
else:
    model_ml = LogisticRegression(max_iter=1000, random_state=42)

y_pred_ml = cross_val_predict(model_ml, X_ml_scaled, y, cv=cv, method='predict_proba')[:, 1]
auc_ml = roc_auc_score(y, y_pred_ml)
print(f"   AUC = {auc_ml:.3f}")
results['ML without ECNP'] = auc_ml

# Model 2: ECNP alone
print("\n2. ECNP ALONE")
auc_ecnp = roc_auc_score(y, df_valid['ecnp_z'].values)
print(f"   AUC = {auc_ecnp:.3f}")
results['ECNP alone'] = auc_ecnp

# Model 3: Full model (ML + ECNP)
print("\n3. FULL MODEL (ML + ECNP)")
X_full_scaled = StandardScaler().fit_transform(X_full)

if XGBOOST_AVAILABLE:
    model_full = XGBClassifier(n_estimators=100, max_depth=3, random_state=42, eval_metric='logloss', verbosity=0)
else:
    model_full = LogisticRegression(max_iter=1000, random_state=42)

y_pred_full = cross_val_predict(model_full, X_full_scaled, y, cv=cv, method='predict_proba')[:, 1]
auc_full = roc_auc_score(y, y_pred_full)
print(f"   AUC = {auc_full:.3f}")
results['Full model (ML + ECNP)'] = auc_full

# Model 4: LogP baseline
print("\n4. LogP ALONE (baseline)")
auc_logp = roc_auc_score(y, df_valid['logp'].values)
print(f"   AUC = {auc_logp:.3f}")
results['LogP alone'] = auc_logp

# =============================================================================
# 5. COMPARISON TO ORIGINAL 202-DRUG MODEL
# =============================================================================
print("\n" + "="*70)
print("COMPARISON TO ORIGINAL 202-DRUG MODEL")
print("="*70)

original_auc = 0.870  # From two_tier_training.txt
print(f"""
Original 202-Drug Model:
  - AUC: {original_auc:.3f}
  - Drugs: 202 (with full ECNP coverage)

This 706-Drug Model:
  - AUC: {auc_full:.3f}
  - Drugs: {len(df_valid)} (with valid ECNP + chemistry)

Delta: {auc_full - original_auc:+.3f}
""")

# =============================================================================
# 6. SUMMARY
# =============================================================================
print("="*70)
print("SUMMARY: ABLATION RESULTS")
print("="*70)

for name, auc in sorted(results.items(), key=lambda x: x[1], reverse=True):
    print(f"  {name:25s}: AUC = {auc:.3f}")

ecnp_contrib = results['Full model (ML + ECNP)'] - results['ML without ECNP']
print(f"\nECNP contribution: {ecnp_contrib:+.3f}")

# Save results
output = ROOT / 'research' / 'ecnp-generalization' / 'results' / 'model_706_vs_202.csv'
pd.DataFrame({
    'model': list(results.keys()) + ['Original 202-drug'],
    'auc': list(results.values()) + [original_auc],
    'dataset': ['706'] * len(results) + ['202']
}).to_csv(output, index=False)
print(f"\nSaved: {output}")

"""
Two-Tier DILI Model - Full Training
====================================

Tier 1: Chemistry-only on 901 drugs (full DILIrank with SMILES)
Tier 2: Chemistry + ECNP on 202 drugs (network subset)

Compare performance and document for publication.
"""
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
    RDKIT_AVAILABLE = True
except:
    RDKIT_AVAILABLE = False
    print("RDKit not available!")
    exit()

ROOT = Path(r'v:\new\h-perforatum-network-tox')

print("="*70)
print("TWO-TIER DILI MODEL TRAINING")
print("="*70)

# =============================================================================
# LOAD TIER 1 DATA (901 drugs)
# =============================================================================

print("\n--- Loading Tier 1: Full DILIrank (901 drugs) ---")

df_full = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'results' / 'dilirank_full_smiles.csv')
print(f"Loaded: {len(df_full)} drugs")
print(f"DILI+: {df_full['is_dili'].sum()}, DILI-: {(df_full['is_dili'] == 0).sum()}")

# =============================================================================
# GENERATE ECFP + DESCRIPTORS FOR TIER 1
# =============================================================================

print("\n--- Generating ECFP4 + Physicochemical Features ---")

def get_features(smiles):
    """Generate ECFP4 + physicochemical descriptors."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        # ECFP4 (radius=2)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
        fp_array = np.array(fp)
        
        # Physicochemical
        logp = Descriptors.MolLogP(mol)
        mw = Descriptors.MolWt(mol)
        tpsa = Descriptors.TPSA(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        
        return {
            'fp': fp_array,
            'logp': logp,
            'mw': mw,
            'tpsa': tpsa,
            'hbd': hbd,
            'hba': hba
        }
    except:
        return None

# Process Tier 1
fps_tier1 = []
desc_tier1 = []
valid_idx_tier1 = []

for idx, row in df_full.iterrows():
    result = get_features(row['smiles'])
    if result:
        fps_tier1.append(result['fp'])
        desc_tier1.append({k: v for k, v in result.items() if k != 'fp'})
        valid_idx_tier1.append(idx)

X_fp_tier1 = np.array(fps_tier1)
df_tier1 = df_full.loc[valid_idx_tier1].copy()
desc_df = pd.DataFrame(desc_tier1, index=valid_idx_tier1)
for col in desc_df.columns:
    df_tier1[col] = desc_df[col]

y_tier1 = df_tier1['is_dili'].values
print(f"Valid Tier 1 compounds: {len(df_tier1)}")

# =============================================================================
# TRAIN TIER 1: CHEMISTRY ONLY
# =============================================================================

print("\n" + "="*70)
print("TIER 1: CHEMISTRY-ONLY (ECFP4 + LogP/MW/TPSA/HBD/HBA)")
print("="*70)

tabular_cols = ['logp', 'mw', 'tpsa', 'hbd', 'hba']
X_tab_tier1 = df_tier1[tabular_cols].values

scaler = StandardScaler()
X_tab_tier1_scaled = scaler.fit_transform(X_tab_tier1)
X_tier1 = np.hstack([X_fp_tier1, X_tab_tier1_scaled])

print(f"Training samples: {len(y_tier1)}")
print(f"Features: {X_tier1.shape[1]} (ECFP: 1024, Tabular: {len(tabular_cols)})")
print(f"Class balance: DILI+ {y_tier1.mean():.1%}, DILI- {1-y_tier1.mean():.1%}")

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
model = RandomForestClassifier(n_estimators=100, max_depth=10, random_state=42, n_jobs=-1)

y_pred_tier1 = cross_val_predict(model, X_tier1, y_tier1, cv=cv, method='predict_proba')[:, 1]
auc_tier1 = roc_auc_score(y_tier1, y_pred_tier1)

print(f"\n*** TIER 1 AUC: {auc_tier1:.3f} ***")

# =============================================================================
# LOAD TIER 2 DATA (202 drugs with ECNP)
# =============================================================================

print("\n" + "="*70)
print("TIER 2: CHEMISTRY + ECNP (Network Subset)")
print("="*70)

df_tier2 = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'results' / 'ecfp_model_results.csv')
print(f"Loaded: {len(df_tier2)} drugs with ECNP")

# Process Tier 2
fps_tier2 = []
valid_idx_tier2 = []

for idx, row in df_tier2.iterrows():
    result = get_features(row['smiles'])
    if result:
        fps_tier2.append(result['fp'])
        valid_idx_tier2.append(idx)

X_fp_tier2 = np.array(fps_tier2)
df_tier2_valid = df_tier2.loc[valid_idx_tier2].copy()
y_tier2 = df_tier2_valid['is_dili'].values

# Features with ECNP
tabular_cols_ecnp = ['logp', 'mw', 'hbd', 'hba', 'tpsa', 'ecnp_z']
tabular_cols_ecnp = [c for c in tabular_cols_ecnp if c in df_tier2_valid.columns]

X_tab_tier2 = df_tier2_valid[tabular_cols_ecnp].values
X_tab_tier2_scaled = scaler.fit_transform(X_tab_tier2)
X_tier2 = np.hstack([X_fp_tier2, X_tab_tier2_scaled])

print(f"Training samples: {len(y_tier2)}")
print(f"Features: {X_tier2.shape[1]} (ECFP: 1024, Tabular+ECNP: {len(tabular_cols_ecnp)})")

y_pred_tier2 = cross_val_predict(model, X_tier2, y_tier2, cv=cv, method='predict_proba')[:, 1]
auc_tier2 = roc_auc_score(y_tier2, y_pred_tier2)

print(f"\n*** TIER 2 AUC: {auc_tier2:.3f} ***")

# =============================================================================
# TIER 1 ON TIER 2 SUBSET (for fair comparison)
# =============================================================================

print("\n" + "="*70)
print("TIER 1 (CHEMISTRY-ONLY) ON TIER 2 SUBSET (Fair Comparison)")
print("="*70)

# Remove ECNP from Tier 2 features
tabular_cols_no_ecnp = [c for c in tabular_cols_ecnp if c != 'ecnp_z']
X_tab_tier2_no_ecnp = df_tier2_valid[tabular_cols_no_ecnp].values
X_tab_tier2_no_ecnp_scaled = scaler.fit_transform(X_tab_tier2_no_ecnp)
X_tier2_no_ecnp = np.hstack([X_fp_tier2, X_tab_tier2_no_ecnp_scaled])

y_pred_tier2_no_ecnp = cross_val_predict(model, X_tier2_no_ecnp, y_tier2, cv=cv, method='predict_proba')[:, 1]
auc_tier2_no_ecnp = roc_auc_score(y_tier2, y_pred_tier2_no_ecnp)

ecnp_contribution = auc_tier2 - auc_tier2_no_ecnp

print(f"Without ECNP: AUC = {auc_tier2_no_ecnp:.3f}")
print(f"With ECNP: AUC = {auc_tier2:.3f}")
print(f"ECNP contribution: {ecnp_contribution:+.3f}")

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "="*70)
print("FINAL SUMMARY")
print("="*70)

print(f"""
MODEL COMPARISON:

  Tier 1 (Full DILIrank - Chemistry Only)
    Drugs:    {len(y_tier1)}
    Features: ECFP4 + LogP/MW/TPSA/HBD/HBA
    AUC:      {auc_tier1:.3f}

  Tier 2 (Network Subset - Chemistry + ECNP)
    Drugs:    {len(y_tier2)}
    Features: ECFP4 + LogP/MW/TPSA/HBD/HBA + ECNP Z-score
    AUC:      {auc_tier2:.3f}

  Fair Comparison (Tier 2 subset, chemistry-only):
    AUC without ECNP: {auc_tier2_no_ecnp:.3f}
    AUC with ECNP:    {auc_tier2:.3f}
    ECNP lift:        {ecnp_contribution:+.3f}

KEY FINDINGS:
  - Tier 1 (901 drugs): AUC {auc_tier1:.3f}
  - Tier 2 (202 drugs): AUC {auc_tier2:.3f}  
  - ECNP adds {ecnp_contribution:+.3f} AUC on network subset
""")

# Save results
results = {
    'tier1_drugs': len(y_tier1),
    'tier1_auc': auc_tier1,
    'tier2_drugs': len(y_tier2),
    'tier2_auc_with_ecnp': auc_tier2,
    'tier2_auc_no_ecnp': auc_tier2_no_ecnp,
    'ecnp_contribution': ecnp_contribution
}
pd.DataFrame([results]).to_csv(
    ROOT / 'research' / 'ecnp-generalization' / 'results' / 'two_tier_final.csv',
    index=False
)
print(f"Saved: two_tier_final.csv")

# Save Tier 1 predictions
df_tier1['y_pred'] = y_pred_tier1
df_tier1.to_csv(
    ROOT / 'research' / 'ecnp-generalization' / 'results' / 'tier1_predictions.csv',
    index=False
)
print(f"Saved: tier1_predictions.csv")

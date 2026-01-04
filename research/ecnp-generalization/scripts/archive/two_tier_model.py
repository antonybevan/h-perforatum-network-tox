"""
Two-Tier DILI Prediction Model
==============================

Tier 1: Chemistry-first (Full DILIrank)
- ECFP4 fingerprints + PhysChem features
- "Is this compound risky at all?"

Tier 2: Network-augmented refinement (706-drug subset)
- Tier 1 features + ECNP Z-score
- "Is this risk driven by biologically coherent target positioning?"

ECNP as conditional enhancer, not universal requirement.
"""
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
import warnings
warnings.filterwarnings('ignore')

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, Lipinski
    RDKIT_AVAILABLE = True
except:
    RDKIT_AVAILABLE = False
    print("ERROR: RDKit required for ECFP")
    exit(1)

ROOT = Path(r'v:\new\h-perforatum-network-tox')

print("="*70)
print("TWO-TIER DILI MODEL WITH ECNP")
print("="*70)

# =============================================================================
# LOAD DATA
# =============================================================================
print("\n--- Loading Data ---")

# DILIrank full (for Tier 1)
dilirank = pd.read_csv(ROOT / 'research/ecnp-generalization/results/dilirank_full_smiles.csv')
print(f"DILIrank full: {len(dilirank)} drugs")

# 706-drug subset with ECNP (for Tier 2)
df_706 = pd.read_csv(ROOT / 'research/ecnp-generalization/results/dilirank_706_with_ecnp.csv')
print(f"706-drug ECNP subset: {len(df_706)} drugs")

# =============================================================================
# FEATURE FUNCTIONS
# =============================================================================

def compute_ecfp(smiles, radius=2, bits=1024):
    """Compute ECFP4 fingerprint"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=bits)
            return list(fp)
    except:
        pass
    return [np.nan] * bits

def compute_physchem(smiles):
    """Compute physicochemical properties"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return {
                'logp': Descriptors.MolLogP(mol),
                'mw': Descriptors.MolWt(mol),
                'tpsa': Descriptors.TPSA(mol),
                'hbd': Lipinski.NumHDonors(mol),
                'hba': Lipinski.NumHAcceptors(mol),
            }
    except:
        pass
    return {k: np.nan for k in ['logp', 'mw', 'tpsa', 'hbd', 'hba']}

# =============================================================================
# TIER 1: CHEMISTRY-FIRST MODEL (Full DILIrank)
# =============================================================================
print("\n" + "="*70)
print("TIER 1: CHEMISTRY-FIRST MODEL")
print("="*70)

# Compute features for all DILIrank
print("Computing ECFP + PhysChem for DILIrank...")
ecfp_list = []
physchem_list = []

for smiles in dilirank['smiles']:
    ecfp_list.append(compute_ecfp(smiles))
    physchem_list.append(compute_physchem(smiles))

ecfp_df = pd.DataFrame(ecfp_list, columns=[f'ecfp_{i}' for i in range(1024)])
physchem_df = pd.DataFrame(physchem_list)

tier1_df = pd.concat([dilirank[['dilirank_name', 'is_dili']], ecfp_df, physchem_df], axis=1)

# Filter valid rows
tier1_valid = tier1_df.dropna()
print(f"Valid Tier 1 samples: {len(tier1_valid)}")
print(f"  DILI+: {tier1_valid.is_dili.sum()}, DILI-: {(tier1_valid.is_dili==0).sum()}")

# Feature matrix
feature_cols_t1 = [c for c in tier1_valid.columns if c.startswith('ecfp_') or c in ['logp','mw','tpsa','hbd','hba']]
X_t1 = tier1_valid[feature_cols_t1].values
y_t1 = tier1_valid['is_dili'].values

# Train/Evaluate
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
model_t1 = LogisticRegression(class_weight='balanced', max_iter=1000, random_state=42)
y_pred_t1 = cross_val_predict(model_t1, X_t1, y_t1, cv=cv, method='predict_proba')[:, 1]

auc_t1 = roc_auc_score(y_t1, y_pred_t1)
pr_auc_t1 = average_precision_score(y_t1, y_pred_t1)

print(f"\nTIER 1 RESULTS:")
print(f"  ROC-AUC: {auc_t1:.3f}")
print(f"  PR-AUC:  {pr_auc_t1:.3f}")

# Store predictions
tier1_valid = tier1_valid.copy()
tier1_valid['tier1_score'] = y_pred_t1

# =============================================================================
# TIER 2: NETWORK-AUGMENTED MODEL (706-drug subset)
# =============================================================================
print("\n" + "="*70)
print("TIER 2: NETWORK-AUGMENTED MODEL (ECNP)")
print("="*70)

# Compute features for 706 subset
print("Computing ECFP + PhysChem for 706-drug subset...")
ecfp_list = []
physchem_list = []

for smiles in df_706['smiles']:
    ecfp_list.append(compute_ecfp(smiles))
    physchem_list.append(compute_physchem(smiles))

ecfp_df = pd.DataFrame(ecfp_list, columns=[f'ecfp_{i}' for i in range(1024)])
physchem_df = pd.DataFrame(physchem_list)

tier2_df = pd.concat([
    df_706[['dilirank_name', 'is_dili', 'ecnp_z', 'n_targets_mapped', 'pool_size']],
    ecfp_df, 
    physchem_df
], axis=1)

# Filter: need valid ECNP
tier2_valid = tier2_df[tier2_df['ecnp_z'].notna()].dropna()
print(f"Valid Tier 2 samples: {len(tier2_valid)}")
print(f"  DILI+: {tier2_valid.is_dili.sum()}, DILI-: {(tier2_valid.is_dili==0).sum()}")

# Features WITHOUT ECNP (for ablation)
feature_cols_base = [c for c in tier2_valid.columns if c.startswith('ecfp_') or c in ['logp','mw','tpsa','hbd','hba']]

# Features WITH ECNP
feature_cols_ecnp = feature_cols_base + ['ecnp_z']

X_base = tier2_valid[feature_cols_base].values
X_ecnp = tier2_valid[feature_cols_ecnp].values
y_t2 = tier2_valid['is_dili'].values

# Train/Evaluate WITHOUT ECNP
model_base = LogisticRegression(class_weight='balanced', max_iter=1000, random_state=42)
y_pred_base = cross_val_predict(model_base, X_base, y_t2, cv=cv, method='predict_proba')[:, 1]
auc_base = roc_auc_score(y_t2, y_pred_base)

# Train/Evaluate WITH ECNP
model_ecnp = LogisticRegression(class_weight='balanced', max_iter=1000, random_state=42)
y_pred_ecnp = cross_val_predict(model_ecnp, X_ecnp, y_t2, cv=cv, method='predict_proba')[:, 1]
auc_ecnp = roc_auc_score(y_t2, y_pred_ecnp)
pr_auc_ecnp = average_precision_score(y_t2, y_pred_ecnp)

delta = auc_ecnp - auc_base

print(f"\nTIER 2 RESULTS:")
print(f"  WITHOUT ECNP: AUC = {auc_base:.3f}")
print(f"  WITH ECNP:    AUC = {auc_ecnp:.3f}")
print(f"  ECNP lift:    {delta:+.3f}")
print(f"  PR-AUC:       {pr_auc_ecnp:.3f}")

tier2_valid = tier2_valid.copy()
tier2_valid['tier2_score'] = y_pred_ecnp

# =============================================================================
# SUMMARY
# =============================================================================
print("\n" + "="*70)
print("SUMMARY")
print("="*70)
print(f"""
TIER 1: Chemistry-First (General Screen)
  - Drugs: {len(tier1_valid)}
  - Features: ECFP4 + LogP/MW/TPSA/HBD/HBA
  - ROC-AUC: {auc_t1:.3f}

TIER 2: Network-Augmented (Mechanistic Refinement)
  - Drugs: {len(tier2_valid)}
  - Features: Tier 1 + ECNP Z-score
  - ROC-AUC without ECNP: {auc_base:.3f}
  - ROC-AUC with ECNP:    {auc_ecnp:.3f}
  - ECNP contribution:    {delta:+.3f}
""")

# Save results
results = pd.DataFrame({
    'tier': ['Tier 1', 'Tier 2 (no ECNP)', 'Tier 2 (with ECNP)'],
    'n_drugs': [len(tier1_valid), len(tier2_valid), len(tier2_valid)],
    'auc': [auc_t1, auc_base, auc_ecnp],
    'ecnp_contribution': [None, None, delta]
})
output = ROOT / 'research/ecnp-generalization/results/two_tier_model_results.csv'
results.to_csv(output, index=False)
print(f"Saved: {output}")

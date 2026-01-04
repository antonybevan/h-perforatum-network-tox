"""
Script: test_ecnp_mechanism.py
Purpose: Validate the hypothesis that ECNP improves discrimination specifically for toxicity mediated by network-proximal pathways.

Methodology:
1. Load Gold Standard Data (202 drugs, Tier 2b).
2. Compute baseline AUC (Chemistry Only) vs Enhanced AUC (+ECNP).
3. Calculate 'Lift' in stratified subgroups:
    - Direct Hitters (n_dili_hits > 0) vs Indirect (n_dili_hits == 0)
    - Transporter Dominant vs Cascade Dominant
    - High vs Low Target Count
"""
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
import warnings
warnings.filterwarnings('ignore')

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Lipinski

ROOT = Path(r'v:\new\h-perforatum-network-tox')
DATA_PATH = ROOT / 'research/ecnp-generalization/results/ecfp_model_results.csv'

# =============================================================================
# 1. LOAD DATA & FEATURES
# =============================================================================
print("Loading data...")
df = pd.read_csv(DATA_PATH)

# Ensure ECNP column exists
ecnp_col = 'ecnp_z' if 'ecnp_z' in df.columns else 'Z'
if ecnp_col not in df.columns:
    raise ValueError(f"ECNP column not found in {df.columns}")

# COMPUTE ECFP (if not present or to be safe)
print("Computing chemical features...")
ecfp = []
physchem = []
for s in df['smiles']:
    mol = Chem.MolFromSmiles(s)
    if mol:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
        ecfp.append(list(fp))
        physchem.append([
            Descriptors.MolLogP(mol),
            Descriptors.MolWt(mol),
            Descriptors.TPSA(mol),
            Lipinski.NumHDonors(mol),
            Lipinski.NumHAcceptors(mol)
        ])
    else:
        ecfp.append([0]*1024)
        physchem.append([0]*5)

X_chem = np.hstack([np.array(ecfp), np.array(physchem)])
X_ecnp = df[[ecnp_col]].values.reshape(-1, 1)
y = df['is_dili'].values

# =============================================================================
# 2. DEFINING SUBGROUPS
# =============================================================================
print("Defining subgroups...")

# A. DIRECT VS INDIRECT (The 'Hidden Trap')
# Direct: Hits known DILI genes directly
# Indirect: No direct DILI hits (but is_dili=1 or 0)
mask_direct = df['n_dili_hits'] > 0
mask_indirect = df['n_dili_hits'] == 0

# B. MECHANISM TYPE
# Transporter: Hits transporters but NO other targets (very specific)
# Cascade: Hits targets but NOT transporters (signaling)
mask_transporter = (df['n_transporter_hits'] > 0) & (df['n_transporter_hits'] == df['n_dili_hits'])
mask_cascade = (df['n_targets'] > 0) & (df['n_transporter_hits'] == 0)

# C. SPECIFICITY
# Specific: <= 2 targets
# Promiscuous: >= 5 targets
mask_specific = df['n_targets'] <= 2
mask_promiscuous = df['n_targets'] >= 5

subgroups = {
    'All Drugs': slice(None),
    'Direct DILI Hitters': mask_direct,
    'Indirect (Hidden) Risk': mask_indirect,
    'Transporter Dominant': mask_transporter,
    'Cascade (Signaling) Dominant': mask_cascade,
    'High Specificity': mask_specific,
    'Promiscuous (Dirty)': mask_promiscuous
}

# =============================================================================
# 3. EVALUATION FUNCTION
# =============================================================================
def get_auc(X_subset, y_subset):
    if len(np.unique(y_subset)) < 2:
        return np.nan
    
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    model = LogisticRegression(class_weight='balanced', max_iter=1000)
    
    try:
        y_pred = cross_val_predict(model, X_subset, y_subset, cv=cv, method='predict_proba')[:, 1]
        return roc_auc_score(y_subset, y_pred)
    except:
        return np.nan

# =============================================================================
# 4. RUN ANALYSIS
# =============================================================================
results = []
print("\nRUNNING STRATIFIED ANALYSIS")
print(f"{'Subgroup':<30} | {'N':<5} | {'Base AUC':<8} | {'+ECNP AUC':<9} | {'Lift':<6}")
print("-" * 75)

for name, mask in subgroups.items():
    # Filter data
    y_sub = y[mask]
    X_chem_sub = X_chem[mask]
    X_ecnp_sub = np.hstack([X_chem_sub, X_ecnp[mask]])
    
    n = len(y_sub)
    if n < 10: # Skip tiny groups
        continue
        
    auc_base = get_auc(X_chem_sub, y_sub)
    auc_ecnp = get_auc(X_ecnp_sub, y_sub)
    
    lift = auc_ecnp - auc_base
    
    results.append({
        'Subgroup': name,
        'N': n,
        'Base_AUC': auc_base,
        'ECNP_AUC': auc_ecnp,
        'Lift': lift
    })
    
    print(f"{name:<30} | {n:<5} | {auc_base:.3f}    | {auc_ecnp:.3f}     | {lift:+.3f}")

# Save
df_res = pd.DataFrame(results)
df_res.to_csv(ROOT / 'research/ecnp-generalization/results/mechanistic_test_results.csv', index=False)
print("\nSaved to results/mechanistic_test_results.csv")

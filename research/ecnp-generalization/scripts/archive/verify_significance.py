"""
Script: verify_significance.py
Purpose: 
1. Calculate statistical significance (p-value) of the ECNP lift using Bootstrapping.
2. Identify "Rescued Drugs": Specific compounds where Chemistry failed but ECNP succeeded.

Context: AUC 0.871 -> 0.875 (+0.004) might seem small, but is it statistically significant?
"""
import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Lipinski
from pathlib import Path

ROOT = Path(r'v:\new\h-perforatum-network-tox')
DATA_PATH = ROOT / 'research/ecnp-generalization/results/dilirank_706_with_ecnp.csv'

# =============================================================================
# 1. LOAD & PREP DATA (Tier 2a Expanded)
# =============================================================================
print("Loading Tier 2a Expanded data...")
df = pd.read_csv(DATA_PATH)
# Filter for valid ECNP
df = df[df['ecnp_z'].notna()].copy()
print(f"N = {len(df)}")

# Re-generate clean features to match exactly what we just trained
print("Regenerating features...")
chem_features = []
for s in df['smiles']:
    mol = Chem.MolFromSmiles(s)
    if mol:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
        chem_features.append(list(fp) + [
            Descriptors.MolLogP(mol),
            Descriptors.MolWt(mol),
            Descriptors.TPSA(mol)
        ])
    else:
        chem_features.append([0]*1027)

X_chem = np.array(chem_features)

# Decomposed ECNP Features
X_ecnp = []

if 'network_fraction' not in df.columns:
    # Compute for Tier 2a
    if 'n_targets_uniprot' in df.columns:
        df['network_fraction'] = df['k'] / df['n_targets_uniprot'].replace(0, 1)
    else:
        df['network_fraction'] = 0

for idx, row in df.iterrows():
    it = row.get('I_T', 0)
    # Tier 2a has 'ecnp_mu_T', Tier 2b has 'mu_T'
    mu = row.get('ecnp_mu_T', row.get('mu_T', 0))
    k = row.get('k', 0)
    pool = row.get('pool_size', 0)
    frac = row.get('network_fraction', 0)
    
    # Derive sigma
    sigma = row.get('ecnp_sigma_T', row.get('sigma_T', 0))
    if pd.isna(sigma) or sigma == 0:
        z = row.get('ecnp_z', row.get('Z', 0))
        if not pd.isna(z) and z != 0:
            sigma = (it - mu) / z
        else:
            sigma = 0
            
    X_ecnp.append([it, mu, sigma, k, pool, frac])

X_ecnp = np.nan_to_num(np.array(X_ecnp))
X_full = np.hstack([X_chem, X_ecnp])
y = df['is_dili'].values

# =============================================================================
# 2. GET PREDICTIONS (5-Fold CV)
# =============================================================================
print("Running Cross-Validation...")
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

model_chem = make_pipeline(StandardScaler(), LogisticRegression(class_weight='balanced', max_iter=1000, C=0.1))
model_full = make_pipeline(StandardScaler(), LogisticRegression(class_weight='balanced', max_iter=1000, C=0.1))

# Get probas
y_pred_chem = cross_val_predict(model_chem, X_chem, y, cv=cv, method='predict_proba')[:, 1]
y_pred_full = cross_val_predict(model_full, X_full, y, cv=cv, method='predict_proba')[:, 1]

auc_chem = roc_auc_score(y, y_pred_chem)
auc_full = roc_auc_score(y, y_pred_full)
print(f"Base AUC: {auc_chem:.4f}")
print(f"Full AUC: {auc_full:.4f}")
print(f"Lift:     {auc_full - auc_chem:.4f}")

# =============================================================================
# 3. BOOTSTRAP SIGNIFICANCE TEST
# =============================================================================
print("\nRunning Bootstrap Test (n=1000)...")
n_bootstraps = 1000
rng = np.random.RandomState(42)
boot_diffs = []

for i in range(n_bootstraps):
    # Bootstrap indices
    indices = rng.randint(0, len(y), len(y))
    if len(np.unique(y[indices])) < 2:
        continue
        
    auc_b_chem = roc_auc_score(y[indices], y_pred_chem[indices])
    auc_b_full = roc_auc_score(y[indices], y_pred_full[indices])
    boot_diffs.append(auc_b_full - auc_b_chem)

boot_diffs = np.array(boot_diffs)
p_value = (boot_diffs <= 0).mean()
ci_lower = np.percentile(boot_diffs, 2.5)
ci_upper = np.percentile(boot_diffs, 97.5)

print(f"p-value: {p_value:.4f}")
print(f"95% CI: [{ci_lower:.4f}, {ci_upper:.4f}]")

# =============================================================================
# 4. RESCUED DRUGS ANALYSIS
# =============================================================================
print("\n--- RESCUED DRUGS (False Negative -> True Positive) ---")
# Use a threshold of 0.5 for classification
df['pred_chem'] = (y_pred_chem > 0.5).astype(int)
df['pred_full'] = (y_pred_full > 0.5).astype(int)

# Identify Name Column
name_col = 'drug_name' if 'drug_name' in df.columns else 'dilirank_name'

# Identify cases where Full got it RIGHT and Chem got it WRONG
rescued = df[ (df['is_dili'] == 1) & (df['pred_chem'] == 0) & (df['pred_full'] == 1) ]
print(f"Total Rescued (FN -> TP): {len(rescued)}")
if len(rescued) > 0:
    print(rescued[[name_col, 'I_T', 'network_fraction', 'pred_chem', 'pred_full']].to_string(index=False))

print("\n--- SAVED DRUGS (False Positive -> True Negative) ---")
saved = df[ (df['is_dili'] == 0) & (df['pred_chem'] == 1) & (df['pred_full'] == 0) ]
print(f"Total Saved (FP -> TN): {len(saved)}")
if len(saved) > 0:
    print(saved[[name_col, 'I_T', 'network_fraction', 'pred_chem', 'pred_full']].to_string(index=False))

"""
Script: optimize_ecnp_features.py
Purpose: Test if decomposing ECNP into its raw components (I_T, mu_T, k) improves predictive performance compared to the single Z-score.

Hypothesis: The Z-score normalization ((I_T - mu_T) / sigma) might be removing valuable signal. "Expected" impact (mu_T) might still be toxic.

Models to compare:
1. Baseline: Chemistry only
2. Standard ECNP: Chemistry + Z
3. Decomposed ECNP: Chemistry + Z + I_T + mu_T + k + network_fraction
"""
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Lipinski

ROOT = Path(r'v:\new\h-perforatum-network-tox')
DATA_PATH = ROOT / 'research/ecnp-generalization/results/ecfp_model_results.csv'

# =============================================================================
# 1. LOAD DATA
# =============================================================================
print("Loading data...")
df = pd.read_csv(DATA_PATH)

# =============================================================================
# 2. FEATURE ENGINEERING
# =============================================================================
print("Generating features...")

# A. Chemistry (ECFP + PhysChem)
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

# B. ECNP Standard
# Prefer 'ecnp_z' if available, else 'Z'
z_col = 'ecnp_z' if 'ecnp_z' in df.columns else 'Z'
X_z = df[[z_col]].fillna(0).values

# C. ECNP Decomposed
# Keys: I_T (Impact), mu_T (Expected Impact), k (LCC Size), network_fraction
ecnp_cols = [z_col, 'I_T', 'mu_T', 'k', 'network_fraction', 'pool_size']
# Fill NaNs with 0
X_decomposed = df[ecnp_cols].fillna(0).values

y = df['is_dili'].values

# =============================================================================
# 3. EVALUATION
# =============================================================================
def evaluate(X, y, name):
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    # Use LR for consistency, but scale features first for regularization equity
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    model = LogisticRegression(class_weight='balanced', max_iter=1000, C=0.1) 
    y_pred = cross_val_predict(model, X_scaled, y, cv=cv, method='predict_proba')[:, 1]
    return roc_auc_score(y, y_pred)

print("\nRESULTS (5-Fold CV AUC):")
print("-" * 40)

# 1. Baseline
auc_base = evaluate(X_chem, y, "Chemistry Only")
print(f"Chemistry Only:        {auc_base:.4f}")

# 2. Standard ECNP
X_standard = np.hstack([X_chem, X_z])
auc_standard = evaluate(X_standard, y, "Standard ECNP (+Z)")
print(f"Standard ECNP (+Z):    {auc_standard:.4f} (Lift: {auc_standard-auc_base:+.4f})")

# 3. Decomposed ECNP
X_full = np.hstack([X_chem, X_decomposed])
auc_full = evaluate(X_full, y, "Decomposed ECNP")
print(f"Decomposed ECNP:       {auc_full:.4f} (Lift: {auc_full-auc_base:+.4f})")

# 4. RF Check (Non-linear interactions?)
print("\n--- Non-Linear Check (Random Forest) ---")
rf = RandomForestClassifier(n_estimators=100, class_weight='balanced', random_state=42, max_depth=5)
# Baseline RF
y_pred_rf_base = cross_val_predict(rf, X_chem, y, cv=5, method='predict_proba')[:, 1]
auc_rf_base = roc_auc_score(y, y_pred_rf_base)
print(f"RF Chem Only:          {auc_rf_base:.4f}")

# Full RF
y_pred_rf_full = cross_val_predict(rf, X_full, y, cv=5, method='predict_proba')[:, 1]
auc_rf_full = roc_auc_score(y, y_pred_rf_full)
print(f"RF Decomposed ECNP:    {auc_rf_full:.4f} (Lift: {auc_rf_full-auc_rf_base:+.4f})")

# =============================================================================
# 4. FEATURE IMPORTANCE (Coefficients)
# =============================================================================
print("\n--- Feature Importance (LR Coefficients from Full Model) ---")
# Fit simple LR on just the ECNP features to see directionality
X_ecnp_only = df[ecnp_cols].fillna(0).values
scaler = StandardScaler()
X_ecnp_scaled = scaler.fit_transform(X_ecnp_only)
model_simple = LogisticRegression(class_weight='balanced').fit(X_ecnp_scaled, y)

for name, coef in zip(ecnp_cols, model_simple.coef_[0]):
    print(f"{name:<20}: {coef:+.4f}")

"""
Pipeline Step 03: Model Training
================================
Trains the Trusted Multi-View Classifier (TMC) with Evidential Deep Learning.

Inputs:
    - results/tmc_features_706.csv
Outputs:
    - results/final_dili_predictions.csv

Usage:
    python pipeline/03_train_models.py
"""
import sys
sys.path.insert(0, str(__file__).replace('pipeline\\03_train_models.py', ''))

import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, average_precision_score

from src.config import TMC_FEATURES, RESULTS_DIR, RANDOM_SEED
from src.models.tmc import TrustedMultiViewClassifier, train_view_model, predict_evidence, softplus
from sklearn.preprocessing import StandardScaler

np.random.seed(RANDOM_SEED)

print("="*60)
print("STEP 03: MODEL TRAINING")
print("="*60)

# --- 1. Load Features ---
df = pd.read_csv(TMC_FEATURES)
print(f"Loaded {len(df)} drugs.")

# Chemical Features (SVD-reduced ECFP + PhysChem)
chem_cols = [c for c in df.columns if c.startswith('ecfp_svd_') or c in ['logp', 'mw', 'tpsa', 'hbd', 'hba']]
X_chem = df[chem_cols].values

# Network Features
net_cols = ['I_T', 'ecnp_mu_T', 'ecnp_sigma_T', 'ecnp_z']
X_net = df[net_cols].values

# Quality Gating Signal
quality = df['k_log'].values
y = df['is_dili'].values

# Handle NaN values (impute with column mean)
from numpy import nanmean
X_chem = np.where(np.isnan(X_chem), np.nanmean(X_chem, axis=0), X_chem)
X_net = np.where(np.isnan(X_net), np.nanmean(X_net, axis=0), X_net)
print(f"  Imputed NaNs. Chem shape: {X_chem.shape}, Net shape: {X_net.shape}")

# --- 2. Cross-Validation Training ---
print("Running 5-Fold Stratified Cross-Validation...")
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_SEED)

all_preds = np.zeros(len(y))

for fold, (train_idx, test_idx) in enumerate(cv.split(X_chem, y)):
    clf = TrustedMultiViewClassifier(quality_threshold=np.percentile(quality[train_idx], 25))
    clf.fit(X_chem[train_idx], X_net[train_idx], y[train_idx])
    probs = clf.predict_proba(X_chem[test_idx], X_net[test_idx], quality[test_idx])
    all_preds[test_idx] = probs[:, 1]

# --- 3. Evaluate ---
auc = roc_auc_score(y, all_preds)
pr_auc = average_precision_score(y, all_preds)

print(f"\nTMC Model Results:")
print(f"  ROC-AUC: {auc:.3f}")
print(f"  PR-AUC:  {pr_auc:.3f}")

# --- 4. Save ---
output = df[['dilirank_name', 'smiles', 'is_dili']].copy()
output['predicted_prob'] = all_preds
output.to_csv(RESULTS_DIR / 'final_dili_predictions.csv', index=False)
print(f"Saved: {RESULTS_DIR / 'final_dili_predictions.csv'}")

print("\n[STEP 03 COMPLETE]")

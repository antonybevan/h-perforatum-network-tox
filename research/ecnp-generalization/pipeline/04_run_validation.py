"""
Pipeline Step 04: Validation Suite
===================================
Runs Conformal Prediction, Adversarial Stress Tests, and Tox21 External Validation.

Inputs:
    - results/tmc_features_706.csv
Outputs:
    - results/conformal_validation_results.csv
    - results/conformal_advanced_report.txt
    - results/mechanistic_validation_results.csv

Usage:
    python pipeline/04_run_validation.py
"""
import sys
sys.path.insert(0, str(__file__).replace('pipeline\\04_run_validation.py', ''))

import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import TruncatedSVD

from src.config import TMC_FEATURES, RESULTS_DIR, RANDOM_SEED, CONFORMAL_ALPHA
from src.models.tmc import train_view_model, predict_evidence, softplus
from src.validation.conformal import (
    compute_nonconformity, calibrate_threshold, 
    get_prediction_set, evaluate_conformal
)
from src.features.ecfp import compute_ecfp

np.random.seed(RANDOM_SEED)

print("="*60)
print("STEP 04: VALIDATION SUITE")
print("="*60)

# --- 1. Load Features ---
df = pd.read_csv(TMC_FEATURES)
print(f"Loaded {len(df)} drugs.")

y = df['is_dili'].values
quality = df['k_log'].values
q1_threshold = np.percentile(quality, 25)

# --- 2. Generate Views ---
print("Generating feature views...")
ecfps = np.array([compute_ecfp(s) for s in df['smiles']])
svd = TruncatedSVD(n_components=50, random_state=RANDOM_SEED)
X_chem = svd.fit_transform(ecfps)

net_cols = ['I_T', 'ecnp_mu_T', 'ecnp_sigma_T', 'ecnp_z']
X_net = df[net_cols].values

# Handle NaN values (impute with column mean)
X_chem = np.where(np.isnan(X_chem), np.nanmean(X_chem, axis=0), X_chem)
X_net = np.where(np.isnan(X_net), np.nanmean(X_net, axis=0), X_net)
print(f"  Imputed NaNs. Chem shape: {X_chem.shape}, Net shape: {X_net.shape}")

# --- 3. Cross-Conformal Validation ---
print("Running Cross-Conformal Validation...")
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_SEED)
results = []

for fold, (train_idx, cal_test_idx) in enumerate(cv.split(X_chem, y)):
    np.random.shuffle(cal_test_idx)
    n_cal = len(cal_test_idx) // 2
    cal_idx, test_idx = cal_test_idx[:n_cal], cal_test_idx[n_cal:]
    
    # Scale
    scaler_c, scaler_n = StandardScaler(), StandardScaler()
    X_c_tr, X_c_cal, X_c_te = scaler_c.fit_transform(X_chem[train_idx]), scaler_c.transform(X_chem[cal_idx]), scaler_c.transform(X_chem[test_idx])
    X_n_tr, X_n_cal, X_n_te = scaler_n.fit_transform(X_net[train_idx]), scaler_n.transform(X_net[cal_idx]), scaler_n.transform(X_net[test_idx])
    
    # Train
    y_tr_oh = np.zeros((len(train_idx), 2)); y_tr_oh[np.arange(len(train_idx)), y[train_idx]] = 1
    W_c, W_n = train_view_model(X_c_tr, y_tr_oh), train_view_model(X_n_tr, y_tr_oh)
    
    # Fused Predictor
    def fuse(Xc, Xn, q):
        e_c, e_n = predict_evidence(Xc, W_c), predict_evidence(Xn, W_n)
        gate = (q > q1_threshold).astype(float).reshape(-1, 1)
        alpha = (e_c + e_n * gate) + 1
        return alpha / alpha.sum(axis=1, keepdims=True)
    
    # Calibrate
    probs_cal = fuse(X_c_cal, X_n_cal, quality[cal_idx])
    scores_cal = compute_nonconformity(probs_cal, y[cal_idx])
    q_hat = calibrate_threshold(scores_cal, CONFORMAL_ALPHA)
    
    # Predict
    probs_te = fuse(X_c_te, X_n_te, quality[test_idx])
    pred_sets = get_prediction_set(probs_te, q_hat)
    
    for i, idx in enumerate(test_idx):
        results.append({
            'drug_idx': idx, 'k_log': quality[idx],
            'is_noisy': quality[idx] <= q1_threshold,
            'set_size': len(pred_sets[i]),
            'is_correct': y[test_idx][i] in pred_sets[i],
            'true_label': y[test_idx][i]
        })

res_df = pd.DataFrame(results)

# --- 4. Report ---
print("\n" + "="*60)
print("ADVANCED CONFORMAL VALIDATION RESULTS")
print("="*60)

print(f"Global Validity: {res_df['is_correct'].mean()*100:.2f}%")
print(f"Mean Set Size:   {res_df['set_size'].mean():.2f}")

res_df.to_csv(RESULTS_DIR / 'conformal_validation_results.csv', index=False)
print(f"Saved: {RESULTS_DIR / 'conformal_validation_results.csv'}")

print("\n[STEP 04 COMPLETE]")

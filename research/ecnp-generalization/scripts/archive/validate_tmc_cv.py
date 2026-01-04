import pandas as pd
import numpy as np
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, average_precision_score
from scipy.optimize import minimize
from scipy.special import digamma, gammaln

# --- 1. Math Utils (Identical to train_tmc_scipy.py) ---

def softplus(x):
    return np.log1p(np.exp(np.clip(x, -20, 20)))

def edl_loss_numpy(weights, X, y_onehot, num_classes):
    D = X.shape[1]
    K = num_classes
    W = weights.reshape(D, K)
    logits = X @ W
    evidence = softplus(logits)
    alpha = evidence + 1
    S = np.sum(alpha, axis=1, keepdims=True)
    err = np.sum( (y_onehot - (alpha / S))**2, axis=1, keepdims=True )
    var = np.sum( (alpha * (S - alpha)) / (S*S*(S+1)), axis=1, keepdims=True )
    loss = np.mean(err + var) 
    reg = 0.01 * np.sum(W**2)
    return loss + reg

# --- 2. Helper Functions ---

def generate_ecfp(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
            return list(fp)
        return [0]*1024
    except:
        return [0]*1024

def train_view_model(X, y_onehot):
    D = X.shape[1]
    K = 2
    initial_weights = np.random.randn(D * K) * 0.01
    res = minimize(
        fun=edl_loss_numpy,
        x0=initial_weights,
        args=(X, y_onehot, K),
        method='L-BFGS-B',
        options={'maxiter': 200, 'disp': False} # Suppress per-fold output
    )
    return res.x.reshape(D, K)

def predict_view(X, W):
    logits = X @ W
    return softplus(logits) # Evidence

# --- 3. CV Loop ---

def main():
    print("Loading data for Integrity Validation (5-Fold CV)...")
    input_path = r"v:\new\h-perforatum-network-tox\research\ecnp-generalization\results\tmc_features_706.csv"
    if not os.path.exists(input_path):
        print("Data not found.")
        return
        
    df = pd.read_csv(input_path)
    y = df['is_dili'].values
    
    # --- Feature Prep (Consistent with Training) ---
    print("Generating Features...")
    ecfps = np.array([generate_ecfp(s) for s in df['smiles']])
    svd = TruncatedSVD(n_components=50, random_state=42)
    X_chem_all = svd.fit_transform(ecfps)
    
    net_cols = ['network_coverage', 'n_targets_log', 'k_log', 'I_T', 'ecnp_mu_T', 'ecnp_sigma_T', 'ecnp_z']
    X_net_all = df[net_cols].values
    
    # Store OOF Predictions
    oof_preds = np.zeros(len(df))
    oof_uncertainty = np.zeros(len(df))
    
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    
    print("\nStarting Cross-Validation...")
    
    for fold, (train_idx, test_idx) in enumerate(skf.split(X_chem_all, y)):
        print(f"  Fold {fold+1}/5...")
        
        # Split
        X_c_train, X_c_test = X_chem_all[train_idx], X_chem_all[test_idx]
        X_n_train, X_n_test = X_net_all[train_idx], X_net_all[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]
        
        # Scaling (Must fit on Train, transform Test to prevent leakage)
        scaler_c = StandardScaler()
        X_c_train = scaler_c.fit_transform(X_c_train)
        X_c_test = scaler_c.transform(X_c_test)
        
        scaler_n = StandardScaler()
        X_n_train = scaler_n.fit_transform(X_n_train)
        X_n_test = scaler_n.transform(X_n_test)
        
        # Prep Y
        y_train_oh = np.zeros((len(y_train), 2))
        y_train_oh[np.arange(len(y_train)), y_train] = 1
        
        # Train
        W_chem = train_view_model(X_c_train, y_train_oh)
        W_net = train_view_model(X_n_train, y_train_oh)
        
        # Predict (Get Evidence)
        e_c = predict_view(X_c_test, W_chem)
        e_n = predict_view(X_n_test, W_net)
        
        # Fusion
        e_fused = e_c + e_n
        alpha_fused = e_fused + 1
        alpha_n = e_n + 1
        
        # Compute Prob and Uncertainty
        probs = alpha_fused[:, 1] / np.sum(alpha_fused, axis=1)
        u_net = 2 / np.sum(alpha_n, axis=1)
        
        # Store
        oof_preds[test_idx] = probs
        oof_uncertainty[test_idx] = u_net
        
    # --- Final validated Analysis ---
    df['cv_prob'] = oof_preds
    df['cv_uncertainty'] = oof_uncertainty
    
    print("\n" + "="*50)
    print("INTEGRITY VALIDATION RESULTS (5-Fold CV)")
    print("="*50)
    
    auc_global = roc_auc_score(y, oof_preds)
    ap_global = average_precision_score(y, oof_preds)
    
    print(f"Global Validation AUC:    {auc_global:.4f}")
    print(f"Global Validation PR-AUC: {ap_global:.4f}")
    
    # Stratified Analysis (Using same logic as before)
    # Use k_log for bins
    col = 'k_log' if 'k_log' in df.columns else 'n_targets_mapped'
    try:
        df['target_bin'] = pd.qcut(df[col], q=4, labels=['Q1 (Low)', 'Q2', 'Q3', 'Q4 (High)'])
        print(f"\nStratified Performance by {col}:")
        for label in ['Q1 (Low)', 'Q2', 'Q3', 'Q4 (High)']:
            subset = df[df['target_bin'] == label]
            if len(subset) > 0 and subset['is_dili'].nunique() > 1:
                auc_sub = roc_auc_score(subset['is_dili'], subset['cv_prob'])
                print(f"  {label} Targets: AUC = {auc_sub:.4f}")
    except:
        pass
        
    # Validation of Uncertainty vs Quality (In CV mode)
    print("\nRegime Detection Integrity check (CV Predictions):")
    # Low Coverage should still be uncertain even in CV
    df['coverage_bin'] = pd.cut(df['network_coverage'], bins=[-0.1, 0.3, 0.7, 1.1], labels=['Low', 'Medium', 'High'])
    print(df.groupby('coverage_bin', observed=False)[['cv_uncertainty']].mean())

    # Save
    df.to_csv(r"v:\new\h-perforatum-network-tox\research\ecnp-generalization\results\tmc_cv_validation_results.csv", index=False)

if __name__ == "__main__":
    main()

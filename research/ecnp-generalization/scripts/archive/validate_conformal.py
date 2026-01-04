import pandas as pd
import numpy as np
import os
from sklearn.model_selection import StratifiedKFold
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
from scipy.optimize import minimize
from scipy.special import xlogy
from rdkit import Chem
from rdkit.Chem import AllChem

# --- 1. Math & Model Utils ---

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
    reg = 0.01 * np.sum(W**2)
    return np.mean(err + var) + reg

def train_view_model(X, y_onehot):
    D = X.shape[1]
    K = 2
    initial_weights = np.random.randn(D * K) * 0.01
    res = minimize(
        fun=edl_loss_numpy,
        x0=initial_weights,
        args=(X, y_onehot, K),
        method='L-BFGS-B',
        options={'maxiter': 200, 'disp': False}
    )
    return res.x.reshape(D, K)

def predict_evidence(X, W):
    return softplus(X @ W)

# --- 2. Conformal Prediction Utils ---

def compute_nonconformity(probs, y_true_idx):
    # Standard conformity score for classification: 1 - softmax_prob_of_true_class
    # Or simpler: probability of true class. 
    # Standard ICP: non_conformity = 1 - f(x)_y
    # We want a score where HIGH means "Strange/Wrong".
    return 1.0 - probs[np.arange(len(probs)), y_true_idx]

def get_prediction_set(probs, q_hat):
    # Return set of classes where non_conformity <= q_hat
    # i.e., 1 - prob <= q_hat  =>  prob >= 1 - q_hat
    # So include class C if P(C) >= 1 - q_hat
    threshold = 1.0 - q_hat
    
    sets = []
    for p_vec in probs:
        classes = np.where(p_vec >= threshold)[0]
        sets.append(classes)
    return sets

# --- 3. Main Symbolic Logic ---

def generate_ecfp(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
            return list(fp)
        return [0]*1024
    except:
        return [0]*1024

def main():
    # Set Global Seed for Reproducibility
    np.random.seed(42)
    
    print("Initializing Conformal Integrity Protocol...")
    input_path = r"v:\new\h-perforatum-network-tox\research\ecnp-generalization\results\tmc_features_706.csv"
    if not os.path.exists(input_path):
        print("Data not found.")
        return
        
    df = pd.read_csv(input_path)
    y = df['is_dili'].values
    
    # Stratification Bin (Approximated by Target Quality)
    # Tier 2b (Curated) vs Tier 2a (Noisy)
    # We define "Noisy/Sparse" as k_log quartiles Q1 (Very few targets)
    # And "Noisy/Promiscuous" as Q4 (Too many targets, potentially)
    # Let's use k_log as the proxy for quality.
    col = 'k_log' if 'k_log' in df.columns else 'n_targets_mapped'
    df['quality_score'] = df[col]
    
    # --- Feature Prep ---
    print("Generating Features...")
    ecfps = np.array([generate_ecfp(s) for s in df['smiles']])
    svd = TruncatedSVD(n_components=50, random_state=42)
    X_chem_all = svd.fit_transform(ecfps)
    
    net_cols = ['I_T', 'ecnp_mu_T', 'ecnp_sigma_T', 'ecnp_z'] # Note: removed epistemic features from INPUT, used in GATE
    # We treat the network model as "blind" to quality, so the Symbolic Gate does the work.
    X_net_all = df[net_cols].values
    
    # Hard Gate Parameter
    # Refined Definition: "Noisy" = Low Information Content (Sparse Targets)
    # We define Noisy as the bottom quartile (Q1) of target counts.
    # Tier 2a is characterized by having very few reliable targets.
    
    k_log_values = df[col].values
    q1_threshold = np.percentile(k_log_values, 25)
    
    print(f"Defining Noisy Regime: k_log <= {q1_threshold:.4f} (Bottom 25%)")
    
    # 5-Fold Stratified CV for Conformal Calibration
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    
    results = []
    
    print("\nRunning Cross-Conformal Validation...")
    
    for fold, (train_idx, cal_test_idx) in enumerate(skf.split(X_chem_all, y)):
        # Split Calibration+Test into Cal and Test (50/50)
        np.random.shuffle(cal_test_idx)
        n_cal = len(cal_test_idx) // 2
        cal_idx = cal_test_idx[:n_cal]
        test_idx = cal_test_idx[n_cal:]
        
        # Scaling
        scaler_c = StandardScaler()
        X_c_tr = scaler_c.fit_transform(X_chem_all[train_idx])
        X_c_cal = scaler_c.transform(X_chem_all[cal_idx])
        X_c_te = scaler_c.transform(X_chem_all[test_idx])
        
        scaler_n = StandardScaler()
        X_n_tr = scaler_n.fit_transform(X_net_all[train_idx])
        X_n_cal = scaler_n.transform(X_net_all[cal_idx])
        X_n_te = scaler_n.transform(X_net_all[test_idx])
        
        y_tr = y[train_idx]
        y_tr_oh = np.zeros((len(y_tr), 2))
        y_tr_oh[np.arange(len(y_tr)), y_tr] = 1
        
        # Train Models (Blind to Quality)
        W_chem = train_view_model(X_c_tr, y_tr_oh)
        W_net = train_view_model(X_n_tr, y_tr_oh)
        
        # --- Define Predictor Function (With Symbolic Gate) ---
        def get_fused_probs(Xc, Xn, quality_vector):
            # 1. Get Raw Evidence
            e_c = predict_evidence(Xc, W_chem)
            e_n = predict_evidence(Xn, W_net)
            
            # 2. Apply Symbolic Gate
            # Gate: 1 if quality > Q1_threshold, else 0
            gate = (quality_vector > q1_threshold).astype(float).reshape(-1, 1)
            e_n_gated = e_n * gate
            
            # 3. Fuse
            e_fused = e_c + e_n_gated
            alpha = e_fused + 1
            return alpha / np.sum(alpha, axis=1, keepdims=True)
            
        # --- Conformal Calibration ---
        
        # 1. Get Probs for Calibration Set
        probs_cal = get_fused_probs(X_c_cal, X_n_cal, k_log_values[cal_idx])
        y_cal = y[cal_idx]
        
        # 2. Compute Scores
        scores_cal = compute_nonconformity(probs_cal, y_cal)
        
        # 3. Compute Quantile (1 - alpha)
        alpha_conf = 0.1
        n = len(scores_cal)
        q_val = np.quantile(scores_cal, np.ceil((n+1)*(1-alpha_conf))/n, method='higher')
        
        # --- Prediction on Test Set ---
        probs_test = get_fused_probs(X_c_te, X_n_te, k_log_values[test_idx])
        y_test = y[test_idx]
        
        pred_sets = get_prediction_set(probs_test, q_val)
        
        # --- Record Results ---
        for i, idx in enumerate(test_idx):
            set_size = len(pred_sets[i])
            is_correct = y_test[i] in pred_sets[i]
            qual_val = k_log_values[idx]
            
            results.append({
                'drug_idx': idx,
                'k_log': qual_val,
                'is_noisy': qual_val <= q1_threshold, # Defined as Sparse
                'set_size': set_size,
                'is_correct': is_correct,
                'true_label': y_test[i],
                'q_hat': q_val
            })
            
    # --- Analysis ---
    res_df = pd.DataFrame(results)
    
    print("\n" + "="*60)
    print("ADVANCED CONFORMAL VALIDATION RESULTS (Angelopoulos Standards)")
    print("="*60)
    
    # 1. Global Metrics
    global_validity = res_df['is_correct'].mean()
    mean_set_size = res_df['set_size'].mean()
    print(f"\n[Global] Validity:      {global_validity*100:.2f}%  (Target: 90%)")
    print(f"[Global] Mean Set Size: {mean_set_size:.2f}    (Lower is Better)")

    # 2. Set Size Histogram
    print("\n[Distribution] Set Size Histogram:")
    counts = res_df['set_size'].value_counts().sort_index()
    total = len(res_df)
    for size, count in counts.items():
        pct = (count / total) * 100
        print(f"  Size {size}: {count:4d} ({pct:5.1f}%) | {'#' * int(pct/2)}")

    # 3. Size-Stratified Coverage (SSC)
    print("\n[Metric] Size-Stratified Coverage (SSC):")
    print("  *Checking if validity holds across different uncertainty levels*")
    ssc = res_df.groupby('set_size')['is_correct'].agg(['mean', 'count'])
    ssc.columns = ['Coverage', 'Count']
    print(ssc)

    # 4. Label-Conditional Validity
    print("\n[Metric] Label-Conditional Validity:")
    print("  *Checking if validity holds for both Toxic (1) and Safe (0)*")
    lcv = res_df.groupby('true_label')['is_correct'].agg(['mean', 'count'])
    lcv.index = ['Safe (0)', 'Toxic (1)']
    lcv.columns = ['Coverage', 'Count']
    print(lcv)

    # 5. Regime-Conditional Efficiency (MSS)
    print("\n[Metric] Regime-Conditional Efficiency:")
    print("  *Do we produce tighter sets for High-Quality data?*")
    regime = res_df.groupby('is_noisy')[['set_size', 'is_correct']].agg(['mean', 'count'])
    regime.index = ['Clean (Tier 2b)', 'Noisy (Tier 2a)']
    regime.columns = ['MSS', 'Valid %', 'N', 'Valid N'] # N is duplicated, clean up
    print(regime[['MSS', 'Valid %', 'N']])
    
    # Save
    res_df.to_csv(r"v:\new\h-perforatum-network-tox\research\ecnp-generalization\results\conformal_validation_results.csv", index=False)

if __name__ == "__main__":
    main()

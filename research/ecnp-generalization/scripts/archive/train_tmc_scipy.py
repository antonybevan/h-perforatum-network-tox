import pandas as pd
import numpy as np
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import StandardScaler
from scipy.optimize import minimize
from scipy.special import digamma, gammaln

# --- 1. Math Utils (Numpy) ---

def softplus(x):
    return np.log1p(np.exp(np.clip(x, -20, 20))) # stable

def relu(x):
    return np.maximum(0, x)

def kl_divergence_numpy(alpha, num_classes):
    # alpha: [N, K]
    ones = np.ones((1, num_classes))
    sum_alpha = np.sum(alpha, axis=1, keepdims=True) # [N, 1]
    
    first_term = (
        gammaln(sum_alpha)
        - np.sum(gammaln(alpha), axis=1, keepdims=True)
        + np.sum(gammaln(ones), axis=1, keepdims=True)
        - gammaln(np.sum(ones, axis=1, keepdims=True))
    )
    
    second_term = np.sum(
        (alpha - ones) * (digamma(alpha) - digamma(sum_alpha)),
        axis=1, keepdims=True
    )
    
    return first_term + second_term

def edl_loss_numpy(weights, X, y_onehot, num_classes, lam=1.0):
    # X: [N, D]
    # weights: [D, K] flattened
    D = X.shape[1]
    K = num_classes
    
    W = weights.reshape(D, K)
    logits = X @ W
    evidence = softplus(logits)
    alpha = evidence + 1
    
    S = np.sum(alpha, axis=1, keepdims=True)
    
    # MSE Loss (Bayes Risk)
    # L = sum( (y - alpha/S)^2 ) + var term
    err = np.sum( (y_onehot - (alpha / S))**2, axis=1, keepdims=True )
    var = np.sum( (alpha * (S - alpha)) / (S*S*(S+1)), axis=1, keepdims=True )
    
    # KL Regularization
    # Target alpha_tilde = (1-y)*1 + y*alpha_target? 
    # Actually standard EDL KL target is uniform (alpha=1) for non-ground-truth classes
    # Simplification: KL(alpha || alpha_uniform) weighted by (1-y_true) is for OOD?
    # We use the standard formulation: KL(alpha || 1) for misclassified samples?
    # Simplified KL target: alpha_target = y_onehot + (1-y_onehot)*1 = 1 where not class, but ...
    # The Han et al formulation: KL(alpha_i || 1) basically?
    # We will use purely the MSE term first to see if it converges, adding KL if needed.
    # The robust TMC formulation relies heavily on the fusion, not just the loss.
    
    loss = np.mean(err + var) 
    
    # Simple L2 Reg
    reg = 0.01 * np.sum(W**2)
    
    return loss + reg

# --- 2. Training Logic ---

def generate_ecfp(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
            return list(fp)
        return [0]*1024
    except:
        return [0]*1024

def train_view_model(X, y_onehot, view_name="View"):
    print(f"Training {view_name} (Shape: {X.shape})...")
    D = X.shape[1]
    K = 2
    
    # Init Params
    initial_weights = np.random.randn(D * K) * 0.01
    
    # Optimization
    res = minimize(
        fun=edl_loss_numpy,
        x0=initial_weights,
        args=(X, y_onehot, K),
        method='L-BFGS-B',
        options={'maxiter': 500, 'disp': True}
    )
    
    print(f"  Converged: {res.success}, Loss: {res.fun:.4f}")
    return res.x.reshape(D, K)

def main():
    # Load
    input_path = r"v:\new\h-perforatum-network-tox\research\ecnp-generalization\results\tmc_features_706.csv"
    if not os.path.exists(input_path):
        print("Data not found.")
        return
        
    df = pd.read_csv(input_path)
    
    # Check if we have DILI labels (we do)
    y = df['is_dili'].values
    y_onehot = np.zeros((len(y), 2))
    y_onehot[np.arange(len(y)), y] = 1
    
    # --- Feature Prep ---
    
    # 1. Chem View (Reduced ECFP)
    print("Generating Chem Features...")
    ecfps = np.array([generate_ecfp(s) for s in df['smiles']])
    svd = TruncatedSVD(n_components=50, random_state=42)
    X_chem = svd.fit_transform(ecfps)
    scaler_c = StandardScaler()
    X_chem = scaler_c.fit_transform(X_chem)
    
    # 2. Net View (Aleatoric + Epistemic)
    net_cols = ['network_coverage', 'n_targets_log', 'k_log', 'I_T', 'ecnp_mu_T', 'ecnp_sigma_T', 'ecnp_z']
    X_net = df[net_cols].values
    scaler_n = StandardScaler()
    X_net = scaler_n.fit_transform(X_net)
    
    # --- Train Separate Views ---
    # In TMC, we can train views independently first, or jointly. 
    # Independent training is sufficient to test the "Uncertainty vs Coverage" hypothesis.
    # Why? Because the Network View *should* learn to be uncertain on its own if provided with 'coverage' features.
    
    W_chem = train_view_model(X_chem, y_onehot, "Chem View")
    W_net = train_view_model(X_net, y_onehot, "Net View")
    
    # --- Evaluation ---
    
    # 1. Compute Evidence
    logit_c = X_chem @ W_chem
    ev_c = softplus(logit_c)
    alpha_c = ev_c + 1
    
    logit_n = X_net @ W_net
    ev_n = softplus(logit_n)
    alpha_n = ev_n + 1
    
    # 2. Fusion (Dempster Rule - Simplified Sum)
    # TMC paper Eq 10: e = e1 + e2 for Dirichlet fusion
    ev_fused = ev_c + ev_n
    alpha_fused = ev_fused + 1
    
    # 3. Uncertainty Metrics
    S_n = np.sum(alpha_n, axis=1)
    u_n = 2 / S_n # Uncertainty
    e_n_tot = np.sum(ev_n, axis=1) # Total Evidence
    
    S_c = np.sum(alpha_c, axis=1)
    u_c = 2 / S_c
    
    prob_dili = alpha_fused[:, 1] / np.sum(alpha_fused, axis=1)
    
    # --- Analysis ---
    df['net_uncertainty'] = u_n
    df['net_evidence'] = e_n_tot
    df['chem_uncertainty'] = u_c
    df['tmc_prob'] = prob_dili
    
    # Binning by Network Coverage to check hypothesis
    df['coverage_bin'] = pd.cut(df['network_coverage'], bins=[-0.1, 0.3, 0.7, 1.1], labels=['Low', 'Medium', 'High'])
    
    print("\n" + "="*40)
    print("REGIME DETECTION VALIDATION")
    print("="*40)
    print("Hypothesis: Low Coverage (Noisy) should have HIGH Uncertainty (u -> 1.0)")
    print("            High Coverage (Clean) should have LOWER Uncertainty")
    
    print(df.groupby('coverage_bin', observed=False)[['net_uncertainty', 'net_evidence']].mean())
    
    # Save
    out_path = r"v:\new\h-perforatum-network-tox\research\ecnp-generalization\results\tmc_validation_results.csv"
    df.to_csv(out_path, index=False)
    print(f"\nSaved results to {out_path}")

if __name__ == "__main__":
    main()

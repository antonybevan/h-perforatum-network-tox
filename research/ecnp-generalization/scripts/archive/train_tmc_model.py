import pandas as pd
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, brier_score_loss, precision_score, recall_score
import matplotlib.pyplot as plt
import seaborn as sns
import os

# --- 1. Evidential Deep Learning Utils (Han et al. 2021) ---

def relu_evidence(y):
    return torch.nn.functional.relu(y)

def exp_evidence(y):
    return torch.exp(torch.clamp(y, -10, 10))

def softplus_evidence(y):
    return torch.nn.functional.softplus(y)

def kl_divergence(alpha, num_classes, device=None):
    ones = torch.ones([1, num_classes], dtype=torch.float32, device=device)
    sum_alpha = torch.sum(alpha, dim=1, keepdim=True)
    first_term = (
        torch.lgamma(sum_alpha)
        - torch.lgamma(alpha).sum(dim=1, keepdim=True)
        + torch.lgamma(ones).sum(dim=1, keepdim=True)
        - torch.lgamma(ones.sum(dim=1, keepdim=True))
    )
    second_term = (
        (alpha - ones)
        .mul(torch.digamma(alpha) - torch.digamma(sum_alpha))
        .sum(dim=1, keepdim=True)
    )
    return first_term + second_term

def loglikelihood_loss(y, alpha, device=None):
    y = y.to(device)
    alpha = alpha.to(device)
    S = torch.sum(alpha, dim=1, keepdim=True)
    loglikelihood_err = torch.sum((y - (alpha / S)) ** 2, dim=1, keepdim=True)
    loglikelihood_var = torch.sum(
        (alpha * (S - alpha) / (S * S * (S + 1))), dim=1, keepdim=True
    )
    return loglikelihood_err + loglikelihood_var

def edl_mse_loss(output, target, epoch_num, num_classes, annealing_step, device=None):
    evidence = relu_evidence(output)
    alpha = evidence + 1
    loss = torch.mean(
        loglikelihood_loss(target, alpha, device)
    )
    
    # KL Regularization: Forces high uncertainty (flat alpha) unless data says otherwise
    kl_alpha = (alpha - 1) * (1 - target) + 1
    kl_div = kl_divergence(kl_alpha, num_classes, device=device)
    annealing_coef = min(1, epoch_num / annealing_step)
    
    total_loss = loss + annealing_coef * torch.mean(kl_div)
    return total_loss

# --- 2. Dempster-Shafer Fusion Layer ---

class DS_Fusion(nn.Module):
    def __init__(self, num_classes):
        super(DS_Fusion, self).__init__()
        self.num_classes = num_classes

    def forward(self, views_evidence):
        # views_evidence: List of [batch_size, num_classes] tensors
        # Each view outputs 'evidence' (relu output) > 0
        
        # Calculate belief and uncertainty for each view
        # alpha = evidence + 1
        # S = sum(alpha)
        # b = evidence / S
        # u = num_classes / S
        
        dim = 1
        # Initialize fused belief/uncertainty
        # We start fusion logic based on TMC paper Eq 6 & 7 (simplified for 2 views)
        
        # For implementation stability, we compute in mass space
        count = 0
        for v_evidence in views_evidence:
            alpha = v_evidence + 1
            S = torch.sum(alpha, dim=dim, keepdim=True)
            b = v_evidence / S
            u = self.num_classes / S
            
            if count == 0:
                M_fused = b # Belief Mass
                u_fused = u # Uncertainty Mass
            else:
                # Dempster's Rule of Combination (Binary case simplified)
                # M = (M1*u2 + M2*u1 + M1*M2) / (1 - C) ... this is complex for multi-class
                # Han et al. 2021 propose a closed form for Dirichlet fusion:
                # alpha_fused = alpha1 + alpha2 - 1 is WRONG (that's naive summing)
                # TMC Paper: "The fused evidence e = e1 + e2"
                # Wait, TMC explicitly simplifies this! 
                # If we assume independent evidence sources, e_total = e1 + e2. 
                # This corresponds to alpha_total = alpha1 + alpha2 - 1.
                # Let's verify: alpha = e + 1. 
                # alpha_new = (e1 + e2) + 1 = (alpha1-1) + (alpha2-1) + 1 = alpha1 + alpha2 - 1.
                # Yes, for Dirichlet paramertization under DST, evidences sum up.
                pass
            count += 1
            
        # The TMC paper "Trusted Multi-View Classification" (Eq 10) confirms:
        # "To combine the evidence from different views... we can simply add the evidence..."
        # e_fused = sum(e_v)
        
        e_fused = torch.zeros_like(views_evidence[0])
        for v in views_evidence:
            e_fused += v
            
        return e_fused

# --- 3. The Model Architecture ---

class TrustedMultiViewClassifier(nn.Module):
    def __init__(self, chem_input_dim, net_input_dim, num_classes=2):
        super(TrustedMultiViewClassifier, self).__init__()
        
        # View 1: Chemistry (Aleatoric Only - we assume Chem data is always "reliable" signals)
        self.chem_view = nn.Sequential(
            nn.Linear(chem_input_dim, 128),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(128, 64),
            nn.ReLU(),
            nn.Linear(64, num_classes) # Output: Logits (will be relu'd for evidence)
        )
        
        # View 2: Network (Regime-Aware)
        # Input includes Aleatoric (I_T, mu, sigma) AND Epistemic (Coverage, N_targets)
        # The network should learn to output Evidence=0 if Epistemic features are bad
        self.net_view = nn.Sequential(
            nn.Linear(net_input_dim, 64),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(64, 32),
            nn.ReLU(),
            nn.Linear(32, num_classes) # Output: Logits
        )
        
        self.ds_fusion = DS_Fusion(num_classes)

    def forward(self, x_chem, x_net):
        # Step 1: Compute Evidence for each view
        logits_chem = self.chem_view(x_chem)
        logits_net = self.net_view(x_net)
        
        evidence_chem = relu_evidence(logits_chem)
        evidence_net = relu_evidence(logits_net)
        
        # Step 2: Fuse Evidence
        evidence_fused = self.ds_fusion([evidence_chem, evidence_net])
        
        # Determine Alpha for loss
        alpha_chem = evidence_chem + 1
        alpha_net = evidence_net + 1
        alpha_fused = evidence_fused + 1
        
        return alpha_fused, alpha_chem, alpha_net

# --- 4. Training Script ---

def generate_ecfp(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
            return list(fp)
        else:
            return [0]*1024
    except:
        return [0]*1024

def train_tmc():
    # Load Data
    df = pd.read_csv(r"v:\new\h-perforatum-network-tox\research\ecnp-generalization\results\tmc_features_706.csv")
    
    # 1. Prepare Chemistry Features (ECFP4)
    print("Generating ECFP4 fingerprints...")
    ecfp_list = [generate_ecfp(s) for s in df['smiles']]
    X_chem = np.array(ecfp_list)
    
    # 2. Prepare Network Features (Aleatoric + Epistemic)
    # Our 'tmc_features_706.csv' has columns: 'network_coverage', 'n_targets_log', 'k_log', 'I_T', 'ecnp_mu_T', 'ecnp_sigma_T', 'ecnp_z'
    net_cols = ['network_coverage', 'n_targets_log', 'k_log', 'I_T', 'ecnp_mu_T', 'ecnp_sigma_T', 'ecnp_z']
    X_net = df[net_cols].values
    
    # Normalize Network Features (Crucial for Neural Nets)
    scaler = StandardScaler()
    X_net = scaler.fit_transform(X_net)
    
    # 3. Targets
    y = df['is_dili'].values
    
    # One-hot encode targets for EDL loss
    y_onehot = np.zeros((y.size, 2))
    y_onehot[np.arange(y.size), y] = 1
    
    # Convert to Tensors
    X_chem_t = torch.FloatTensor(X_chem)
    X_net_t = torch.FloatTensor(X_net)
    y_t = torch.FloatTensor(y_onehot)
    
    # 4. Train/Test Split (Simple Random for Arch Validation, will do Stratified later)
    # Actually, for "Regime Validation", we want to see if the model treats Tier 2a differently.
    # We will use the 'n_targets' column later to analyze results.
    # Let's train on 80% of data.
    
    dataset = TensorDataset(X_chem_t, X_net_t, y_t)
    train_size = int(0.8 * len(dataset))
    test_size = len(dataset) - train_size
    train_dataset, test_dataset = torch.utils.data.random_split(dataset, [train_size, test_size])
    
    train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)
    test_loader = DataLoader(test_dataset, batch_size=32, shuffle=False)
    
    # 5. Model Setup
    model = TrustedMultiViewClassifier(chem_input_dim=1024, net_input_dim=len(net_cols), num_classes=2)
    optimizer = optim.Adam(model.parameters(), lr=0.001)
    
    epochs = 50
    annealing_epoch = 10 # KL annealing
    
    print("\nStarting TMC Training...")
    
    for epoch in range(epochs):
        model.train()
        total_loss = 0
        for xc, xn, target in train_loader:
            optimizer.zero_grad()
            alpha_fused, alpha_c, alpha_n = model(xc, xn)
            
            # Loss is sum of losses from all views (to force each view to learn) + fused loss
            loss_c = edl_mse_loss(alpha_c - 1, target, epoch, 2, annealing_epoch)
            loss_n = edl_mse_loss(alpha_n - 1, target, epoch, 2, annealing_epoch)
            loss_f = edl_mse_loss(alpha_fused - 1, target, epoch, 2, annealing_epoch)
            
            loss = loss_c + loss_n + loss_f
            loss.backward()
            optimizer.step()
            total_loss += loss.item()
            
        if epoch % 10 == 0:
            print(f"Epoch {epoch}: Loss = {total_loss / len(train_loader):.4f}")
            
    # 6. Evaluation & Investigation
    print("\nEvaluating Model & Extracting Uncertainty...")
    model.eval()
    results = []
    
    # We run on the FULL dataset to analyze "Regime Detection" properties across all drugs
    with torch.no_grad():
        alpha_f, alpha_c, alpha_n = model(X_chem_t, X_net_t)
        
        # Calculate Evidence & Uncertainty
        # u = K / S, where S = sum(alpha)
        # Evidence = alpha - 1
        
        # Network View Metrics
        S_n = torch.sum(alpha_n, dim=1)
        u_n = 2 / S_n # Uncertainty of Network View
        e_n = torch.sum(alpha_n - 1, dim=1) # Total Evidence of Network View
        
        # Chem View Metrics
        S_c = torch.sum(alpha_c, dim=1)
        u_n_chem = 2 / S_c
        
        # Predictions (Expected Probability)
        prob_dili = alpha_f[:, 1] / torch.sum(alpha_f, dim=1)
        
    # Create Analysis DataFrame
    analysis_df = df.copy()
    analysis_df['tmc_prob_dili'] = prob_dili.numpy()
    analysis_df['net_uncertainty'] = u_n.numpy()
    analysis_df['net_evidence'] = e_n.numpy()
    analysis_df['chem_uncertainty'] = u_n_chem.numpy()
    
    # Save Results
    analysis_df.to_csv(r"v:\new\h-perforatum-network-tox\research\ecnp-generalization\results\tmc_validation_results.csv", index=False)
    
    # --- Check: Does Network Uncertainty correlate with Data Quality? ---
    # Low Coverage -> Expect High Uncertainty (u_n approx 1.0)
    # High Coverage -> Expect Low Uncertainty (u_n < 0.5)
    
    print("\n--- Regime Detection Check ---")
    # Bin by coverage
    analysis_df['coverage_bin'] = pd.cut(analysis_df['network_coverage'], bins=[-0.1, 0.2, 0.8, 1.1], labels=['Low', 'Medium', 'High'])
    print(analysis_df.groupby('coverage_bin')[['net_uncertainty', 'net_evidence']].mean())
    
    print("\nResults saved to: research/ecnp-generalization/results/tmc_validation_results.csv")

if __name__ == "__main__":
    train_tmc()

"""
Final Publication DILI Model
============================

This script trains and evaluates the complete hierarchy of DILI prediction models for the manuscript:

1. TIER 1: BROAD SCREENING (901 drugs)
   - Features: Chemistry Only (ECFP4 + PhysChem)
   - Goal: Baseline performance on maximum available data.

2. TIER 2a: GENERALIZATION VALIDATION (456 drugs)
   - Features: Chemistry + ECNP (Automated Mapping)
   - Goal: Prove ECNP signal persists in broader, noisier dataset.

3. TIER 2b: HIGH-CONFIDENCE MECHANISTIC MODEL (202 drugs)
   - Features: Chemistry + ECNP (Curated Targets)
   - Goal: The "Gold Standard" model for publication (High AUC).

Outputs:
- final_model_results.csv: Consolidated metrics.
- final_summary.txt: Formatted results for the manuscript.
"""
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
import warnings
warnings.filterwarnings('ignore')

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, Lipinski
    RDKIT_AVAILABLE = True
except:
    RDKIT_AVAILABLE = False
    print("WARNING: RDKit not found. Will use cached ECFP if available or fail.")

ROOT = Path(r'v:\new\h-perforatum-network-tox')

print("="*70)
print("FINAL PUBLICATION DILI MODEL TRAINING")
print("="*70)

# =============================================================================
# DATA LOADING
# =============================================================================
print("\n--- Loading Datasets ---")

# 1. Tier 1 (Full DILIrank)
df_tier1 = pd.read_csv(ROOT / 'research/ecnp-generalization/results/dilirank_full_smiles.csv')
print(f"Tier 1 (Full): {len(df_tier1)} drugs")

# 2. Tier 2a (Expanded ECNP)
df_tier2a = pd.read_csv(ROOT / 'research/ecnp-generalization/results/dilirank_706_with_ecnp.csv')
# Filter for valid ECNP
df_tier2a = df_tier2a[df_tier2a['ecnp_z'].notna()].copy()
print(f"Tier 2a (Expanded): {len(df_tier2a)} drugs with valid ECNP")

# 3. Tier 2b (Curated 202)
df_tier2b = pd.read_csv(ROOT / 'research/ecnp-generalization/results/ecfp_model_results.csv')
print(f"Tier 2b (Curated): {len(df_tier2b)} drugs (Gold Standard)")

# =============================================================================
# FEATURE GENERATION AND EVALUATION FUNCTIONS
# =============================================================================

def get_chem_features(smiles_list):
    """Compute ECFP4 and PhysChem features"""
    ecfp = []
    physchem = []
    
    if not RDKIT_AVAILABLE:
        return None, None

    print(f"Computing features for {len(smiles_list)} compounds...")
    for i, s in enumerate(smiles_list):
        try:
            mol = Chem.MolFromSmiles(s)
            if mol:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
                ecfp.append(list(fp))
                physchem.append({
                    'logp': Descriptors.MolLogP(mol),
                    'mw': Descriptors.MolWt(mol),
                    'tpsa': Descriptors.TPSA(mol),
                    'hbd': Lipinski.NumHDonors(mol),
                    'hba': Lipinski.NumHAcceptors(mol)
                })
            else:
                ecfp.append([np.nan]*1024)
                physchem.append({k:np.nan for k in ['logp','mw','tpsa','hbd','hba']})
        except:
            ecfp.append([np.nan]*1024)
            physchem.append({k:np.nan for k in ['logp','mw','tpsa','hbd','hba']})
            
    return pd.DataFrame(ecfp, columns=[f'ecfp_{i}' for i in range(1024)]), pd.DataFrame(physchem)

def evaluate_model_tier1(X, y, name="Model"):
    """Train and evaluate logistic regression with 5-fold CV for Tier 1"""
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    model = LogisticRegression(class_weight='balanced', max_iter=1000, random_state=42)
    
    # Handle NaN
    if np.isnan(X).any():
        print(f"Warning: NaNs in {name} features. Filling with 0.")
        X = np.nan_to_num(X)
        
    try:
        y_pred = cross_val_predict(model, X, y, cv=cv, method='predict_proba')[:, 1]
        auc = roc_auc_score(y, y_pred)
        pr_auc = average_precision_score(y, y_pred)
        return auc, pr_auc
    except Exception as e:
        print(f"Error evaluating {name}: {e}")
        return 0.0, 0.0

def evaluate_model_tier2(df, name="Model"):
    """
    Train and evaluate logistic regression with 5-fold CV for Tier 2 models,
    using decomposed ECNP features and StandardScaler.
    """
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    
    # Feature Engineering
    print("  Generating features...")
    
    # 1. Chemistry Features (Base)
    chem_features = []
    for s in df['smiles']:
        mol = Chem.MolFromSmiles(s)
        if mol:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
            logp = Descriptors.MolLogP(mol)
            mw = Descriptors.MolWt(mol)
            tpsa = Descriptors.TPSA(mol)
            chem_features.append(list(fp) + [logp, mw, tpsa])
        else:
            chem_features.append([0]*1027) # 1024 ECFP + 3 PhysChem
    
    X_chem = np.array(chem_features)
    
    # 2. Decomposed ECNP Features
    # Recalculate network_fraction if missing
    if 'network_fraction' not in df.columns:
        # For Tier 2a (706 set), k is intersection, n_targets_uniprot is total
        # Avoid division by zero
        sizes = df['n_targets_uniprot'].replace(0, 1)
        df['network_fraction'] = df['k'] / sizes
    
    # Select Decomposed Features
    # I_T, mu_T, sigma_T, k, pool_size, network_fraction
    
    # Handle column name mismatches between datasets
    X_ecnp = []
    for idx, row in df.iterrows():
        it = row.get('I_T', 0)
        mu = row.get('ecnp_mu_T', row.get('mu_T', 0))
        k = row.get('k', 0)
        pool = row.get('pool_size', 0)
        frac = row.get('network_fraction', 0)
        
        # Derive sigma_T if feasible
        sigma = row.get('ecnp_sigma_T', row.get('sigma_T', 0))
        if pd.isna(sigma) or sigma == 0:
            # Try to derive from Z definition: Z = (I - mu) / sigma => sigma = (I - mu) / Z
            z = row.get('ecnp_z', row.get('Z', 0))
            if not pd.isna(z) and z != 0:
                sigma = (it - mu) / z
            else:
                sigma = 0
                
        X_ecnp.append([it, mu, sigma, k, pool, frac])
            
    X_ecnp = np.nan_to_num(np.array(X_ecnp))
    
    # Combine
    X_combined = np.hstack([X_chem, X_ecnp])
    y = df['is_dili'].values
    
    # Evaluation with Scaling for LR
    # Use Pipeline for proper scaling inside CV
    model_chem = make_pipeline(StandardScaler(), LogisticRegression(class_weight='balanced', max_iter=1000, random_state=42))
    model_full = make_pipeline(StandardScaler(), LogisticRegression(class_weight='balanced', max_iter=1000, random_state=42))
    
    # 1. Baseline (Chem Only)
    y_pred_base = cross_val_predict(model_chem, X_chem, y, cv=cv, method='predict_proba')[:, 1]
    auc_base = roc_auc_score(y, y_pred_base)
    pr_auc_base = average_precision_score(y, y_pred_base)
    
    # 2. Full (Chem + Decomposed ECNP)
    y_pred_full = cross_val_predict(model_full, X_combined, y, cv=cv, method='predict_proba')[:, 1]
    auc_full = roc_auc_score(y, y_pred_full)
    pr_auc_full = average_precision_score(y, y_pred_full)
    
    lift = auc_full - auc_base
    
    print(f"  [Chemistry Only] AUC: {auc_base:.3f}")
    print(f"  [+ Decomposed]   AUC: {auc_full:.3f}")
    print(f"  [Lift]           {lift:+.3f}")
    
    return auc_base, pr_auc_base, auc_full, pr_auc_full, lift

# =============================================================================
# MODEL TRAINING
# =============================================================================

results = []

# --- MODEL 1: TIER 1 (CHEMISTRY-ONLY BASELINE) ---
print("\n" + "-"*50)
print("1. Training Tier 1 (901 drugs)...")
df_t1_chem, df_t1_phys = get_chem_features(df_tier1['smiles'])
df_t1 = pd.concat([df_tier1, df_t1_chem, df_t1_phys], axis=1).dropna()

X_t1 = df_t1[[c for c in df_t1.columns if 'ecfp_' in c or c in ['logp','mw','tpsa','hbd','hba']]].values
y_t1 = df_t1['is_dili'].values

auc_t1, pr_t1 = evaluate_model_tier1(X_t1, y_t1, "Tier 1")
results.append({'Model': 'Tier 1 (All)', 'N': len(y_t1), 'Features': 'Chemistry', 'AUC': auc_t1, 'PR-AUC': pr_t1, 'Lift': 0.0})
print(f"  AUC: {auc_t1:.3f}")

# --- MODEL 2: TIER 2a (EXPANDED) ---
print("\n" + "-"*50)
print("2. Training Tier 2a (456 Expanded) - Decomposed ECNP...")
# Use the new evaluate_model_tier2 which handles feature engineering internally
auc_base_2a, pr_base_2a, auc_full_2a, pr_full_2a, lift_2a = evaluate_model_tier2(df_tier2a, "Tier 2a")

results.append({'Model': 'Tier 2a (Expanded)', 'N': len(df_tier2a), 'Features': 'Chemistry', 'AUC': auc_base_2a, 'PR-AUC': pr_base_2a, 'Lift': 0.0})
results.append({'Model': 'Tier 2a (Expanded)', 'N': len(df_tier2a), 'Features': 'Chem + Decomposed ECNP', 'AUC': auc_full_2a, 'PR-AUC': pr_full_2a, 'Lift': lift_2a})

# --- MODEL 3: TIER 2b (CURATED GOLD STANDARD) ---
print("\n" + "-"*50)
print("3. Training Tier 2b (202 Curated) - Decomposed ECNP...")
# Ensure ECNP column exists for old checks, but Tier 2 function handles it.
# We just pass the dataframe.
auc_base_2b, pr_base_2b, auc_full_2b, pr_full_2b, lift_2b = evaluate_model_tier2(df_tier2b, "Tier 2b")

results.append({'Model': 'Tier 2b (Curated)', 'N': len(df_tier2b), 'Features': 'Chemistry', 'AUC': auc_base_2b, 'PR-AUC': pr_base_2b, 'Lift': 0.0})
results.append({'Model': 'Tier 2b (Curated)', 'N': len(df_tier2b), 'Features': 'Chem + Decomposed ECNP', 'AUC': auc_full_2b, 'PR-AUC': pr_full_2b, 'Lift': lift_2b})

# =============================================================================
# SAVE RESULTS
# =============================================================================
df_res = pd.DataFrame(results)
output_csv = ROOT / 'research/ecnp-generalization/results/final_model_results.csv'
df_res.to_csv(output_csv, index=False)

print("\n" + "="*70)
print("FINAL PUBLICATION RESULTS (Decomposed ECNP)")
print("="*70)
print(df_res.to_string(index=False))
print("\nSaved to:", output_csv)
print("Done.")

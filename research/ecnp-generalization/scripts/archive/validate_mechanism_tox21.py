
import pandas as pd
import numpy as np
from rdkit import Chem
from scipy.stats import pointbiserialr
import os

# --- PATHS ---
MECH_PATH = r"v:\new\h-perforatum-network-tox\research\ecnp-generalization\results\mechanism_features.csv"
TOX21_PATH = r"v:\new\h-perforatum-network-tox\research\ecnp-generalization\data\external\tox21\tox21.csv.gz"
DILI_SMILES_PATH = r"v:\new\h-perforatum-network-tox\research\ecnp-generalization\results\dilirank_full_smiles.csv" # Need SMILES for DILIrank drugs

OUTPUT_PATH = r"v:\new\h-perforatum-network-tox\research\ecnp-generalization\results\mechanistic_validation_results.csv"

def canonicalize(smiles):
    if pd.isna(smiles): return None
    try:
        mol = Chem.MolFromSmiles(str(smiles))
        if mol:
            return Chem.MolToSmiles(mol, isomericSmiles=False) # Loose matching
    except:
        pass
    return None

def main():
    print("==================================================")
    print("PHASE 12: QUANTITATIVE MECHANISTIC VALIDATION")
    print("==================================================")

    # 1. Load Mechanism Features (Z-scores)
    print(f"Loading Mechanism Features: {MECH_PATH}")
    df_mech = pd.read_csv(MECH_PATH)
    # df_mech has 'dilirank_name'
    
    # 2. Load DILIrank SMILES to bridge the gap
    print(f"Loading DILIrank SMILES: {DILI_SMILES_PATH}")
    if not os.path.exists(DILI_SMILES_PATH):
        print("Error: DILIrank SMILES file not found. Cannot link keys.")
        return
    df_smiles = pd.read_csv(DILI_SMILES_PATH)
    # Assume cols 'dilirank_name', 'SMILES' or similar
    
    # Merge Z-scores with SMILES
    df_merged = pd.merge(df_mech, df_smiles, on='dilirank_name', how='inner')
    print(f"Merged DILIrank drugs with SMILES: {len(df_merged)}")
    
    # Canonicalize DILI SMILES
    print("Canonicalizing DILIrank SMILES...")
    df_merged['canon_smiles'] = df_merged['smiles'].apply(canonicalize)
    df_merged = df_merged.dropna(subset=['canon_smiles'])
    
    # 3. Load Tox21 Data
    print(f"Loading Tox21 Data: {TOX21_PATH}")
    df_tox21 = pd.read_csv(TOX21_PATH, compression='gzip')
    print(f"Tox21 Rows: {len(df_tox21)}")
    
    # Canonicalize Tox21 SMILES
    print("Canonicalizing Tox21 SMILES (this may take a moment)...")
    df_tox21['canon_smiles'] = df_tox21['smiles'].apply(canonicalize)
    df_tox21 = df_tox21.dropna(subset=['canon_smiles'])
    
    # 4. Match Datasets
    print("Matching datasets on SMILES...")
    # Using merge on canon_smiles
    df_final = pd.merge(df_merged, df_tox21, on='canon_smiles', how='inner')
    print(f"Overlap Size: {len(df_final)} compounds found in both DILIrank and Tox21")
    
    if len(df_final) < 10:
        print("Warning: Overlap is too small for significant validation.")
    
    # 5. Correlation Analysis
    # Z_mito vs SR-MMP
    # Z_reactive vs SR-ARE
    # Z_reactive vs SR-p53
    # Z_bile (No direct Tox21 assay, maybe NR-FXR equivalent? No, NR-AhR/PPAR are different)
    
    assays = [
        ('Z_mito', 'SR-MMP'),
        ('Z_reactive', 'SR-ARE'),
        ('Z_reactive', 'SR-p53'),
        ('Z_mito', 'SR-p53') # Mitochondrial tox often leads to apoptosis/p53
    ]
    
    results = []
    
    print("\n--- Correlation Results ---")
    for z_col, assay_col in assays:
        if assay_col not in df_final.columns:
            print(f"Skipping {assay_col} (Not in Tox21 keys)")
            continue
            
        # Drop NaNs for this specific pair
        sub = df_final[[z_col, assay_col]].dropna()
        if len(sub) < 10:
            print(f"Skipping {z_col} vs {assay_col}: Too few valid points ({len(sub)})")
            continue
            
        # Point Biserial: Continuous vs Binary
        corr, pval = pointbiserialr(sub[assay_col], sub[z_col]) # x=binary, y=continuous
        
        # Also compute AUC (ROC) - Can Z-score predict the assay outcome?
        from sklearn.metrics import roc_auc_score
        try:
            auc = roc_auc_score(sub[assay_col], sub[z_col])
        except:
            auc = 0.5
            
        print(f"{z_col} vs {assay_col}: r={corr:.3f} (p={pval:.4f}), AUC={auc:.3f} (n={len(sub)})")
        
        results.append({
            'Hypothesis': f"{z_col} -> {assay_col}",
            'Correlation': corr,
            'P_Value': pval,
            'AUC': auc,
            'N': len(sub)
        })
        
    # Save results
    pd.DataFrame(results).to_csv(OUTPUT_PATH, index=False)
    print(f"\nSaved validation results to {OUTPUT_PATH}")

if __name__ == "__main__":
    main()


import pandas as pd
import numpy as np
import os

def compute_features():
    # Load data
    input_path = r"v:\new\h-perforatum-network-tox\research\ecnp-generalization\results\dilirank_706_with_ecnp.csv"
    output_path = r"v:\new\h-perforatum-network-tox\research\ecnp-generalization\results\tmc_features_706.csv"
    
    if not os.path.exists(input_path):
        print(f"Error: File not found at {input_path}")
        return

    df = pd.read_csv(input_path)
    
    print(f"Loaded {len(df)} compounds.")

    # --- Epistemic Features (Data Quality / Reliability) ---
    # 1. Network Coverage: What fraction of the drug's known biology is in our network?
    # High coverage -> High Trust. Low coverage -> We are missing info -> High Uncertainty.
    # Handle division by zero if n_targets_uniprot is 0
    df['network_coverage'] = df.apply(
        lambda row: row['n_targets_mapped'] / row['n_targets_uniprot'] if row['n_targets_uniprot'] > 0 else 0.0, 
        axis=1
    )
    
    # 2. Target Count (Log): Very low counts (1-2) are often less reliable than moderate counts (5-20).
    # Extremely high counts (>100) might be "dirty" scraped data.
    # We use log1p to compress the scale.
    df['n_targets_log'] = np.log1p(df['n_targets_uniprot'])
    
    # 3. Mapped Count (Log): Raw count of nodes in the network.
    df['k_log'] = np.log1p(df['n_targets_mapped'])

    # --- Feature Sets ---
    
    # Set 1: Epistemic (Gate/Trust)
    epistemic_cols = ['network_coverage', 'n_targets_log', 'k_log']
    
    # Set 2: Aleatoric (Network Signal)
    # Note: 'network_fraction' is effectively k / NetworkSize. 
    # We include raw ECNP components.
    network_cols = ['I_T', 'ecnp_mu_T', 'ecnp_sigma_T', 'ecnp_z']
    
    # Fill missing values (some drugs calculated ECNP might have NaN if k=0)
    for col in network_cols:
        df[col] = df[col].fillna(0.0)
        
    # Standardize Epistemic Features? 
    # Neural Networks like standardized inputs. We will handle normalization in the training script
    # to avoid data leakage (fit scaler on train, transform test).
    # But ensuring no NaNs here is important.
    for col in epistemic_cols:
        df[col] = df[col].fillna(0.0)

    # Save
    cols_to_save = ['dilirank_name', 'drugbank_id', 'is_dili'] + epistemic_cols + network_cols
    if 'smiles' in df.columns:
        cols_to_save.append('smiles')
        
    output_df = df[cols_to_save]
    output_df.to_csv(output_path, index=False)
    
    print("-" * 30)
    print("Feature Engineering Complete")
    print(f"Saved to: {output_path}")
    print("-" * 30)
    print("Preview:")
    print(output_df[['dilirank_name', 'network_coverage', 'n_targets_log']].head())

if __name__ == "__main__":
    compute_features()


import pandas as pd
import numpy as np
import os

def main():
    print("==================================================")
    print("PHASE 11: MECHANISM DECLARATION (LOGIC LAYER)")
    print("==================================================")
    
    # Load Features
    mech_path = r"v:\new\h-perforatum-network-tox\research\ecnp-generalization\results\mechanism_features.csv"
    if not os.path.exists(mech_path):
        print("Mechanism Features not found.")
        return
        
    df = pd.read_csv(mech_path)
    
    # We also need the Probability of DILI to filter (only explain Toxic drugs)
    # Load main results
    res_path = r"v:\new\h-perforatum-network-tox\research\ecnp-generalization\results\tmc_validation_results.csv"
    # Actually, main results might not align row-by-row if indices changed.
    # But mechanism_features has 'dilirank_name'. We should merge.
    
    if os.path.exists(res_path):
        df_res = pd.read_csv(res_path)
        # Check merge key
        if 'dilirank_name' in df_res.columns:
            # Only keep columns that exist
            cols_to_use = ['dilirank_name', 'is_dili', 'tmc_prob']
            if 'cv_prob' in df_res.columns:
                cols_to_use.append('cv_prob')
            
            df = df.merge(df_res[cols_to_use], on='dilirank_name', how='left')
        else:
            print("Warning: Could not merge outcome probabilities. Assuming all relevant.")
    
    # Logic Rules
    # If Z > 1.96 (p < 0.05), the mechanism is "Active".
    # Drivers:
    threshold = 1.96
    
    mechanisms = []
    
    for i, row in df.iterrows():
        drivers = []
        if row['Z_mito'] > threshold: drivers.append(f"Mitochondrial (Z={row['Z_mito']:.1f})")
        if row['Z_bile'] > threshold: drivers.append(f"Cholestatic (Z={row['Z_bile']:.1f})")
        if row['Z_reactive'] > threshold: drivers.append(f"Reactive Met. (Z={row['Z_reactive']:.1f})")
        
        if len(drivers) == 0:
            mech = "Unexplained / Idiosyncratic"
        elif len(drivers) == 1:
            mech = drivers[0]
        else:
            mech = "Multi-Factorial: " + " + ".join(drivers)
            
        mechanisms.append(mech)
        
    df['Mechanism_Declaration'] = mechanisms
    
    # View some examples
    print("\n--- EXAMPLE DECLARATIONS ---")
    examples = df[df['Mechanism_Declaration'] != "Unexplained / Idiosyncratic"].head(10)
    for i, row in examples.iterrows():
        print(f"Drug: {row['dilirank_name']}")
        print(f"  Outcome: {'Toxic' if row.get('is_dili',0)==1 else 'Safe'}")
        print(f"  Declaration: {row['Mechanism_Declaration']}")
        print("-" * 30)
        
    # Save
    out_path = r"v:\new\h-perforatum-network-tox\research\ecnp-generalization\results\final_mechanism_reports.csv"
    df.to_csv(out_path, index=False)
    print(f"\nSaved Final Mechanism Reports to {out_path}")

if __name__ == "__main__":
    main()

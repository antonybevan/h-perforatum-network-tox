
import pandas as pd
import numpy as np
import os
import re

# --- ECNP UTILS (Simplified for brevity) ---
# In a real pipeline, we'd import the core ECNP function.
# Here, we will implement the simplified "Z-score of Subgraph" logic.

def compute_subgraph_z(targets, ontology_set, universe_size=20000):
    # targets: list of genes hit by drug
    # ontology_set: list of genes in the mechanism
    # Z = (k_obs - k_exp) / sigma_exp
    # k_obs: hits in ontology
    # k_exp: k_total * (ontology_size / universe_size)
    
    hits = [t for t in targets if t in ontology_set]
    k_obs = len(hits)
    k_total = len(targets)
    
    if k_total == 0:
        return 0.0, 0
    
    p = len(ontology_set) / universe_size
    k_exp = k_total * p
    var = k_total * p * (1 - p)
    sigma = np.sqrt(var)
    
    if sigma == 0:
        return 0.0, k_obs
        
    z = (k_obs - k_exp) / sigma
    return z, k_obs

def main():
    print("==================================================")
    print("PHASE 11: MECHANISM DECLARATION LAYER")
    print("==================================================")
    
    # 1. Load Data
    mapping_path = r"v:\new\h-perforatum-network-tox\research\ecnp-generalization\results\dilirank_drugbank_targets.csv" # Adjusted path based on 'find' result in 'results'
    # Wait, the find said 'results\dilirank_drugbank_targets.csv', let's double check path.
    # Actually, previous 'find' showed it in 'results'.
    
    if not os.path.exists(mapping_path):
        # Fallback to search if path is wrong (it might be in root or data)
        print(f"Warning: {mapping_path} not found. Checking alternate...")
        mapping_path = r"v:\new\h-perforatum-network-tox\research\ecnp-generalization\results\dilirank_drugbank_targets.csv"
    
    print(f"Loading Mapping: {mapping_path}")
    try:
        df_map = pd.read_csv(mapping_path)
    except:
        print("Error loading file.")
        return

    # Assume columns: 'DILIrank_name', 'Gene_Symbol' (or similar)
    print("Columns:", df_map.columns.tolist())
    
    # Cleaning
    # Let's assume standard 'Gene' column is present. If distinct rows per target.
    
    # Logic to identify columns
    if 'targets' in df_map.columns:
        targets_col = 'targets'
        drug_col = 'dilirank_name'
    elif 'target' in df_map.columns:
        targets_col = 'target'
        drug_col = 'dilirank_name'
    else:
        # Fallback
        targets_col = df_map.columns[-1]
        drug_col = df_map.columns[0]
        
    print(f"Grouping by {drug_col} using target col {targets_col}")
    
    drugs = df_map[drug_col].unique()
    
    # Build Universe (Observed Genes)
    # We need to handle the list strings here too for universe size if needed, 
    # but for simplicity let's assume universe is large or just use count of unique strings
    # Correct approach:
    universe_set = set()
    raw_all = df_map[targets_col].dropna().values
    for t in raw_all:
        if str(t).startswith("['") and str(t).endswith("']"):
            import ast
            try:
                parsed = ast.literal_eval(t)
                universe_set.update(parsed)
            except:
                pass
        else:
            universe_set.add(t)
            
    universe_size = len(universe_set)
    print(f"Universe Size: {universe_size}")
    
    # 2. Define Ontologies (Hardcoded Uniprot IDs)
    # We use high-confidence drivers for these mechanisms.
    
    # MITOCHONDRIAL (Complex I/V, SOD2, POLG, VDAC)
    # Source: Uniprot Human
    mito_genes = {
        'P03886', # MT-ND1
        'P03891', # MT-ND2
        'P00395', # MT-CO1
        'P04179', # SOD2 (Mangnese SOD)
        'P54098', # POLG
        'Q00059', # TFAM
        'P21796', # VDAC1
        'P45880', # VDAC2
        'Q9Y277', # VDAC3
        'P10415', # BCL2 (Mito apoptosis)
        'P55957', # BID (Mito apoptosis)
    }
    
    # CHOLESTATIC (BSEP, MDR3, MRP2, FXR, NTCP)
    bile_genes = {
        'O95342', # ABCB11 (BSEP) - Critical
        'P21439', # ABCB4 (MDR3)
        'Q92887', # ABCC2 (MRP2)
        'Q96RI1', # NR1H4 (FXR) - Master Regulator
        'Q14973', # SLC10A1 (NTCP)
        'P43358', # ABCC3 (MRP3)
        'Q9H2M9', # ABCC4 (MRP4)
    }
    
    # REACTIVE (CYP, UGT, GST) - Bioactivation Engines
    reactive_genes = {
        'P08684', # CYP3A4
        'P11712', # CYP2C9
        'P11509', # CYP2C19
        'P10635', # CYP2D6
        'P05177', # CYP1A2
        'P05181', # CYP2E1 (Toxicophore generator)
        'P22309', # UGT1A1
        'P09210', # GSTA1
        'P09488', # GSTM1
        'P08263', # GSTP1
    }
    
    print(f"  Mitochondrial Ontology Size: {len(mito_genes)}")
    print(f"  Cholestatic Ontology Size: {len(bile_genes)}")
    print(f"  Reactive Ontology Size: {len(reactive_genes)}")
    
    results = [] # Initialize results list
    for drug in drugs:
        sub = df_map[df_map[drug_col] == drug]
        # Targets are in 'target' column (already IDs)
        # If the file format uses list strings (e.g. "['P123', 'P456']"), we need to parse.
        # But 'dilirank_drugbank_targets.csv' looked like one target per row.
        # Let's handle both.
        
        raw_targets = sub[targets_col].dropna().values
        drug_targets = set()
        
        for t in raw_targets:
            # Check if it looks like a list string "['A', 'B']"
            if str(t).startswith("['") and str(t).endswith("']"):
                # Clean parse
                import ast
                try:
                    parsed = ast.literal_eval(t)
                    drug_targets.update(parsed)
                except:
                    pass
            else:
                drug_targets.add(t)
        
        z_mito, k_mito = compute_subgraph_z(list(drug_targets), mito_genes, universe_size)
        
        # NOTE: BSEP is rare. If Z is low but k_bile > 0, we should flag it.
        # Modified logic for Bile: Small set, so even 1 hit is significant.
        z_bile, k_bile = compute_subgraph_z(list(drug_targets), bile_genes, universe_size)
        
        z_reac, k_reac = compute_subgraph_z(list(drug_targets), reactive_genes, universe_size)
        
        # Dominant Mechanism
        # Special rule: Any Direct BSEP/MDR3 hit is Cholestatic regardless of Z (High Specificity)
        mech = "Unknown"
        
        drivers = []
        if z_mito > 1.0 or k_mito > 0: drivers.append("Mitochondrial")
        if z_bile > 1.0 or k_bile > 0: drivers.append("Cholestatic")
        if z_reac > 1.6: drivers.append("Reactive") # Higher bar for CYPs as they are common
        
        if not drivers:
            mech = "Unexplained"
        else:
            mech = "+".join(drivers)
            
        results.append({
            'dilirank_name': drug,
            'Z_mito': z_mito,
            'k_mito': k_mito,
            'Z_bile': z_bile,
            'k_bile': k_bile,
            'Z_reactive': z_reac,
            'k_reactive': k_reac,
            'Proposed_Mechanism': mech,
            'Max_Z_Mech': max(z_mito, z_bile, z_reac)
        })
        
    df_mech = pd.DataFrame(results)
    output_path = r"v:\new\h-perforatum-network-tox\research\ecnp-generalization\results\mechanism_features.csv"
    df_mech.to_csv(output_path, index=False)
    print(f"Saved Mechanism Features to {output_path}")
    print(df_mech['Proposed_Mechanism'].value_counts())

if __name__ == "__main__":
    main()

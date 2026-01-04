
import pandas as pd
import numpy as np

def audit_integrity():
    input_path = r"v:\new\h-perforatum-network-tox\research\ecnp-generalization\results\tmc_features_706.csv"
    try:
        df = pd.read_csv(input_path)
    except FileNotFoundError:
        print("Dataset not found!")
        return

    print("="*60)
    print("FORENSIC DATA INTEGRITY AUDIT")
    print("="*60)
    print(f"Dataset Size: {len(df)} rows")
    
    # 1. Null Check
    print("\n[1] Missing Value Check:")
    nulls = df.isnull().sum()
    if nulls.sum() == 0:
        print("    PASS: No missing values found.")
    else:
        print("    FAIL: Missing values detected:")
        print(nulls[nulls > 0])

    # 2. Known Control Check (Spot Check)
    # List of known DILI drugs to verify alignment
    # Acetaminophen: vMost-DILI-Concern -> 1
    # Amoxicillin: vMost-DILI-Concern -> 1
    # Metformin: Ambiguous/No-DILI-Concern -> 0? (Need to check specific mapping)
    
    # Let's check a few specific names if they exist
    controls = [
        {'name': 'acetaminophen', 'expect_dili': 1, 'expect_high_targets': True},
        {'name': 'amoxicillin', 'expect_dili': 1, 'expect_high_targets': False}, 
        {'name': 'ibuprofen', 'expect_dili': 1, 'expect_high_targets': True}
    ]
    
    # Normalize names for matching
    df['name_lower'] = df['dilirank_name'].astype(str).str.lower()
    
    print("\n[2] Control Spot Checks (Row Alignment):")
    alignment_issues = 0
    for ctrl in controls:
        match = df[df['name_lower'] == ctrl['name']]
        if len(match) == 0:
            # Fuzzy match?
            match = df[df['name_lower'].str.contains(ctrl['name'])]
            
        if len(match) > 0:
            row = match.iloc[0]
            label = row['is_dili']
            
            # Use k_log (log1p of mapped targets)
            k_log = row.get('k_log', 0)
            n_targets_est = np.expm1(k_log) # Reverse log1p
            
            status = "OK"
            if label != ctrl['expect_dili']:
                status = "MISMATCH (Label)"
                alignment_issues += 1
            
            # Target sanity check
            tgt_status = "OK"
            if ctrl['expect_high_targets'] and n_targets_est < 1:
                 tgt_status = "ZERO_TARGETS_WARNING"
                 # Only flag as 'issue' if we strictly expect high targets
                 # But for integrity, main thing is row alignment
            
            print(f"    {ctrl['name'].upper()}: Label={label} (Exp: {ctrl['expect_dili']}) | k_log={k_log:.2f} (Appx Targets={n_targets_est:.1f}) [{status}] {tgt_status if tgt_status != 'OK' else ''}")
        else:
            print(f"    {ctrl['name'].upper()}: Not found in dataset")

    if alignment_issues == 0:
        print("    PASS: Known drugs are aligned with labels.")
    else:
        print("    FAIL: Label mismatches detected.")

    # 3. Distribution Consistency
    print("\n[3] Feature Distribution Check:")
    # Check if 'network_coverage' is valid [0, 1]
    cov_min, cov_max = df['network_coverage'].min(), df['network_coverage'].max()
    print(f"    Network Coverage range: [{cov_min:.4f}, {cov_max:.4f}]")
    if cov_min < 0 or cov_max > 1.01:
        print("    FAIL: Coverage out of bounds!")
    else:
        print("    PASS: Coverage is valid.")
        
    # Check if 'is_dili' is binary
    labels = df['is_dili'].unique()
    print(f"    Labels present: {labels}")
    
    print("\nAudit Complete.")

if __name__ == "__main__":
    audit_integrity()

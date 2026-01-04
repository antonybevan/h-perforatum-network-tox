"""
Script: analyze_noise_quality.py
Purpose: Investigate 'why' the model fails on Tier 2a (Noisy) vs Tier 2b (Clean).
Hypothesis: Tier 2a contains "Promiscuous Targets" (e.g., CYPs, Transporters) that inflate ECNP scores for safe drugs, destroying the signal.
"""
import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import pearsonr

ROOT = Path(r'v:\new\h-perforatum-network-tox')
TIER2A_PATH = ROOT / 'research/ecnp-generalization/results/dilirank_706_with_ecnp.csv'
TIER2B_PATH = ROOT / 'research/ecnp-generalization/results/ecfp_model_results.csv'

def analyze_dataset(path, name):
    print(f"\n--- Analyzing {name} ---")
    df = pd.read_csv(path)
    
    # Filter for valid data
    if 'ecnp_z' in df.columns:
        df = df[df['ecnp_z'].notna()]
    elif 'Z' in df.columns:
        df = df[df['Z'].notna()]
        
    # 1. Target Counts (Degree of Toxicity?)
    if 'n_targets_uniprot' in df.columns:
        mean_targets = df['n_targets_uniprot'].mean()
        col = 'n_targets_uniprot'
    elif 'n_targets' in df.columns:
        mean_targets = df['n_targets'].mean()
        col = 'n_targets'
    else:
        mean_targets = 0
        col = None
        
    print(f"Mean Targets per Drug: {mean_targets:.2f}")
    
    # 2. Correlation: Network Impact vs DILI?
    # Does 'I_T' (Impact) actually correlate with toxicity?
    if 'is_dili' in df.columns and 'I_T' in df.columns:
        corr_it, _ = pearsonr(df['is_dili'], df['I_T'])
        print(f"Corr(I_T, DILI):       {corr_it:.4f}")
    
    # 3. Correlation: Network Fraction vs DILI?
    # Tier 2b had fraction ~0.33 correlation. Check 2a.
    if 'network_fraction' not in df.columns and col:
        # Approximate if missing
        df['network_fraction'] = df['k'] / df[col].replace(0, 1)
        
    if 'network_fraction' in df.columns:
        corr_frac, _ = pearsonr(df['is_dili'], df['network_fraction'])
        print(f"Corr(Fraction, DILI):  {corr_frac:.4f}")
        
    # 4. False High Scores? (Safe drugs with High Impact)
    # Safe drugs with Top 20% Impact scores
    top_20_imp = df['I_T'].quantile(0.8)
    false_alarms = df[(df['is_dili'] == 0) & (df['I_T'] > top_20_imp)]
    fa_rate = len(false_alarms) / len(df[df['is_dili'] == 0]) 
    print(f"False Alarm Rate (Safe drugs w/ High Impact): {fa_rate:.1%}")
    
    return df

# Main
df_2a = analyze_dataset(TIER2A_PATH, "Tier 2a (Noisy/Automated)")
df_2b = analyze_dataset(TIER2B_PATH, "Tier 2b (Clean/Curated)")

print("\n--- COMPARISON EXPLANATION ---")
print("If Corr(I_T, DILI) is lower in Tier 2a, the 'Impact' signal is diluted.")
print("If False Alarm Rate is higher in Tier 2a, random targets are 'lying' to the model.")

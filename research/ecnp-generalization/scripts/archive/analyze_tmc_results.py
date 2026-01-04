import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score
import numpy as np

def analyze():
    path = r"v:\new\h-perforatum-network-tox\research\ecnp-generalization\results\tmc_validation_results.csv"
    try:
        df = pd.read_csv(path)
    except FileNotFoundError:
        print("Error: Results file not found.")
        return

    # Binning by Target Count (Approximated by k_log)
    # k_log = log(n_targets_mapped + 1)
    
    # Check if k_log exists, otherwise try n_targets_mapped
    col = 'k_log' if 'k_log' in df.columns else 'n_targets_mapped'
    
    # Bin by quartiles
    try:
        df['target_bin'] = pd.qcut(df[col], q=4, labels=['Q1 (Low)', 'Q2', 'Q3', 'Q4 (High)'])
        print(f"Stratifying by {col}...")
        
        summary = df.groupby('target_bin', observed=False)[['net_uncertainty', 'net_evidence', 'tmc_prob', 'is_dili']].mean()
        print(summary)
        
        print(f"\nCorrelation ({col} vs Uncertainty):")
        print(df[[col, 'net_uncertainty']].corr())
    except Exception as e:
        print(f"Stratification failed: {e}")

    print("\n" + "="*50)
    print("PERFORMANCE METRICS (AUC)")
    print("="*50)

    # Global AUC
    y_true = df['is_dili']
    y_pred = df['tmc_prob']
    if y_true.nunique() > 1:
        auc_global = roc_auc_score(y_true, y_pred)
        ap_global = average_precision_score(y_true, y_pred)
        print(f"Global AUC (All 706 drugs): {auc_global:.4f}")
        print(f"Global PR-AUC:              {ap_global:.4f}")
    
    print("\nStratified Performance by Target Count:")
    try:
        for label in ['Q1 (Low)', 'Q2', 'Q3', 'Q4 (High)']:
            subset = df[df['target_bin'] == label]
            if len(subset) > 0 and subset['is_dili'].nunique() > 1:
                auc_sub = roc_auc_score(subset['is_dili'], subset['tmc_prob'])
                print(f"  {label} Targets (N={len(subset)}): AUC = {auc_sub:.4f}")
    except:
        pass



if __name__ == "__main__":
    analyze()

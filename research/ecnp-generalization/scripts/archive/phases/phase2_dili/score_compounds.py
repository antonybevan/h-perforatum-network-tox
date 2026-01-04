"""
ECNP Scoring at Scale - Phase 2
================================

Score all labeled compounds with ECNP and compute ROC/AUC.
"""
import sys
sys.path.insert(0, r'v:\new\h-perforatum-network-tox\research\ecnp-closed-form\src\core')

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score, roc_curve, confusion_matrix
import ast

from ecnp_optimized import ECNPOptimized
from ecnp_permutation_test import ECNPPermutationTest, PermutationConfig

# Paths
RESEARCH_ROOT = Path(r'v:\new\h-perforatum-network-tox\research')
LABELED_FILE = RESEARCH_ROOT / 'ecnp-generalization' / 'data' / 'curated' / 'labeled_compounds_for_validation.csv'
OUTPUT_DIR = RESEARCH_ROOT / 'ecnp-generalization' / 'results' / 'phase2'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Load ECNP
print("Loading ECNP...")
project_root = Path(r'v:\new\h-perforatum-network-tox')
ecnp = ECNPOptimized(project_root)
print(f"Loaded: {ecnp.n_nodes} nodes, {len(ecnp.dili_genes)} DILI genes")

# Load compounds
print("\nLoading labeled compounds...")
df = pd.read_csv(LABELED_FILE)
print(f"Compounds to score: {len(df)}")

# Parse target lists (stored as string representation of list)
def parse_targets(targets_str):
    """Parse target list from CSV string."""
    try:
        return ast.literal_eval(targets_str)
    except:
        return []

# Score each compound
print("\nScoring compounds with ECNP...")
results = []

for idx, row in df.iterrows():
    targets = parse_targets(row['targets'])
    
    if len(targets) < 3:
        continue
    
    # Layer 1: Fast Z-score
    result = ecnp.compute(targets)
    
    if result['status'].value == 'success':
        results.append({
            'drugbank_id': row['drugbank_id'],
            'drug_name': row['drug_name'],
            'n_targets': row['n_targets'],
            'dili_label': row['dili_label'],
            'is_dili': 1 if row['dili_label'] in ['Most-DILI', 'Less-DILI'] else 0,
            'Z': result['Z'],
            'I_T': result['I_T'],
            'mu_T': result['mu_T'],
            'pool_size': result['pool_size']
        })
    
    if (idx + 1) % 10 == 0:
        print(f"  Scored {idx + 1}/{len(df)}...")

print(f"Successfully scored: {len(results)}")

# Create results dataframe
results_df = pd.DataFrame(results)
results_df.to_csv(OUTPUT_DIR / 'ecnp_scores.csv', index=False)
print(f"Saved: {OUTPUT_DIR / 'ecnp_scores.csv'}")

# ROC/AUC Analysis
print("\n" + "="*50)
print("ROC/AUC ANALYSIS")
print("="*50)

y_true = results_df['is_dili'].values
y_score = results_df['Z'].values

# Check we have both classes
n_pos = y_true.sum()
n_neg = len(y_true) - n_pos
print(f"Positive (DILI): {n_pos}")
print(f"Negative (No-DILI): {n_neg}")

if n_pos > 0 and n_neg > 0:
    auc = roc_auc_score(y_true, y_score)
    print(f"\nAUC: {auc:.3f}")
    
    # ROC curve
    fpr, tpr, thresholds = roc_curve(y_true, y_score)
    
    # Find optimal threshold (Youden's J)
    j_scores = tpr - fpr
    best_idx = np.argmax(j_scores)
    best_threshold = thresholds[best_idx]
    
    print(f"Optimal threshold: Z = {best_threshold:.2f}")
    print(f"At threshold:")
    print(f"  Sensitivity: {tpr[best_idx]:.3f}")
    print(f"  Specificity: {1 - fpr[best_idx]:.3f}")
    
    # Confusion matrix at optimal threshold
    y_pred = (y_score >= best_threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    print(f"\nConfusion Matrix (at Z >= {best_threshold:.2f}):")
    print(f"  TP: {tp}, FP: {fp}")
    print(f"  FN: {fn}, TN: {tn}")
    
    # Save ROC data
    roc_df = pd.DataFrame({'fpr': fpr, 'tpr': tpr, 'threshold': thresholds})
    roc_df.to_csv(OUTPUT_DIR / 'roc_curve.csv', index=False)
else:
    print("WARNING: Need both positive and negative samples for AUC")
    auc = None

# Summary table
print("\n" + "="*50)
print("COMPOUND RANKING")
print("="*50)
print("\nTop 10 by Z-score:")
print(results_df.nlargest(10, 'Z')[['drug_name', 'dili_label', 'Z', 'n_targets']].to_string(index=False))

print("\nBottom 10 by Z-score:")
print(results_df.nsmallest(10, 'Z')[['drug_name', 'dili_label', 'Z', 'n_targets']].to_string(index=False))

# By DILI label
print("\n" + "="*50)
print("Z-SCORE BY DILI LABEL")
print("="*50)
for label in results_df['dili_label'].unique():
    subset = results_df[results_df['dili_label'] == label]
    print(f"\n{label} (n={len(subset)}):")
    print(f"  Mean Z: {subset['Z'].mean():.2f}")
    print(f"  Std Z: {subset['Z'].std():.2f}")
    print(f"  Range: [{subset['Z'].min():.2f}, {subset['Z'].max():.2f}]")

# Final summary
print("\n" + "="*50)
print("PHASE 2 SUMMARY")
print("="*50)
print(f"Compounds scored: {len(results_df)}")
print(f"AUC: {auc:.3f}" if auc else "AUC: N/A (need both classes)")
print(f"Mean Z (DILI+): {results_df[results_df['is_dili']==1]['Z'].mean():.2f}")
print(f"Mean Z (DILI-): {results_df[results_df['is_dili']==0]['Z'].mean():.2f}" if n_neg > 0 else "Mean Z (DILI-): N/A")

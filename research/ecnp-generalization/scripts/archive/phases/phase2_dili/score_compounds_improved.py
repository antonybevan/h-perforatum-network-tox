"""
ECNP Scoring on Improved Dataset
=================================

Score the 207 matched compounds and compute ROC/AUC.
"""
import sys
sys.path.insert(0, r'v:\new\h-perforatum-network-tox\research\ecnp-closed-form\src\core')

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score, roc_curve, confusion_matrix
import ast

from ecnp_optimized import ECNPOptimized

RESEARCH_ROOT = Path(r'v:\new\h-perforatum-network-tox\research')
LABELED_FILE = RESEARCH_ROOT / 'ecnp-generalization' / 'data' / 'curated' / 'labeled_compounds_improved.csv'
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
print(f"Columns: {df.columns.tolist()}")

# Check dilirank values
print(f"\nDILIrank distribution:")
print(df['dilirank'].value_counts())

# Create binary labels
# vMost-DILI-concern and vLess-DILI-concern = positive
# vNo-DILI-concern = negative
def get_binary_label(x):
    if pd.isna(x):
        return np.nan
    x = str(x).lower()
    if 'most' in x or 'less' in x:
        return 1
    elif 'no' in x:
        return 0
    else:
        return np.nan

df['is_dili'] = df['dilirank'].apply(get_binary_label)

# Filter to labeled
df_labeled = df[df['is_dili'].notna()].copy()
df_labeled['is_dili'] = df_labeled['is_dili'].astype(int)

print(f"\nAfter binary labeling:")
print(f"  Total: {len(df_labeled)}")
print(f"  DILI+: {(df_labeled['is_dili']==1).sum()}")
print(f"  DILI-: {(df_labeled['is_dili']==0).sum()}")

# Parse target lists
def parse_targets(targets_str):
    try:
        return ast.literal_eval(targets_str)
    except:
        return []

# Score each compound
print("\nScoring compounds with ECNP...")
results = []

for idx, row in df_labeled.iterrows():
    targets = parse_targets(row['targets'])
    
    if len(targets) < 3:
        continue
    
    result = ecnp.compute(targets)
    
    if result['status'].value == 'success':
        results.append({
            'drugbank_id': row['drugbank_id'],
            'drug_name': row['drug_name'],
            'dilirank_name': row['dilirank_name'],
            'n_targets': row['n_targets'],
            'dilirank': row['dilirank'],
            'is_dili': row['is_dili'],
            'Z': result['Z'],
            'I_T': result['I_T'],
            'mu_T': result['mu_T'],
            'pool_size': result['pool_size']
        })
    
    if (idx + 1) % 50 == 0:
        print(f"  Progress: {len(results)} scored...")

print(f"Successfully scored: {len(results)}")

results_df = pd.DataFrame(results)
results_df.to_csv(OUTPUT_DIR / 'ecnp_scores_improved.csv', index=False)
print(f"Saved: {OUTPUT_DIR / 'ecnp_scores_improved.csv'}")

# ROC/AUC Analysis
print("\n" + "="*60)
print("ROC/AUC ANALYSIS")
print("="*60)

y_true = results_df['is_dili'].values
y_score = results_df['Z'].values

n_pos = y_true.sum()
n_neg = len(y_true) - n_pos
print(f"Positive (DILI+): {n_pos}")
print(f"Negative (DILI-): {n_neg}")

if n_pos > 0 and n_neg > 0:
    auc = roc_auc_score(y_true, y_score)
    print(f"\nAUC: {auc:.3f}")
    
    fpr, tpr, thresholds = roc_curve(y_true, y_score)
    
    # Optimal threshold (Youden's J)
    j_scores = tpr - fpr
    best_idx = np.argmax(j_scores)
    best_threshold = thresholds[best_idx]
    
    print(f"Optimal threshold: Z = {best_threshold:.2f}")
    print(f"At threshold:")
    print(f"  Sensitivity: {tpr[best_idx]:.3f}")
    print(f"  Specificity: {1 - fpr[best_idx]:.3f}")
    
    # Confusion matrix
    y_pred = (y_score >= best_threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    print(f"\nConfusion Matrix:")
    print(f"  TP: {tp}, FP: {fp}")
    print(f"  FN: {fn}, TN: {tn}")
    
    # Save ROC
    roc_df = pd.DataFrame({'fpr': fpr, 'tpr': tpr, 'threshold': thresholds})
    roc_df.to_csv(OUTPUT_DIR / 'roc_curve_improved.csv', index=False)
else:
    auc = None
    print("ERROR: Need both classes for AUC")

# Summary by DILI category
print("\n" + "="*60)
print("Z-SCORE BY DILI CATEGORY")
print("="*60)

for dilirank_val, label in [(1.0, 'Most-DILI'), (2.0, 'Less-DILI'), (3.0, 'No-DILI')]:
    subset = results_df[results_df['dilirank'] == dilirank_val]
    if len(subset) > 0:
        print(f"\n{label} (n={len(subset)}):")
        print(f"  Mean Z: {subset['Z'].mean():.2f}")
        print(f"  Std Z: {subset['Z'].std():.2f}")
        print(f"  Range: [{subset['Z'].min():.2f}, {subset['Z'].max():.2f}]")

# Final summary
print("\n" + "="*60)
print("PHASE 2 SUMMARY (IMPROVED)")
print("="*60)
print(f"Compounds scored: {len(results_df)}")
print(f"AUC: {auc:.3f}" if auc else "AUC: N/A")
print(f"Target n>=150: {'YES' if len(results_df) >= 150 else 'NO'} ({len(results_df)})")

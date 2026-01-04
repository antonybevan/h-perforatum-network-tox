"""
Re-run ECNP Validation with Full DILI Module (506 genes)
=========================================================
"""
import sys
sys.path.insert(0, r'v:\new\h-perforatum-network-tox\research\ecnp-closed-form\src\core')

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score
import ast

ROOT = Path(r'v:\new\h-perforatum-network-tox')
RESEARCH_DIR = ROOT / 'research' / 'ecnp-closed-form' / 'data'

# Load influence matrix and new DILI vector
print("Loading data with FULL DILI module (506 genes)...")
npz = np.load(RESEARCH_DIR / 'influence_matrix_900.npz', allow_pickle=True)
M = npz['M']
node_list = npz['node_list'].tolist()
node_to_idx = {g: i for i, g in enumerate(node_list)}
n_nodes = len(node_list)

# Load NEW DILI influence vector
m_df = pd.read_csv(RESEARCH_DIR / 'dili_influence_vector_900_full.csv')
m_vector = m_df.set_index('gene')['dili_influence']
m_array = np.array([m_vector.get(g, 0) for g in node_list])

# Load NEW DILI genes
dili_df = pd.read_csv(RESEARCH_DIR / 'dili_genes_full.csv')
dili_genes = set(dili_df['gene_name'].tolist())
print(f"DILI module: {len(dili_genes)} genes")

# Load labeled compounds
compounds_df = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'data' / 'curated' / 'labeled_compounds_improved.csv')
print(f"Compounds to score: {len(compounds_df)}")

# Parse target lists
def parse_targets(t):
    try:
        return ast.literal_eval(t) if isinstance(t, str) else []
    except:
        return []

# Simple ECNP computation
def compute_ecnp_simple(targets, m_array, node_to_idx, node_list):
    """Simplified ECNP: just sum DILI influence of targets, normalize."""
    idx = [node_to_idx[t] for t in targets if t in node_to_idx]
    if len(idx) < 3:
        return None
    
    I_T = sum(m_array[i] for i in idx)
    k = len(idx)
    
    # Compute pool mean (all nodes in network)
    pool_mean = np.mean(m_array)
    pool_std = np.std(m_array)
    
    mu_T = k * pool_mean
    sigma_T = np.sqrt(k) * pool_std
    
    Z = (I_T - mu_T) / sigma_T if sigma_T > 0 else 0
    
    return {
        'Z': Z,
        'I_T': I_T,
        'k': k
    }

# Create binary labels
def get_binary_label(x):
    x = str(x).lower()
    if 'most' in x or 'less' in x:
        return 1
    elif 'no' in x:
        return 0
    else:
        return np.nan

compounds_df['is_dili'] = compounds_df['dilirank'].apply(get_binary_label)
labeled = compounds_df[compounds_df['is_dili'].notna()].copy()
labeled['is_dili'] = labeled['is_dili'].astype(int)

# Score all compounds
print("\nScoring compounds with 506-gene DILI module...")
results = []
for idx, row in labeled.iterrows():
    targets = parse_targets(row['targets'])
    result = compute_ecnp_simple(targets, m_array, node_to_idx, node_list)
    
    if result:
        # Count direct hits (in new DILI module)
        n_direct_hits = len(set(targets) & dili_genes)
        
        results.append({
            'drug_name': row['drug_name'],
            'is_dili': row['is_dili'],
            'Z': result['Z'],
            'k': result['k'],
            'n_direct_hits': n_direct_hits,
            'pct_direct_hits': n_direct_hits / len(targets) * 100
        })

results_df = pd.DataFrame(results)
print(f"Scored: {len(results_df)}")

# Compare direct hits with new module
print("\n" + "="*60)
print("DIRECT HIT ANALYSIS (506-gene module)")
print("="*60)
print(f"\nMean direct hits: {results_df['n_direct_hits'].mean():.2f}")
print(f"Max direct hits: {results_df['n_direct_hits'].max()}")

# By DILI status
for status in [1, 0]:
    label = "DILI+" if status == 1 else "DILI-"
    subset = results_df[results_df['is_dili'] == status]
    print(f"\n{label} (n={len(subset)}):")
    print(f"  Mean direct hits: {subset['n_direct_hits'].mean():.2f}")
    print(f"  Mean Z: {subset['Z'].mean():.2f}")

# ROC/AUC
print("\n" + "="*60)
print("ROC/AUC ANALYSIS")
print("="*60)

y_true = results_df['is_dili'].values
y_score = results_df['Z'].values

n_pos = y_true.sum()
n_neg = len(y_true) - n_pos
print(f"Positive: {n_pos}, Negative: {n_neg}")

if n_pos > 0 and n_neg > 0:
    auc = roc_auc_score(y_true, y_score)
    print(f"\nAUC: {auc:.3f}")
    
    # Compare to old module
    print(f"\nOld AUC (82 genes): 0.624")
    print(f"New AUC (506 genes): {auc:.3f}")
    print(f"Improvement: {(auc - 0.624):.3f}")
    
    if auc >= 0.75:
        print("\n** TARGET MET: AUC >= 0.75 **")
    elif auc >= 0.70:
        print("\n** ACCEPTABLE: AUC >= 0.70 **")

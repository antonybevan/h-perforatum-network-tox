"""
Scope-Gated ECNP Validation
============================

ECNP applies to the NETWORK-MEDIATED regime only.

Filter criteria:
1. >=1 direct DILI gene hit, OR
2. >=1 CYP regulator / transcription factor target

Then compute AUC on this subset only.
"""
import sys
sys.path.insert(0, r'v:\new\h-perforatum-network-tox\research\ecnp-closed-form\src\core')

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score, roc_curve
import ast

ROOT = Path(r'v:\new\h-perforatum-network-tox')

# Load 82-gene DILI module (the original, validated one)
dili_82 = pd.read_csv(ROOT / 'data' / 'processed' / 'dili_900_lcc.csv')
dili_genes_82 = set(dili_82['gene_name'].tolist())
print(f"DILI module: {len(dili_genes_82)} genes")

# CYP and transcription factor genes (network regulators)
# These are canonical drug metabolism and regulatory genes
CYP_GENES = {
    'CYP1A1', 'CYP1A2', 'CYP1B1', 'CYP2A6', 'CYP2B6', 'CYP2C8', 'CYP2C9',
    'CYP2C19', 'CYP2D6', 'CYP2E1', 'CYP3A4', 'CYP3A5', 'CYP3A7',
    'CYP2J2', 'CYP4A11', 'CYP4F2', 'CYP4F3', 'CYP4F11', 'CYP4F12',
    'CYP7A1', 'CYP8B1', 'CYP27A1', 'CYP51A1'
}

# Key liver transcription factors and nuclear receptors
TF_REGULATORS = {
    'PXR', 'NR1I2',  # PXR
    'CAR', 'NR1I3',  # CAR
    'HNF4A', 'HNF1A', 'HNF1B',  # Hepatocyte nuclear factors
    'PPARA', 'PPARG', 'PPARD',  # PPARs
    'RXRA', 'RXRB', 'RXRG',  # RXRs
    'LXR', 'NR1H3', 'NR1H2',  # LXRs
    'FXR', 'NR1H4',  # FXR
    'AHR',  # Aryl hydrocarbon receptor
    'CEBPA', 'CEBPB',  # C/EBPs
    'FOXA1', 'FOXA2', 'FOXA3',  # FOXAs
    'SREBF1', 'SREBF2',  # SREBPs
    'XBP1',  # UPR
}

# Transporter genes (for cholestatic regime tagging)
TRANSPORTERS = {
    'ABCB1', 'ABCB4', 'ABCB11',  # MDR1, MDR3, BSEP
    'ABCC1', 'ABCC2', 'ABCC3', 'ABCC4',  # MRPs
    'ABCG2',  # BCRP
    'SLC22A1', 'SLC22A6', 'SLC22A7', 'SLC22A8',  # OATs
    'SLCO1B1', 'SLCO1B3', 'SLCO2B1',  # OATPs
    'SLC10A1', 'SLC10A2',  # NTCPs
}

NETWORK_MEDIATED_GENES = dili_genes_82 | CYP_GENES | TF_REGULATORS
print(f"Network-mediated gene set: {len(NETWORK_MEDIATED_GENES)} genes")

# Load scored compounds
scores_df = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'results' / 'phase2' / 'ecnp_scores_improved.csv')
compounds_df = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'data' / 'curated' / 'labeled_compounds_improved.csv')

def parse_targets(t):
    try:
        return ast.literal_eval(t) if isinstance(t, str) else []
    except:
        return []

compounds_df['target_list'] = compounds_df['targets'].apply(parse_targets)

# Count hits in each regime
def count_hits(targets, gene_set):
    return len(set(targets) & gene_set)

compounds_df['n_dili_hits'] = compounds_df['target_list'].apply(lambda t: count_hits(t, dili_genes_82))
compounds_df['n_cyp_tf_hits'] = compounds_df['target_list'].apply(lambda t: count_hits(t, CYP_GENES | TF_REGULATORS))
compounds_df['n_transporter_hits'] = compounds_df['target_list'].apply(lambda t: count_hits(t, TRANSPORTERS))
compounds_df['n_network_hits'] = compounds_df['target_list'].apply(lambda t: count_hits(t, NETWORK_MEDIATED_GENES))

# Merge with scores
merged = scores_df.merge(compounds_df[['drugbank_id', 'n_dili_hits', 'n_cyp_tf_hits', 'n_transporter_hits', 'n_network_hits']], 
                          on='drugbank_id', how='left')

print(f"\nTotal compounds scored: {len(merged)}")

# SCOPE FILTER: Network-mediated regime
# >=1 DILI hit OR >=1 CYP/TF hit
scope_mask = (merged['n_dili_hits'] >= 1) | (merged['n_cyp_tf_hits'] >= 1)
in_scope = merged[scope_mask].copy()
out_of_scope = merged[~scope_mask].copy()

print(f"\n{'='*60}")
print("SCOPE-GATED ANALYSIS")
print('='*60)
print(f"\nIN SCOPE (network-mediated): {len(in_scope)} compounds")
print(f"OUT OF SCOPE: {len(out_of_scope)} compounds")

# By DILI status
for subset, name in [(in_scope, 'IN SCOPE'), (out_of_scope, 'OUT OF SCOPE')]:
    n_pos = (subset['is_dili'] == 1).sum()
    n_neg = (subset['is_dili'] == 0).sum()
    print(f"\n{name}:")
    print(f"  DILI+: {n_pos}")
    print(f"  DILI-: {n_neg}")

# ROC/AUC for IN-SCOPE compounds
print(f"\n{'='*60}")
print("AUC: IN-SCOPE vs ALL")
print('='*60)

# All compounds (baseline)
y_true_all = merged['is_dili'].values
y_score_all = merged['Z'].values
auc_all = roc_auc_score(y_true_all, y_score_all)
print(f"\nALL compounds (n={len(merged)}): AUC = {auc_all:.3f}")

# In-scope only
y_true_scope = in_scope['is_dili'].values
y_score_scope = in_scope['Z'].values
n_pos = y_true_scope.sum()
n_neg = len(y_true_scope) - n_pos

if n_pos >= 3 and n_neg >= 3:
    auc_scope = roc_auc_score(y_true_scope, y_score_scope)
    print(f"IN-SCOPE only (n={len(in_scope)}): AUC = {auc_scope:.3f}")
    print(f"\nImprovement: {auc_scope - auc_all:+.3f}")
    
    if auc_scope >= 0.75:
        print("\n** TARGET MET: AUC >= 0.75 in network-mediated regime! **")
    elif auc_scope >= 0.70:
        print("\n** ACCEPTABLE: AUC >= 0.70 in network-mediated regime **")
    
    # Details for in-scope
    print(f"\nIN-SCOPE breakdown:")
    print(f"  Mean Z (DILI+): {in_scope[in_scope['is_dili']==1]['Z'].mean():.2f}")
    print(f"  Mean Z (DILI-): {in_scope[in_scope['is_dili']==0]['Z'].mean():.2f}")
    print(f"  Mean DILI hits (DILI+): {in_scope[in_scope['is_dili']==1]['n_dili_hits'].mean():.2f}")
    print(f"  Mean DILI hits (DILI-): {in_scope[in_scope['is_dili']==0]['n_dili_hits'].mean():.2f}")
else:
    print(f"IN-SCOPE: Too few samples (pos={n_pos}, neg={n_neg})")

# Out-of-scope AUC (should be poor)
y_true_out = out_of_scope['is_dili'].values
y_score_out = out_of_scope['Z'].values
n_pos_out = y_true_out.sum()
n_neg_out = len(y_true_out) - n_pos_out

if n_pos_out >= 3 and n_neg_out >= 3:
    auc_out = roc_auc_score(y_true_out, y_score_out)
    print(f"\nOUT-OF-SCOPE (n={len(out_of_scope)}): AUC = {auc_out:.3f}")
    print("  (Expected to be poor - ECNP doesn't apply here)")

# Stratify further by hit count
print(f"\n{'='*60}")
print("AUC BY DILI HIT COUNT")
print('='*60)

for min_hits in [0, 1, 2, 3]:
    subset = merged[merged['n_dili_hits'] >= min_hits]
    if len(subset) < 20:
        continue
    
    y_true = subset['is_dili'].values
    y_score = subset['Z'].values
    n_p = y_true.sum()
    n_n = len(y_true) - n_p
    
    if n_p >= 3 and n_n >= 3:
        auc = roc_auc_score(y_true, y_score)
        print(f">=  {min_hits} DILI hits (n={len(subset)}, pos={n_p}, neg={n_n}): AUC = {auc:.3f}")

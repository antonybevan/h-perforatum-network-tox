"""
Complete Regime-Aware DILI Prediction Pipeline
==============================================

Implements the 4-layer architecture:
- Layer 0: Exposure gate (LogP, dose, half-life)
- Layer 1: Soft mechanism classifier (regime probabilities)
- Layer 2: Regime-specific models (ECNP for network regime)
- Layer 3: Max-risk aggregation with explanation

Success metrics: Regime-specific AUC (primary), Global AUC (secondary)
"""
import sys
sys.path.insert(0, r'v:\new\h-perforatum-network-tox\research\ecnp-closed-form\src\core')

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
import ast
import warnings
warnings.filterwarnings('ignore')

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

ROOT = Path(r'v:\new\h-perforatum-network-tox')

# =============================================================================
# DATA LOADING
# =============================================================================

print("Loading data...")
compounds = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'data' / 'curated' / 'labeled_compounds_improved.csv')
scores = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'results' / 'phase2' / 'ecnp_scores_improved.csv')

# Load DILI and target gene sets
dili_genes = set(pd.read_csv(ROOT / 'data' / 'processed' / 'dili_900_lcc.csv')['gene_name'])

CYP_GENES = {'CYP1A1', 'CYP1A2', 'CYP1B1', 'CYP2A6', 'CYP2B6', 'CYP2C8', 'CYP2C9',
             'CYP2C19', 'CYP2D6', 'CYP2E1', 'CYP3A4', 'CYP3A5', 'CYP3A7'}

TF_GENES = {'NR1I2', 'NR1I3', 'HNF4A', 'HNF1A', 'PPARA', 'PPARG', 'RXRA', 'AHR', 
            'CEBPA', 'CEBPB', 'FOXA1', 'FOXA2', 'SREBF1', 'SREBF2'}

TRANSPORTERS = {'ABCB1', 'ABCB4', 'ABCB11', 'ABCC1', 'ABCC2', 'ABCC3', 
                'ABCG2', 'SLC22A1', 'SLCO1B1', 'SLCO1B3', 'SLC10A1'}

NETWORK_GENES = dili_genes | CYP_GENES | TF_GENES

print(f"Compounds: {len(compounds)}, Scores: {len(scores)}")

# Merge
df = scores.merge(compounds[['drugbank_id', 'smiles', 'targets']], on='drugbank_id', how='left')
df = df[df['is_dili'].notna()].copy()
df['is_dili'] = df['is_dili'].astype(int)
print(f"Labeled compounds: {len(df)}")

def parse_targets(t):
    try:
        return ast.literal_eval(t) if isinstance(t, str) else []
    except:
        return []

df['target_list'] = df['targets'].apply(parse_targets)

# =============================================================================
# LAYER 0: EXPOSURE FEATURES
# =============================================================================

print("\nLayer 0: Computing exposure features...")

if RDKIT_AVAILABLE:
    def compute_chemistry(smiles):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                return {
                    'logp': Descriptors.MolLogP(mol),
                    'mw': Descriptors.MolWt(mol),
                    'hbd': Lipinski.NumHDonors(mol),
                    'hba': Lipinski.NumHAcceptors(mol),
                    'tpsa': Descriptors.TPSA(mol)
                }
        except:
            pass
        return {'logp': np.nan, 'mw': np.nan, 'hbd': np.nan, 'hba': np.nan, 'tpsa': np.nan}
    
    chem_features = df['smiles'].apply(compute_chemistry).apply(pd.Series)
    df = pd.concat([df, chem_features], axis=1)
    
    # Exposure flags
    df['high_logp'] = (df['logp'] > 3).astype(int)
    df['lipinski_violations'] = (
        (df['mw'] > 500).astype(int) + 
        (df['logp'] > 5).astype(int) + 
        (df['hbd'] > 5).astype(int) + 
        (df['hba'] > 10).astype(int)
    )

# =============================================================================
# LAYER 1: SOFT MECHANISM CLASSIFIER
# =============================================================================

print("Layer 1: Computing regime probabilities...")

def compute_regime_weights(targets):
    """Output soft probabilities for each regime."""
    target_set = set(targets)
    
    n_network = len(target_set & NETWORK_GENES)
    n_transporter = len(target_set & TRANSPORTERS)
    n_total = len(targets)
    
    # Base weights from target composition
    if n_total == 0:
        return {'network': 0.33, 'reactive': 0.33, 'cholestatic': 0.34}
    
    # Network weight: based on DILI/CYP/TF hits
    w_network = min(1.0, n_network / max(n_total, 1) * 2)
    
    # Cholestatic weight: based on transporter hits
    w_cholestatic = min(1.0, n_transporter / max(n_total, 1) * 3)
    
    # Reactive weight: complement (will be refined with structural alerts)
    w_reactive = max(0.1, 1 - w_network - w_cholestatic)
    
    # Normalize
    total = w_network + w_cholestatic + w_reactive
    return {
        'network': w_network / total,
        'reactive': w_reactive / total,
        'cholestatic': w_cholestatic / total
    }

regime_weights = df['target_list'].apply(compute_regime_weights).apply(pd.Series)
df['w_network'] = regime_weights['network']
df['w_reactive'] = regime_weights['reactive']
df['w_cholestatic'] = regime_weights['cholestatic']

# ECNP eligibility: high network weight
df['ecnp_eligible'] = (df['w_network'] >= 0.3).astype(int)

# =============================================================================
# LAYER 2: REGIME-SPECIFIC SCORES
# =============================================================================

print("Layer 2: Computing regime-specific scores...")

# Network regime: ECNP Z-score (already have)
df['score_network'] = df['Z']

# Reactive regime: LogP + chemistry (simplified structural alert proxy)
if RDKIT_AVAILABLE:
    # High LogP + high molecular weight = reactive metabolite risk
    df['score_reactive'] = (df['logp'] / 5 + df['mw'] / 500) / 2
    df['score_reactive'] = df['score_reactive'].clip(0, 1)
else:
    df['score_reactive'] = 0.5

# Cholestatic regime: transporter hit count normalized
df['n_transporter'] = df['target_list'].apply(lambda t: len(set(t) & TRANSPORTERS))
df['score_cholestatic'] = df['n_transporter'] / 3  # Normalize to ~0-1 range

# =============================================================================
# LAYER 3: MAX-RISK AGGREGATION
# =============================================================================

print("Layer 3: Aggregating risk with regime weights...")

# Weighted combination of regime scores
df['weighted_score'] = (
    df['w_network'] * df['score_network'].clip(0, 1) +
    df['w_reactive'] * df['score_reactive'] +
    df['w_cholestatic'] * df['score_cholestatic']
)

# Max-risk: take highest regime score
df['max_risk_score'] = df[['score_network', 'score_reactive', 'score_cholestatic']].max(axis=1)

# Explanation: which regime triggered
df['triggered_regime'] = df[['score_network', 'score_reactive', 'score_cholestatic']].idxmax(axis=1)
df['triggered_regime'] = df['triggered_regime'].str.replace('score_', '')

# =============================================================================
# VALIDATION
# =============================================================================

print("\n" + "="*60)
print("VALIDATION: REGIME-SPECIFIC AUC")
print("="*60)

valid_mask = df['logp'].notna() if RDKIT_AVAILABLE else np.ones(len(df), dtype=bool)
df_valid = df[valid_mask].copy()

y_true = df_valid['is_dili'].values

# Overall AUC for different scores
auc_results = {}

for score_col, name in [
    ('Z', 'ECNP Z-score alone'),
    ('logp', 'LogP alone'),
    ('weighted_score', 'Weighted regime score'),
    ('max_risk_score', 'Max-risk score'),
]:
    if score_col in df_valid.columns and df_valid[score_col].notna().all():
        y_score = df_valid[score_col].values
        auc = roc_auc_score(y_true, y_score)
        auc_results[name] = auc
        print(f"{name}: AUC = {auc:.3f}")

# Regime-specific AUC
print("\n" + "-"*40)
print("REGIME-SPECIFIC AUC (Primary Metrics)")
print("-"*40)

# Network regime (ECNP eligible compounds)
ecnp_eligible = df_valid[df_valid['ecnp_eligible'] == 1]
if len(ecnp_eligible) >= 20:
    y_true_net = ecnp_eligible['is_dili'].values
    y_score_net = ecnp_eligible['score_network'].values
    n_pos = y_true_net.sum()
    n_neg = len(y_true_net) - n_pos
    if n_pos >= 3 and n_neg >= 3:
        auc_network = roc_auc_score(y_true_net, y_score_net)
        print(f"\nNetwork regime (ECNP eligible, n={len(ecnp_eligible)}):")
        print(f"  AUC = {auc_network:.3f}")
        print(f"  DILI+: {n_pos}, DILI-: {n_neg}")

# By triggered regime
print("\n" + "-"*40)
print("AUC BY TRIGGERED REGIME")
print("-"*40)

for regime in ['network', 'reactive', 'cholestatic']:
    subset = df_valid[df_valid['triggered_regime'] == regime]
    if len(subset) >= 10:
        y_t = subset['is_dili'].values
        y_s = subset['max_risk_score'].values
        n_p = y_t.sum()
        n_n = len(y_t) - n_p
        if n_p >= 3 and n_n >= 3:
            auc = roc_auc_score(y_t, y_s)
            print(f"{regime}: n={len(subset)}, AUC = {auc:.3f}")
        else:
            print(f"{regime}: n={len(subset)}, too few in one class")
    else:
        print(f"{regime}: n={len(subset)}, skipped (< 10)")

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "="*60)
print("PIPELINE SUMMARY")
print("="*60)

print(f"\nCompounds modeled: {len(df_valid)}")
print(f"ECNP eligible: {df_valid['ecnp_eligible'].sum()}")

print(f"\nRegime distribution:")
print(df_valid['triggered_regime'].value_counts())

print(f"\nBest single score: {max(auc_results.items(), key=lambda x: x[1])}")
print(f"\nTarget: Network-regime AUC >= 0.75")

# Save results
output_file = ROOT / 'research' / 'ecnp-generalization' / 'results' / 'pipeline_results.csv'
df_valid.to_csv(output_file, index=False)
print(f"\nSaved: {output_file}")

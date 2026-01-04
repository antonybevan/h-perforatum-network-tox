"""
ML-Enhanced Regime-Aware DILI Pipeline
=======================================

Architecture:
- Layer 0: Exposure gate (rule-based)
- Layer 1: ML regime classifier (learns when to ask ECNP)
- Layer 2: Regime experts (ECNP for network, chemistry for reactive)
- Layer 3: Weighted risk aggregation

Validation:
- Ablation: ML without ECNP, ECNP without ML, Full model
- Stratified by regime
- Permutation tests preserved for ECNP
"""
import sys
sys.path.insert(0, r'v:\new\h-perforatum-network-tox\research\ecnp-closed-form\src\core')

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score, brier_score_loss
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.linear_model import LogisticRegression
import ast
import warnings
warnings.filterwarnings('ignore')

try:
    from xgboost import XGBClassifier
    XGBOOST_AVAILABLE = True
except ImportError:
    XGBOOST_AVAILABLE = False
    print("XGBoost not available, using logistic regression")

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

# Gene sets
dili_genes = set(pd.read_csv(ROOT / 'data' / 'processed' / 'dili_900_lcc.csv')['gene_name'])
CYP_GENES = {'CYP1A1', 'CYP1A2', 'CYP2C9', 'CYP2C19', 'CYP2D6', 'CYP3A4', 'CYP3A5'}
TF_GENES = {'NR1I2', 'NR1I3', 'HNF4A', 'PPARA', 'RXRA', 'AHR'}
TRANSPORTERS = {'ABCB1', 'ABCB11', 'ABCC2', 'SLCO1B1', 'SLC22A1'}
NETWORK_GENES = dili_genes | CYP_GENES | TF_GENES

# Merge
df = scores.merge(compounds[['drugbank_id', 'smiles', 'targets']], on='drugbank_id', how='left')
df = df[df['is_dili'].notna()].copy()
df['is_dili'] = df['is_dili'].astype(int)

def parse_targets(t):
    try:
        return ast.literal_eval(t) if isinstance(t, str) else []
    except:
        return []

df['target_list'] = df['targets'].apply(parse_targets)
print(f"Compounds: {len(df)}")

# =============================================================================
# FEATURE ENGINEERING
# =============================================================================

print("\nEngineering features...")

# Block A: Exposure / Chemistry
if RDKIT_AVAILABLE:
    def get_chem(smiles):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                return {
                    'logp': Descriptors.MolLogP(mol),
                    'mw': Descriptors.MolWt(mol),
                    'hbd': Lipinski.NumHDonors(mol),
                    'hba': Lipinski.NumHAcceptors(mol),
                    'tpsa': Descriptors.TPSA(mol),
                    'rotatable': Descriptors.NumRotatableBonds(mol)
                }
        except:
            pass
        return {k: np.nan for k in ['logp', 'mw', 'hbd', 'hba', 'tpsa', 'rotatable']}
    
    chem = df['smiles'].apply(get_chem).apply(pd.Series)
    df = pd.concat([df, chem], axis=1)

# Block B: Target topology
df['k'] = df['n_targets']
df['n_dili_hits'] = df['target_list'].apply(lambda t: len(set(t) & dili_genes))
df['n_cyp_hits'] = df['target_list'].apply(lambda t: len(set(t) & CYP_GENES))
df['n_tf_hits'] = df['target_list'].apply(lambda t: len(set(t) & TF_GENES))
df['n_transporter_hits'] = df['target_list'].apply(lambda t: len(set(t) & TRANSPORTERS))
df['n_network_hits'] = df['target_list'].apply(lambda t: len(set(t) & NETWORK_GENES))

# Hub fraction
df['network_fraction'] = df['n_network_hits'] / df['k'].clip(lower=1)
df['dili_fraction'] = df['n_dili_hits'] / df['k'].clip(lower=1)

# Block C: Regime indicators (for ML to learn routing)
df['has_network_targets'] = (df['n_network_hits'] >= 1).astype(int)
df['has_transporter_targets'] = (df['n_transporter_hits'] >= 1).astype(int)
df['high_logp'] = (df['logp'] > 3).astype(int) if 'logp' in df.columns else 0

# ECNP features
df['ecnp_z'] = df['Z']
df['ecnp_eligible'] = ((df['n_network_hits'] >= 1) | (df['n_dili_hits'] >= 1)).astype(int)

# =============================================================================
# ML REGIME CLASSIFIER (Layer 1)
# =============================================================================

print("\nBuilding ML regime classifier...")

# Features for regime prediction
feature_cols = [
    # Chemistry (Block A)
    'logp', 'mw', 'hbd', 'hba', 'tpsa',
    # Topology (Block B)
    'k', 'n_dili_hits', 'n_network_hits', 'n_transporter_hits',
    'network_fraction', 'dili_fraction',
    # Regime indicators (Block C)
    'has_network_targets', 'has_transporter_targets', 'high_logp'
]

# Filter to valid rows
valid_mask = df[feature_cols].notna().all(axis=1)
df_valid = df[valid_mask].copy()
print(f"Valid samples: {len(df_valid)}")

X = df_valid[feature_cols].values
y = df_valid['is_dili'].values

# Scale features
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# =============================================================================
# MODEL TRAINING WITH ABLATION
# =============================================================================

print("\n" + "="*60)
print("ABLATION STUDY")
print("="*60)

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
results = {}

# Model 1: ML only (no ECNP features)
print("\n1. ML WITHOUT ECNP")
ml_features = ['logp', 'mw', 'hbd', 'hba', 'tpsa', 'k', 
               'has_network_targets', 'has_transporter_targets', 'high_logp']
X_ml = df_valid[ml_features].values
X_ml_scaled = StandardScaler().fit_transform(X_ml)

if XGBOOST_AVAILABLE:
    model_ml = XGBClassifier(n_estimators=100, max_depth=3, random_state=42, eval_metric='logloss')
else:
    model_ml = LogisticRegression(max_iter=1000, random_state=42)

y_pred_ml = cross_val_predict(model_ml, X_ml_scaled, y, cv=cv, method='predict_proba')[:, 1]
auc_ml = roc_auc_score(y, y_pred_ml)
print(f"   AUC = {auc_ml:.3f}")
results['ML without ECNP'] = auc_ml

# Model 2: ECNP only (network expert alone)
print("\n2. ECNP ALONE")
auc_ecnp = roc_auc_score(y, df_valid['ecnp_z'].values)
print(f"   AUC = {auc_ecnp:.3f}")
results['ECNP alone'] = auc_ecnp

# Model 3: Full model (ML + ECNP integration)
print("\n3. FULL MODEL (ML + ECNP)")
full_features = ml_features + ['ecnp_z', 'n_dili_hits', 'n_network_hits', 'dili_fraction']
X_full = df_valid[full_features].values
X_full_scaled = StandardScaler().fit_transform(X_full)

if XGBOOST_AVAILABLE:
    model_full = XGBClassifier(n_estimators=100, max_depth=3, random_state=42, eval_metric='logloss')
else:
    model_full = LogisticRegression(max_iter=1000, random_state=42)

y_pred_full = cross_val_predict(model_full, X_full_scaled, y, cv=cv, method='predict_proba')[:, 1]
auc_full = roc_auc_score(y, y_pred_full)
print(f"   AUC = {auc_full:.3f}")
results['Full model (ML + ECNP)'] = auc_full

# Model 4: LogP alone (exposure baseline)
print("\n4. LOGP ALONE (Exposure baseline)")
auc_logp = roc_auc_score(y, df_valid['logp'].values)
print(f"   AUC = {auc_logp:.3f}")
results['LogP alone'] = auc_logp

# =============================================================================
# FEATURE IMPORTANCE
# =============================================================================

print("\n" + "="*60)
print("FEATURE IMPORTANCE (Full Model)")
print("="*60)

model_full.fit(X_full_scaled, y)

if XGBOOST_AVAILABLE:
    importances = model_full.feature_importances_
else:
    importances = np.abs(model_full.coef_[0])

feat_imp = sorted(zip(full_features, importances), key=lambda x: x[1], reverse=True)
for name, imp in feat_imp[:10]:
    print(f"  {name:20s}: {imp:.3f}")

# =============================================================================
# REGIME-SPECIFIC PERFORMANCE
# =============================================================================

print("\n" + "="*60)
print("REGIME-SPECIFIC AUC")
print("="*60)

df_valid['y_pred_full'] = y_pred_full

# ECNP-eligible subset
ecnp_eligible = df_valid[df_valid['ecnp_eligible'] == 1]
if len(ecnp_eligible) >= 20:
    y_t = ecnp_eligible['is_dili'].values
    y_p = ecnp_eligible['y_pred_full'].values
    n_pos = y_t.sum()
    n_neg = len(y_t) - n_pos
    if n_pos >= 3 and n_neg >= 3:
        auc_ecnp_regime = roc_auc_score(y_t, y_p)
        print(f"\nECNP-eligible (n={len(ecnp_eligible)}, {n_pos}:{n_neg})")
        print(f"  Full model AUC = {auc_ecnp_regime:.3f}")
        print(f"  ECNP alone AUC = {roc_auc_score(y_t, ecnp_eligible['ecnp_z']):.3f}")

# By target type
for regime, condition in [
    ('Network targets', df_valid['has_network_targets'] == 1),
    ('Transporter targets', df_valid['has_transporter_targets'] == 1),
    ('High LogP', df_valid['high_logp'] == 1),
]:
    subset = df_valid[condition]
    if len(subset) >= 15:
        y_t = subset['is_dili'].values
        y_p = subset['y_pred_full'].values
        n_pos = y_t.sum()
        n_neg = len(y_t) - n_pos
        if n_pos >= 3 and n_neg >= 3:
            auc = roc_auc_score(y_t, y_p)
            print(f"\n{regime} (n={len(subset)}, {n_pos}:{n_neg})")
            print(f"  AUC = {auc:.3f}")

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "="*60)
print("SUMMARY: ABLATION RESULTS")
print("="*60)

for name, auc in sorted(results.items(), key=lambda x: x[1], reverse=True):
    print(f"  {name:25s}: AUC = {auc:.3f}")

# ECNP contribution
ecnp_contrib = results['Full model (ML + ECNP)'] - results['ML without ECNP']
print(f"\nECNP contribution: +{ecnp_contrib:.3f}")

if results['Full model (ML + ECNP)'] > results['ML without ECNP'] and results['Full model (ML + ECNP)'] > results['ECNP alone']:
    print("\nVERDICT: Full model beats both ablations - integration justified")
elif results['ML without ECNP'] > results['Full model (ML + ECNP)']:
    print("\nVERDICT: ML alone better - ECNP not adding value")
else:
    print("\nVERDICT: Mixed results - further investigation needed")

# Save predictions
output_file = ROOT / 'research' / 'ecnp-generalization' / 'results' / 'ml_pipeline_results.csv'
df_valid.to_csv(output_file, index=False)
print(f"\nSaved: {output_file}")

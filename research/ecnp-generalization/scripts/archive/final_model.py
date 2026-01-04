"""
Final DILI Prediction Model with Clinical Risk Flags
=====================================================

Combines:
1. ML model (0.88 AUC) for base prediction
2. HLA genetic risk flags as clinical overlay
3. Regime classification for interpretability
"""
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score
import warnings
warnings.filterwarnings('ignore')

ROOT = Path(r'v:\new\h-perforatum-network-tox')

# =============================================================================
# LOAD DATA
# =============================================================================

print("Loading data...")
df = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'results' / 'ecfp_model_results.csv')
hla = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'data' / 'curated' / 'hla_dili_associations.csv')

print(f"Compounds: {len(df)}")

# =============================================================================
# FINAL RISK SCORE
# =============================================================================

# Use combined ECFP + tabular prediction
if 'y_pred_combined' in df.columns:
    df['risk_score'] = df['y_pred_combined']
elif 'y_pred_ecfp' in df.columns:
    df['risk_score'] = df['y_pred_ecfp']

# =============================================================================
# HLA CLINICAL FLAGS
# =============================================================================

print("\nAdding clinical risk flags...")

# Create HLA lookup
hla_lookup = {}
for _, row in hla.iterrows():
    drug = row['drug_name'].lower().strip()
    if drug not in hla_lookup:
        hla_lookup[drug] = []
    hla_lookup[drug].append({
        'allele': row['hla_allele'],
        'odds_ratio': row['odds_ratio'],
        'dili_type': row['dili_type'],
        'source': row['source']
    })

# Add HLA flags
df['drug_name_lower'] = df['drug_name'].str.lower().str.strip()
df['hla_risk'] = df['drug_name_lower'].apply(lambda x: x in hla_lookup)

def get_hla_warning(drug_name):
    drug = drug_name.lower().strip()
    if drug not in hla_lookup:
        return ""
    alleles = [a['allele'] for a in hla_lookup[drug]]
    return f"HLA screening recommended: {', '.join(alleles)}"

df['hla_warning'] = df['drug_name'].apply(get_hla_warning)

# =============================================================================
# REGIME CLASSIFICATION
# =============================================================================

print("Classifying regimes...")

# Network regime
df['regime'] = 'unknown'
if 'ecnp_eligible' in df.columns:
    df.loc[df['ecnp_eligible'] == 1, 'regime'] = 'network'

# High LogP regime
if 'logp' in df.columns:
    high_logp = (df['logp'] >= 3) & (df['regime'] == 'unknown')
    df.loc[high_logp, 'regime'] = 'intrinsic'

# HLA regime override for clinical flagging
df.loc[df['hla_risk'], 'regime'] = df.loc[df['hla_risk'], 'regime'] + '+HLA'

print(f"Regime distribution:")
print(df['regime'].value_counts())

# =============================================================================
# RISK TIER CLASSIFICATION
# =============================================================================

print("\nClassifying risk tiers...")

def classify_risk(row):
    score = row['risk_score']
    hla = row['hla_risk']
    
    if score >= 0.8:
        return "HIGH_RISK"
    elif score >= 0.5:
        if hla:
            return "HIGH_RISK"  # HLA elevates moderate to high
        return "MODERATE_RISK"
    else:
        if hla:
            return "MODERATE_RISK"  # HLA elevates low to moderate
        return "LOW_RISK"

df['risk_tier'] = df.apply(classify_risk, axis=1)

print(f"Risk tier distribution:")
print(df['risk_tier'].value_counts())

# =============================================================================
# VALIDATION
# =============================================================================

print("\n" + "="*60)
print("FINAL MODEL VALIDATION")
print("="*60)

# AUC on original score (unchanged)
auc = roc_auc_score(df['is_dili'], df['risk_score'])
print(f"AUC: {auc:.3f}")

# Confusion matrix by risk tier
print("\nDILI rate by risk tier:")
for tier in ['HIGH_RISK', 'MODERATE_RISK', 'LOW_RISK']:
    subset = df[df['risk_tier'] == tier]
    if len(subset) > 0:
        dili_rate = subset['is_dili'].mean() * 100
        print(f"  {tier}: {subset['is_dili'].sum()}/{len(subset)} ({dili_rate:.0f}%)")

# HLA rescue of false negatives
df['y_pred_binary'] = (df['risk_score'] > 0.5).astype(int)
df['fn'] = (df['y_pred_binary'] == 0) & (df['is_dili'] == 1)

fn_total = df['fn'].sum()
fn_with_hla = df[df['fn'] & df['hla_risk']].shape[0]
print(f"\nFalse negatives: {fn_total}")
print(f"  With HLA warning: {fn_with_hla} (now flagged)")
print(f"  Without HLA info: {fn_total - fn_with_hla} (irreducible)")

# =============================================================================
# EXAMPLE OUTPUT
# =============================================================================

print("\n" + "="*60)
print("EXAMPLE PREDICTIONS")
print("="*60)

# Show some examples
examples = df.nlargest(5, 'risk_score')[['drug_name', 'risk_score', 'risk_tier', 'regime', 'hla_warning']]
print("\nHIGH RISK examples:")
for _, row in examples.iterrows():
    warning = f" [{row['hla_warning']}]" if row['hla_warning'] else ""
    print(f"  {row['drug_name']}: {row['risk_score']:.2f} ({row['risk_tier']}, {row['regime']}){warning}")

examples = df[df['hla_risk']][['drug_name', 'risk_score', 'risk_tier', 'regime', 'hla_warning']]
print("\nHLA-FLAGGED drugs:")
for _, row in examples.iterrows():
    print(f"  {row['drug_name']}: {row['risk_score']:.2f} ({row['risk_tier']}) -> {row['hla_warning']}")

# =============================================================================
# SAVE FINAL MODEL
# =============================================================================

output_cols = [
    'drug_name', 'drugbank_id', 'smiles', 'is_dili',
    'risk_score', 'risk_tier', 'regime',
    'hla_risk', 'hla_warning',
    'logp', 'mw', 'ecnp_z', 'k'
]
output_cols = [c for c in output_cols if c in df.columns]

output = df[output_cols].copy()
output = output.sort_values('risk_score', ascending=False)

output_path = ROOT / 'research' / 'ecnp-generalization' / 'results' / 'final_dili_predictions.csv'
output.to_csv(output_path, index=False)
print(f"\nSaved: {output_path}")

print("\n" + "="*60)
print("FINAL MODEL SUMMARY")
print("="*60)
print(f"""
Model: ECFP4 + Tabular Features + ECNP Network Score
AUC: {auc:.3f}

Risk Tiers:
  HIGH_RISK: score >= 0.8 OR (score >= 0.5 AND HLA risk)
  MODERATE_RISK: score 0.5-0.8 OR (score < 0.5 AND HLA risk)
  LOW_RISK: score < 0.5 AND no HLA risk

Clinical Flags:
  HLA genetic risk warnings for {df['hla_risk'].sum()} drugs
  {fn_with_hla} false negatives rescued by HLA overlay

Regime Classification:
  network: ECNP-eligible (targets DILI genes)
  intrinsic: High lipophilicity (LogP >= 3)
  unknown: Neither (may be idiosyncratic)
  +HLA: Has known genetic risk factor
""")

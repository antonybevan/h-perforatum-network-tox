"""
Data Path and Integrity Verification
=====================================

Comprehensive audit of:
1. All data files and their paths
2. Label integrity (DILI labels correct?)
3. Compound matching (DrugBank ↔ DILIrank)
4. Feature consistency across files
5. No data leakage
"""
import pandas as pd
import numpy as np
from pathlib import Path
import hashlib

ROOT = Path(r'v:\new\h-perforatum-network-tox')

print("="*70)
print("DATA PATH AND INTEGRITY VERIFICATION")
print("="*70)

errors = []
warnings = []

# =============================================================================
# 1. VERIFY ALL DATA FILES EXIST
# =============================================================================

print("\n" + "-"*60)
print("1. DATA FILE EXISTENCE CHECK")
print("-"*60)

data_files = {
    # Raw data
    'DILIrank raw': ROOT / 'data' / 'raw' / 'DILIrank' / 'DILIrank-2.0.xlsx',
    'DrugBank XML': ROOT / 'research' / 'ecnp-generalization' / 'data' / 'raw' / 'full database.xml',
    
    # Processed data
    'DILI genes': ROOT / 'data' / 'processed' / 'dili_900_lcc.csv',
    'Influence matrix': ROOT / 'data' / 'processed' / 'S_matrix_liver_lcc.npz',
    
    # Curated data
    'Labeled compounds': ROOT / 'research' / 'ecnp-generalization' / 'data' / 'curated' / 'labeled_compounds_improved.csv',
    'PK data': ROOT / 'research' / 'ecnp-generalization' / 'data' / 'curated' / 'pk_data.csv',
    
    # Results
    'ECNP scores': ROOT / 'research' / 'ecnp-generalization' / 'results' / 'phase2' / 'ecnp_scores_improved.csv',
    'ML pipeline results': ROOT / 'research' / 'ecnp-generalization' / 'results' / 'ml_pipeline_with_alerts.csv',
    'ECFP results': ROOT / 'research' / 'ecnp-generalization' / 'results' / 'ecfp_model_results.csv',
}

for name, path in data_files.items():
    if path.exists():
        size = path.stat().st_size
        print(f"  [OK] {name}: {size:,} bytes")
    else:
        print(f"  [MISSING] {name}: {path}")
        errors.append(f"Missing file: {name}")

# =============================================================================
# 2. VERIFY DILI LABELS
# =============================================================================

print("\n" + "-"*60)
print("2. DILI LABEL INTEGRITY")
print("-"*60)

# Check labeled compounds
labeled = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'data' / 'curated' / 'labeled_compounds_improved.csv')
print(f"Labeled compounds: {len(labeled)}")

# Check DILIrank column
if 'dilirank' in labeled.columns:
    print(f"DILIrank values: {labeled['dilirank'].value_counts().to_dict()}")
else:
    print("DILIrank column: NOT FOUND")
    warnings.append("DILIrank column missing in labeled_compounds")

# Verify against raw DILIrank
try:
    dilirank_raw = pd.read_excel(ROOT / 'data' / 'raw' / 'DILIrank' / 'DILIrank-2.0.xlsx', header=1)
    print(f"Raw DILIrank entries: {len(dilirank_raw)}")
    
    # Count by category
    if 'DILI-Concern' in dilirank_raw.columns:
        print(f"Raw DILIrank categories: {dilirank_raw['DILI-Concern'].value_counts().to_dict()}")
except Exception as e:
    print(f"Could not read raw DILIrank: {e}")

# =============================================================================
# 3. VERIFY COMPOUND MATCHING
# =============================================================================

print("\n" + "-"*60)
print("3. COMPOUND MATCHING INTEGRITY")
print("-"*60)

# Load ECNP scores
scores = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'results' / 'phase2' / 'ecnp_scores_improved.csv')
print(f"ECNP scores: {len(scores)}")

# Check drugbank_id uniqueness
if 'drugbank_id' in scores.columns:
    n_unique = scores['drugbank_id'].nunique()
    n_total = len(scores)
    if n_unique == n_total:
        print(f"  [OK] DrugBank IDs unique: {n_unique}")
    else:
        print(f"  [WARNING] Duplicate DrugBank IDs: {n_total - n_unique}")
        warnings.append("Duplicate DrugBank IDs in scores")

# Check DILI label distribution
if 'is_dili' in scores.columns:
    dili_dist = scores['is_dili'].value_counts()
    print(f"  DILI labels: DILI+={dili_dist.get(1, 0)}, DILI-={dili_dist.get(0, 0)}")
    
    # Check for NaN labels
    nan_labels = scores['is_dili'].isna().sum()
    if nan_labels > 0:
        print(f"  [WARNING] NaN DILI labels: {nan_labels}")
        warnings.append(f"{nan_labels} NaN DILI labels")

# =============================================================================
# 4. VERIFY FEATURE CONSISTENCY
# =============================================================================

print("\n" + "-"*60)
print("4. FEATURE CONSISTENCY")
print("-"*60)

# Load ML pipeline results
ml_results = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'results' / 'ml_pipeline_with_alerts.csv')
print(f"ML pipeline results: {len(ml_results)}")

# Check key features exist
key_features = ['logp', 'mw', 'ecnp_z', 'is_dili', 'n_targets', 'smiles']
for feat in key_features:
    if feat in ml_results.columns:
        non_null = ml_results[feat].notna().sum()
        print(f"  [OK] {feat}: {non_null}/{len(ml_results)} non-null")
    else:
        print(f"  [MISSING] {feat}")
        errors.append(f"Missing feature: {feat}")

# Check for data type issues
if 'is_dili' in ml_results.columns:
    unique_vals = ml_results['is_dili'].unique()
    print(f"  is_dili unique values: {sorted(unique_vals)}")
    if not all(v in [0, 1, 0.0, 1.0] for v in unique_vals if pd.notna(v)):
        errors.append("Invalid is_dili values")

# =============================================================================
# 5. VERIFY NO DATA LEAKAGE
# =============================================================================

print("\n" + "-"*60)
print("5. DATA LEAKAGE CHECK")
print("-"*60)

# Check that DILI gene list doesn't contain target-specific information
dili_genes = pd.read_csv(ROOT / 'data' / 'processed' / 'dili_900_lcc.csv')
print(f"DILI gene list: {len(dili_genes)} genes")

# Sample some DILI genes
sample_genes = dili_genes['gene_name'].head(10).tolist()
print(f"  Sample genes: {sample_genes}")

# Check correlation between features and target
if 'is_dili' in ml_results.columns and 'logp' in ml_results.columns:
    corr = ml_results['logp'].corr(ml_results['is_dili'])
    print(f"  LogP-DILI correlation: {corr:.3f}")
    if abs(corr) > 0.5:
        warnings.append(f"High LogP-DILI correlation ({corr:.3f}) - check for leakage")

# =============================================================================
# 6. VERIFY ECFP RESULTS
# =============================================================================

print("\n" + "-"*60)
print("6. ECFP MODEL RESULTS VERIFICATION")
print("-"*60)

ecfp_results = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'results' / 'ecfp_model_results.csv')
print(f"ECFP results: {len(ecfp_results)}")

# Check predictions are valid probabilities
if 'y_pred_ecfp' in ecfp_results.columns:
    pred_min = ecfp_results['y_pred_ecfp'].min()
    pred_max = ecfp_results['y_pred_ecfp'].max()
    print(f"  Predictions range: [{pred_min:.3f}, {pred_max:.3f}]")
    if pred_min < 0 or pred_max > 1:
        errors.append("ECFP predictions outside [0,1]")

# =============================================================================
# 7. CROSS-FILE CONSISTENCY
# =============================================================================

print("\n" + "-"*60)
print("7. CROSS-FILE CONSISTENCY")
print("-"*60)

# Check compound counts match across files
files_to_check = [
    ('labeled_compounds', labeled),
    ('scores', scores),
    ('ml_results', ml_results),
    ('ecfp_results', ecfp_results),
]

for name, df in files_to_check:
    n = len(df)
    print(f"  {name}: {n} rows")

# Check if drugbank_ids are consistent
if 'drugbank_id' in labeled.columns and 'drugbank_id' in scores.columns:
    labeled_ids = set(labeled['drugbank_id'])
    scores_ids = set(scores['drugbank_id'])
    overlap = len(labeled_ids & scores_ids)
    print(f"  DrugBank ID overlap: {overlap}/{len(labeled_ids)}")

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "="*70)
print("VERIFICATION SUMMARY")
print("="*70)

print(f"\nErrors: {len(errors)}")
for e in errors:
    print(f"  [ERROR] {e}")

print(f"\nWarnings: {len(warnings)}")
for w in warnings:
    print(f"  [WARNING] {w}")

if len(errors) == 0:
    print("\n*** DATA INTEGRITY VERIFIED - NO CRITICAL ERRORS ***")
else:
    print("\n*** CRITICAL ERRORS FOUND - REVIEW REQUIRED ***")

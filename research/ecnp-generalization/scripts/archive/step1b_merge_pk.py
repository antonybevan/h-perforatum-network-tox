"""
Merge PK data and build full feature matrix for pipeline
"""
import pandas as pd
import numpy as np
from pathlib import Path

ROOT = Path(r'v:\new\h-perforatum-network-tox')

# Load data
pk = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'data' / 'curated' / 'pk_data.csv')
compounds = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'data' / 'curated' / 'labeled_compounds_improved.csv')
scores = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'results' / 'phase2' / 'ecnp_scores_improved.csv')

print(f'PK data rows: {len(pk)}')
print(f'Target compounds: {len(compounds)}')
print(f'Scores: {len(scores)}')

# Match PK by drug name (case insensitive)
pk['drug_name_lower'] = pk['drug_name'].str.lower().str.strip()
compounds['drug_name_lower'] = compounds['drug_name'].str.lower().str.strip()

# Merge compounds with PK
merged = compounds.merge(pk, on='drug_name_lower', how='left', suffixes=('', '_pk'))
print(f'\nPK matched: {merged["half_life_bin"].notna().sum()}/{len(merged)}')

# Merge with scores
merged = merged.merge(scores[['drugbank_id', 'Z', 'is_dili', 'I_T', 'n_targets']], on='drugbank_id', how='inner')
print(f'With scores: {len(merged)}')

# Create exposure features
merged['half_life_bin'] = merged['half_life_bin'].fillna('unknown')
merged['route'] = merged['route'].fillna('unknown')
merged['high_dose'] = merged['high_dose'].fillna(0)

# Summary
print(f'\n{"="*60}')
print('FEATURE AVAILABILITY SUMMARY')
print('='*60)

print(f'\nHalf-life bins:')
print(merged['half_life_bin'].value_counts())

print(f'\nHigh dose:')
print(merged['high_dose'].value_counts())

print(f'\nRoute:')
print(merged['route'].value_counts())

# Save
output_file = ROOT / 'research' / 'ecnp-generalization' / 'data' / 'curated' / 'compounds_full_features.csv'
merged.to_csv(output_file, index=False)
print(f'\nSaved: {output_file}')
print(f'Total compounds for modeling: {len(merged)}')

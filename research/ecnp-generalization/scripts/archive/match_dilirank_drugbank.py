"""
Match DILIrank to DrugBank
==========================

Match 982 non-ambiguous DILIrank drugs to our DrugBank dataset.
Use normalized name matching with synonym handling.
"""
import pandas as pd
import numpy as np
from pathlib import Path
import re
import warnings
warnings.filterwarnings('ignore')

ROOT = Path(r'v:\new\h-perforatum-network-tox')

print("="*70)
print("MATCHING DILIRANK TO DRUGBANK")
print("="*70)

# =============================================================================
# LOAD DILIRANK
# =============================================================================

print("\n--- Loading DILIrank ---")
dilirank = pd.read_excel(ROOT / 'data' / 'external' / 'DILIrank_2.0.xlsx', header=1)
print(f"Total: {len(dilirank)}")

# Binary labels, remove ambiguous
dilirank['is_dili'] = dilirank['vDILI-Concern'].apply(
    lambda x: 1 if 'Most' in str(x) or 'Less' in str(x) else 0
)
dilirank = dilirank[~dilirank['vDILI-Concern'].str.contains('Ambiguous', na=False)].copy()
print(f"Non-ambiguous: {len(dilirank)}")

# =============================================================================
# LOAD OUR DRUGBANK DATA
# =============================================================================

print("\n--- Loading DrugBank Data ---")

# Find our compounds file with SMILES
compounds_file = ROOT / 'data' / 'processed' / 'hp_compound_targets.csv'
if compounds_file.exists():
    drugbank = pd.read_csv(compounds_file)
else:
    # Try other locations
    for path in [
        ROOT / 'research' / 'ecnp-generalization' / 'results' / 'ecfp_model_results.csv',
        ROOT / 'data' / 'processed' / 'compounds_with_smiles.csv'
    ]:
        if path.exists():
            drugbank = pd.read_csv(path)
            break
    else:
        print("No DrugBank data found!")
        exit()

print(f"DrugBank compounds: {len(drugbank)}")
print(f"Columns: {list(drugbank.columns)}")

# =============================================================================
# NAME NORMALIZATION
# =============================================================================

DRUG_SYNONYMS = {
    'ciclosporin': 'cyclosporine', 'ciclosporine': 'cyclosporine',
    'cyclosporin': 'cyclosporine', 'cyclosporin a': 'cyclosporine',
    'paracetamol': 'acetaminophen', 'rifampin': 'rifampicin',
    'furosemide': 'frusemide', 'sulfasalazine': 'sulphasalazine',
    'sulfamethoxazole': 'sulphamethoxazole',
    'epinephrine': 'adrenaline', 'norepinephrine': 'noradrenaline',
    'aluminum': 'aluminium', 'estrogen': 'oestrogen',
    'levodopa': 'l-dopa', 'thyroxine': 'levothyroxine',
}

def normalize_name(name):
    """Normalize drug name for matching."""
    if pd.isna(name):
        return ''
    name = str(name).lower().strip()
    
    # Remove common suffixes
    suffixes = [' hydrochloride', ' hcl', ' sodium', ' potassium', ' sulfate', 
                ' acetate', ' mesylate', ' maleate', ' tartrate', ' citrate',
                ' phosphate', ' chloride', ' bromide', ' succinate', ' fumarate',
                ' dihydrate', ' monohydrate', ' trihydrate']
    for suffix in suffixes:
        name = name.replace(suffix, '')
    
    # Apply synonyms
    if name in DRUG_SYNONYMS:
        name = DRUG_SYNONYMS[name]
    
    return name.strip()

# Normalize names in both datasets
dilirank['norm_name'] = dilirank['CompoundName'].apply(normalize_name)

# Find drug name column in DrugBank data
name_col = None
for col in ['drug_name', 'name', 'DrugName', 'compound_name']:
    if col in drugbank.columns:
        name_col = col
        break

if name_col:
    drugbank['norm_name'] = drugbank[name_col].apply(normalize_name)
else:
    print("No drug name column found in DrugBank data")
    exit()

# =============================================================================
# MATCHING
# =============================================================================

print("\n--- Matching Names ---")

# Drop duplicates before creating lookup
drugbank_unique = drugbank.drop_duplicates(subset='norm_name', keep='first')
print(f"Unique DrugBank names: {len(drugbank_unique)}")

# Create lookup from DrugBank
drugbank_lookup = drugbank_unique.set_index('norm_name').to_dict('index')
drugbank_names = set(drugbank_unique['norm_name'])

# Match
matches = []
not_matched = []

for idx, row in dilirank.iterrows():
    norm = row['norm_name']
    
    if norm in drugbank_names:
        db_row = drugbank_lookup[norm]
        matches.append({
            'dilirank_name': row['CompoundName'],
            'drugbank_name': db_row.get(name_col, norm),
            'norm_name': norm,
            'is_dili': row['is_dili'],
            'smiles': db_row.get('smiles'),
            'drugbank_id': db_row.get('drugbank_id', db_row.get('db_id'))
        })
    else:
        not_matched.append({
            'dilirank_name': row['CompoundName'],
            'norm_name': norm,
            'is_dili': row['is_dili']
        })

print(f"Matched: {len(matches)}")
print(f"Not matched: {len(not_matched)}")

df_matched = pd.DataFrame(matches)
df_not_matched = pd.DataFrame(not_matched)

# =============================================================================
# EXPANDED DATASET
# =============================================================================

print("\n" + "="*70)
print("EXPANDED DATASET")
print("="*70)

if len(df_matched) > 0:
    # Filter to those with valid SMILES
    df_matched = df_matched[df_matched['smiles'].notna()].copy()
    print(f"With valid SMILES: {len(df_matched)}")
    
    # Label distribution
    print(f"\nLabel distribution:")
    print(f"  DILI+: {df_matched['is_dili'].sum()}")
    print(f"  DILI-: {(df_matched['is_dili'] == 0).sum()}")
    
    # Compare to original
    print(f"\n{'Dataset':<25} {'Drugs':<10} {'Coverage of DILIrank':<20}")
    print("-"*55)
    print(f"{'Original subset':<25} {'202':<10} {'15%':<20}")
    print(f"{'Expanded (matched)':<25} {len(df_matched):<10} {len(df_matched)/982*100:.1f}%")
    
    # Save
    output = ROOT / 'research' / 'ecnp-generalization' / 'results' / 'dilirank_drugbank_matched.csv'
    df_matched.to_csv(output, index=False)
    print(f"\nSaved matched: {output}")
    
    df_not_matched.to_csv(
        ROOT / 'research' / 'ecnp-generalization' / 'results' / 'dilirank_not_matched.csv',
        index=False
    )
    print(f"Saved not matched: dilirank_not_matched.csv")

print("\n" + "="*70)
print("SUMMARY")
print("="*70)
print(f"""
DILIrank total: 1336
Non-ambiguous: 982
Matched to DrugBank: {len(df_matched)}
Coverage: {len(df_matched)/982*100:.1f}%

Not matched: {len(df_not_matched)} (would need PubChem lookup)
""")

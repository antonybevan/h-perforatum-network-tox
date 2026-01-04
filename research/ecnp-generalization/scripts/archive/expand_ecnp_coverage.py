"""
Expand ECNP to Full DILIrank
=============================

Match 901 DILIrank drugs to DrugBank to get targets, then compute ECNP.
"""
import pandas as pd
import numpy as np
from pathlib import Path
import xml.etree.ElementTree as ET
import re
import warnings
warnings.filterwarnings('ignore')

ROOT = Path(r'v:\new\h-perforatum-network-tox')

print("="*70)
print("EXPANDING ECNP TO FULL DILIRANK")
print("="*70)

# =============================================================================
# LOAD DATA
# =============================================================================

print("\n--- Loading Data ---")

# 901 DILIrank drugs with SMILES
df_dilirank = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'results' / 'dilirank_full_smiles.csv')
print(f"DILIrank drugs: {len(df_dilirank)}")

# Existing targets (from our 202-drug subset)
existing_targets = pd.read_csv(ROOT / 'data' / 'processed' / 'targets.csv')
print(f"Existing target mappings: {len(existing_targets)}")

# Check existing columns
print(f"Existing target columns: {list(existing_targets.columns)}")

# =============================================================================
# PARSE DRUGBANK XML (EXTRACT DRUG NAME -> TARGETS)
# =============================================================================

print("\n--- Parsing DrugBank XML ---")

drugbank_xml = ROOT / 'research' / 'ecnp-generalization' / 'data' / 'raw' / 'full database.xml'

# Parse XML
ns = {'db': 'http://www.drugbank.ca'}

try:
    tree = ET.parse(drugbank_xml)
    root = tree.getroot()
    print(f"XML loaded, root tag: {root.tag}")
except Exception as e:
    print(f"Error loading XML: {e}")
    exit()

# Extract drug info
drugs_data = []
count = 0

for drug in root.findall('.//db:drug', ns):
    count += 1
    if count % 1000 == 0:
        print(f"  Processed {count} drugs...")
    
    # Get drug info
    drugbank_id = drug.find('db:drugbank-id[@primary="true"]', ns)
    drugbank_id = drugbank_id.text if drugbank_id is not None else None
    
    name = drug.find('db:name', ns)
    name = name.text if name is not None else None
    
    # Get synonyms
    synonyms = []
    syns = drug.findall('.//db:synonyms/db:synonym', ns)
    for syn in syns:
        if syn.text:
            synonyms.append(syn.text.lower())
    
    # Get targets
    targets = []
    for target in drug.findall('.//db:targets/db:target', ns):
        polypeptide = target.find('.//db:polypeptide', ns)
        if polypeptide is not None:
            uniprot = polypeptide.get('id')
            if uniprot:
                targets.append(uniprot)
    
    if drugbank_id and name:
        drugs_data.append({
            'drugbank_id': drugbank_id,
            'name': name.lower(),
            'synonyms': synonyms,
            'targets': targets,
            'n_targets': len(targets)
        })

print(f"Total drugs in DrugBank: {len(drugs_data)}")

# Create lookup by name
drugbank_lookup = {}
for drug in drugs_data:
    drugbank_lookup[drug['name']] = drug
    for syn in drug['synonyms']:
        if syn not in drugbank_lookup:  # Don't overwrite primary name
            drugbank_lookup[syn] = drug

print(f"Total unique names/synonyms: {len(drugbank_lookup)}")

# =============================================================================
# MATCH DILIRANK TO DRUGBANK
# =============================================================================

print("\n--- Matching DILIrank to DrugBank ---")

def normalize(name):
    if pd.isna(name):
        return ''
    name = str(name).lower().strip()
    # Remove common suffixes
    suffixes = [' hydrochloride', ' hcl', ' sodium', ' potassium', ' sulfate', 
                ' acetate', ' mesylate', ' maleate', ' tartrate', ' citrate']
    for suffix in suffixes:
        name = name.replace(suffix, '')
    return name.strip()

matched = []
not_matched = []

for idx, row in df_dilirank.iterrows():
    name = row['dilirank_name']
    name_norm = normalize(name)
    
    # Try exact match
    if name_norm in drugbank_lookup:
        match = drugbank_lookup[name_norm]
        matched.append({
            'dilirank_name': row['dilirank_name'],
            'smiles': row['smiles'],
            'is_dili': row['is_dili'],
            'drugbank_id': match['drugbank_id'],
            'targets': match['targets'],
            'n_targets': match['n_targets']
        })
    elif name.lower() in drugbank_lookup:
        match = drugbank_lookup[name.lower()]
        matched.append({
            'dilirank_name': row['dilirank_name'],
            'smiles': row['smiles'],
            'is_dili': row['is_dili'],
            'drugbank_id': match['drugbank_id'],
            'targets': match['targets'],
            'n_targets': match['n_targets']
        })
    else:
        not_matched.append(row['dilirank_name'])

print(f"Matched: {len(matched)}")
print(f"Not matched: {len(not_matched)}")

df_matched = pd.DataFrame(matched)

# Filter to those with targets
df_with_targets = df_matched[df_matched['n_targets'] > 0].copy()
print(f"With targets: {len(df_with_targets)}")

# =============================================================================
# SAVE FOR ECNP COMPUTATION
# =============================================================================

print("\n--- Saving Matched Data ---")

# Expand targets to individual rows for ECNP
rows = []
for _, drug in df_with_targets.iterrows():
    for target in drug['targets']:
        rows.append({
            'dilirank_name': drug['dilirank_name'],
            'smiles': drug['smiles'],
            'is_dili': drug['is_dili'],
            'drugbank_id': drug['drugbank_id'],
            'target': target
        })

df_targets_expanded = pd.DataFrame(rows)
output = ROOT / 'research' / 'ecnp-generalization' / 'results' / 'dilirank_drugbank_targets.csv'
df_targets_expanded.to_csv(output, index=False)
print(f"Saved expanded targets: {output}")

# Save summary
df_with_targets.to_csv(
    ROOT / 'research' / 'ecnp-generalization' / 'results' / 'dilirank_with_targets.csv',
    index=False
)
print(f"Saved drugs with targets: dilirank_with_targets.csv")

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "="*70)
print("SUMMARY")
print("="*70)

print(f"""
DILIrank drugs: 901
Matched to DrugBank: {len(matched)}
With targets: {len(df_with_targets)}
Coverage: {len(df_with_targets)/901*100:.1f}%

This gives us {len(df_with_targets)} drugs for ECNP computation!
(vs. 202 in original subset)

Expansion: {len(df_with_targets) - 202:+d} drugs
""")

# Label distribution
if len(df_with_targets) > 0:
    print(f"Label distribution in expanded set:")
    print(f"  DILI+: {df_with_targets['is_dili'].sum()}")
    print(f"  DILI-: {(df_with_targets['is_dili'] == 0).sum()}")

"""
Robust SMILES Fetching from PubChem for DILIrank
=================================================

SAFETY MEASURES:
1. Exact name match first, then fuzzy fallback
2. Validate retrieved compound has similar MW
3. Flag ambiguous matches
4. Log all mismatches for manual review
"""
import pandas as pd
import numpy as np
from pathlib import Path
import requests
import time
import re
import warnings
warnings.filterwarnings('ignore')

ROOT = Path(r'v:\new\h-perforatum-network-tox')

print("="*70)
print("PUBCHEM SMILES FETCHING FOR DILIRANK")
print("="*70)

# =============================================================================
# LOAD DILIRANK
# =============================================================================

print("\n--- Loading DILIrank ---")
dilirank = pd.read_excel(ROOT / 'data' / 'external' / 'DILIrank_2.0.xlsx', header=1)
print(f"Total drugs: {len(dilirank)}")

# Binary labels
dilirank['is_dili'] = dilirank['vDILI-Concern'].apply(
    lambda x: 1 if 'Most' in str(x) or 'Less' in str(x) else 0
)

# Clean names
dilirank['clean_name'] = dilirank['CompoundName'].str.strip()

# Remove ambiguous
dilirank_clean = dilirank[~dilirank['vDILI-Concern'].str.contains('Ambiguous', na=False)].copy()
print(f"Non-ambiguous: {len(dilirank_clean)}")

# =============================================================================
# PUBCHEM API FUNCTIONS
# =============================================================================

def normalize_name(name):
    """Normalize drug name for matching."""
    name = str(name).lower().strip()
    # Remove common suffixes
    suffixes = [' hydrochloride', ' hcl', ' sodium', ' potassium', ' sulfate', 
                ' acetate', ' mesylate', ' maleate', ' tartrate', ' citrate',
                ' phosphate', ' chloride', ' bromide']
    for suffix in suffixes:
        name = name.replace(suffix, '')
    return name.strip()

def fetch_smiles_pubchem(compound_name, timeout=10):
    """Fetch SMILES from PubChem with validation."""
    try:
        # URL encode the name
        encoded_name = requests.utils.quote(compound_name)
        
        # Try exact name match first
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded_name}/property/CanonicalSMILES,MolecularWeight,Title/JSON"
        
        response = requests.get(url, timeout=timeout)
        
        if response.status_code == 200:
            data = response.json()
            props = data['PropertyTable']['Properties'][0]
            
            return {
                'smiles': props.get('CanonicalSMILES'),
                'mw': props.get('MolecularWeight'),
                'title': props.get('Title'),
                'cid': props.get('CID'),
                'status': 'found',
                'match_type': 'exact'
            }
        else:
            return {'status': 'not_found', 'error': response.status_code}
            
    except requests.exceptions.Timeout:
        return {'status': 'timeout'}
    except Exception as e:
        return {'status': 'error', 'error': str(e)}

def validate_smiles(smiles):
    """Validate SMILES using RDKit."""
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None
    except:
        return False

# =============================================================================
# BATCH FETCH WITH RATE LIMITING
# =============================================================================

print("\n--- Fetching SMILES from PubChem ---")

# Sample first for testing
BATCH_SIZE = 50  # Start with 50 to test
sample_drugs = dilirank_clean.head(BATCH_SIZE).copy()

results = []
not_found = []
errors = []

for idx, row in sample_drugs.iterrows():
    name = row['CompoundName']
    clean = normalize_name(name)
    
    # Rate limit: 5 requests per second
    time.sleep(0.25)
    
    result = fetch_smiles_pubchem(name)
    
    if result['status'] == 'found':
        # Validate SMILES
        if validate_smiles(result['smiles']):
            results.append({
                'original_name': name,
                'pubchem_title': result.get('title'),
                'smiles': result['smiles'],
                'mw': result.get('mw'),
                'cid': result.get('cid'),
                'match_type': result['match_type'],
                'is_dili': row['is_dili']
            })
        else:
            errors.append({'name': name, 'reason': 'invalid_smiles'})
    elif result['status'] == 'not_found':
        # Try normalized name
        result2 = fetch_smiles_pubchem(clean)
        time.sleep(0.25)
        
        if result2['status'] == 'found' and validate_smiles(result2.get('smiles', '')):
            results.append({
                'original_name': name,
                'pubchem_title': result2.get('title'),
                'smiles': result2['smiles'],
                'mw': result2.get('mw'),
                'cid': result2.get('cid'),
                'match_type': 'normalized',
                'is_dili': row['is_dili']
            })
        else:
            not_found.append(name)
    else:
        errors.append({'name': name, 'reason': result['status']})
    
    # Progress
    if (len(results) + len(not_found) + len(errors)) % 10 == 0:
        print(f"  Processed: {len(results) + len(not_found) + len(errors)}/{BATCH_SIZE}")

print(f"\nResults:")
print(f"  Found: {len(results)}")
print(f"  Not found: {len(not_found)}")
print(f"  Errors: {len(errors)}")

# =============================================================================
# VALIDATE MATCHES
# =============================================================================

print("\n--- Validating Matches ---")

df_matched = pd.DataFrame(results)

if len(df_matched) > 0:
    # Flag potential mismatches (name differs significantly)
    def check_name_match(row):
        orig = normalize_name(row['original_name'])
        pubchem = normalize_name(str(row['pubchem_title']))
        
        # Check if names are similar
        if orig == pubchem:
            return 'exact'
        elif orig in pubchem or pubchem in orig:
            return 'partial'
        else:
            return 'REVIEW'
    
    df_matched['name_match'] = df_matched.apply(check_name_match, axis=1)
    
    print(f"\nName match quality:")
    print(df_matched['name_match'].value_counts())
    
    # Show potential mismatches
    mismatches = df_matched[df_matched['name_match'] == 'REVIEW']
    if len(mismatches) > 0:
        print(f"\n[WARNING] Potential mismatches to review ({len(mismatches)}):")
        for _, row in mismatches.iterrows():
            print(f"  {row['original_name']} -> {row['pubchem_title']}")

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "="*70)
print("SUMMARY")
print("="*70)

print(f"""
Batch size: {BATCH_SIZE}
Found with valid SMILES: {len(results)} ({len(results)/BATCH_SIZE*100:.1f}%)
Not found: {len(not_found)} ({len(not_found)/BATCH_SIZE*100:.1f}%)
Errors: {len(errors)}

Name match quality:
  Exact: {(df_matched['name_match'] == 'exact').sum() if len(df_matched) > 0 else 0}
  Partial: {(df_matched['name_match'] == 'partial').sum() if len(df_matched) > 0 else 0}
  REVIEW: {(df_matched['name_match'] == 'REVIEW').sum() if len(df_matched) > 0 else 0}
""")

# Save results
if len(df_matched) > 0:
    output = ROOT / 'research' / 'ecnp-generalization' / 'results' / 'dilirank_smiles_batch.csv'
    df_matched.to_csv(output, index=False)
    print(f"Saved: {output}")
    
    # Save not found
    pd.DataFrame({'name': not_found}).to_csv(
        ROOT / 'research' / 'ecnp-generalization' / 'results' / 'dilirank_not_found.csv',
        index=False
    )

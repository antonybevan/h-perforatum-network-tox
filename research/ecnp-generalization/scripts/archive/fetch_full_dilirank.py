"""
Full DILIrank SMILES Fetch
==========================

Using PubChem TXT endpoint which works reliably.
Fetch SMILES for all 982 non-ambiguous DILIrank drugs.
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
print("FETCHING SMILES FOR FULL DILIRANK")
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
# PUBCHEM FETCH FUNCTION
# =============================================================================

def fetch_smiles_pubchem(name, timeout=10):
    """Fetch SMILES from PubChem using TXT endpoint."""
    try:
        encoded = requests.utils.quote(str(name).strip())
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded}/property/CanonicalSMILES/TXT"
        r = requests.get(url, timeout=timeout)
        if r.ok and r.text.strip():
            smiles = r.text.strip().split('\n')[0]  # Take first line
            if smiles and len(smiles) > 3:  # Basic validation
                return {'status': 'ok', 'smiles': smiles}
        return {'status': 'not_found'}
    except requests.exceptions.Timeout:
        return {'status': 'timeout'}
    except Exception as e:
        return {'status': 'error', 'error': str(e)}

# =============================================================================
# BATCH FETCH
# =============================================================================

print("\n--- Fetching SMILES ---")

results = []
not_found = []
errors = []

for i, (idx, row) in enumerate(dilirank.iterrows()):
    name = row['CompoundName']
    
    # Rate limit: 5 req/sec
    time.sleep(0.25)
    
    result = fetch_smiles_pubchem(name)
    
    if result['status'] == 'ok':
        results.append({
            'dilirank_name': name,
            'smiles': result['smiles'],
            'is_dili': row['is_dili'],
            'dili_concern': row['vDILI-Concern'],
            'ltkb_id': row['LTKBID']
        })
    elif result['status'] == 'not_found':
        not_found.append(name)
    else:
        errors.append({'name': name, 'error': result.get('error', 'unknown')})
    
    # Progress
    if (i + 1) % 100 == 0:
        print(f"  Processed: {i+1}/{len(dilirank)} | Found: {len(results)} | Not found: {len(not_found)} | Errors: {len(errors)}")

print(f"\nFinal:")
print(f"  Found: {len(results)}")
print(f"  Not found: {len(not_found)}")
print(f"  Errors: {len(errors)}")

# =============================================================================
# SAVE RESULTS
# =============================================================================

df_results = pd.DataFrame(results)

if len(df_results) > 0:
    # SMILES validation with RDKit
    try:
        from rdkit import Chem
        valid_smiles = []
        for idx, row in df_results.iterrows():
            mol = Chem.MolFromSmiles(row['smiles'])
            valid_smiles.append(mol is not None)
        df_results['valid_smiles'] = valid_smiles
        df_results = df_results[df_results['valid_smiles']].drop('valid_smiles', axis=1)
        print(f"  With valid SMILES (RDKit): {len(df_results)}")
    except:
        print("  RDKit not available, skipping validation")
    
    # Save
    output = ROOT / 'research' / 'ecnp-generalization' / 'results' / 'dilirank_full_smiles.csv'
    df_results.to_csv(output, index=False)
    print(f"\nSaved: {output}")
    
    # Label distribution
    print(f"\nLabel distribution:")
    print(f"  DILI+: {df_results['is_dili'].sum()}")
    print(f"  DILI-: {(df_results['is_dili'] == 0).sum()}")

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "="*70)
print("SUMMARY")
print("="*70)

print(f"""
DILIrank non-ambiguous: 982
SMILES found: {len(df_results)}
Coverage: {len(df_results)/982*100:.1f}%

This gives us {len(df_results)} drugs for Tier 1 (chemistry-only) model!

Not found drugs saved to: dilirank_not_found_pubchem.csv
""")

# Save not found
pd.DataFrame({'name': not_found}).to_csv(
    ROOT / 'research' / 'ecnp-generalization' / 'results' / 'dilirank_not_found_pubchem.csv',
    index=False
)

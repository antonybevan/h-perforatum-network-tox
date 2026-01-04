"""
Alternative SMILES APIs Test
============================

1. PubChem - Use SDF endpoint instead of JSON
2. ChEMBL - Use molecule search API
3. UniChem - Cross-reference multiple databases
"""
import requests
import time
import warnings
warnings.filterwarnings('ignore')

print("="*70)
print("TESTING ALTERNATIVE SMILES APIs")
print("="*70)

drugs = ['Acetaminophen', 'Methotrexate', 'Ibuprofen', 'Amoxicillin', 
         'Cyclosporine', 'Rifampicin', 'Isoniazid', 'Valproic acid']

# =============================================================================
# METHOD 1: PubChem with different endpoint
# =============================================================================

print("\n--- Method 1: PubChem REST (different format) ---")

def pubchem_smiles_v2(name):
    """Try getting SMILES via PubChem CID lookup."""
    try:
        # First get CID
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{requests.utils.quote(name)}/cids/JSON"
        r = requests.get(url, timeout=10)
        if not r.ok:
            return None
        
        cid = r.json()['IdentifierList']['CID'][0]
        
        # Then get SMILES by CID
        url2 = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES/TXT"
        r2 = requests.get(url2, timeout=10)
        if r2.ok:
            return r2.text.strip()
        return None
    except Exception as e:
        return None

success_v2 = 0
for drug in drugs[:4]:
    smiles = pubchem_smiles_v2(drug)
    if smiles:
        print(f"  OK: {drug} -> {smiles[:40]}...")
        success_v2 += 1
    else:
        print(f"  FAIL: {drug}")
    time.sleep(0.3)
print(f"  Success: {success_v2}/4")

# =============================================================================
# METHOD 2: ChEMBL API
# =============================================================================

print("\n--- Method 2: ChEMBL API ---")

def chembl_smiles(name):
    """Get SMILES from ChEMBL by molecule name search."""
    try:
        # Search by synonym
        url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/search.json?q={requests.utils.quote(name)}"
        r = requests.get(url, timeout=15)
        if r.ok:
            data = r.json()
            if data.get('molecules') and len(data['molecules']) > 0:
                mol = data['molecules'][0]
                smiles = mol.get('molecule_structures', {}).get('canonical_smiles')
                return smiles
        return None
    except Exception as e:
        return None

success_chembl = 0
for drug in drugs[:4]:
    smiles = chembl_smiles(drug)
    if smiles:
        print(f"  OK: {drug} -> {smiles[:40]}...")
        success_chembl += 1
    else:
        print(f"  FAIL: {drug}")
    time.sleep(0.5)
print(f"  Success: {success_chembl}/4")

# =============================================================================
# METHOD 3: PubChem via CID and SDF
# =============================================================================

print("\n--- Method 3: PubChem TXT endpoint ---")

def pubchem_txt(name):
    """Get SMILES directly as text."""
    try:
        encoded = requests.utils.quote(name)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded}/property/CanonicalSMILES/TXT"
        r = requests.get(url, timeout=10)
        if r.ok:
            return r.text.strip()
        return None
    except:
        return None

success_txt = 0
for drug in drugs[:4]:
    smiles = pubchem_txt(drug)
    if smiles:
        print(f"  OK: {drug} -> {smiles[:40]}...")
        success_txt += 1
    else:
        print(f"  FAIL: {drug}")
    time.sleep(0.3)
print(f"  Success: {success_txt}/4")

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "="*70)
print("SUMMARY")
print("="*70)

print(f"""
Method 1 (PubChem CID+TXT): {success_v2}/4
Method 2 (ChEMBL):          {success_chembl}/4
Method 3 (PubChem TXT):     {success_txt}/4
""")

best = max([('PubChem CID+TXT', success_v2), ('ChEMBL', success_chembl), ('PubChem TXT', success_txt)], key=lambda x: x[1])
print(f"Best method: {best[0]} ({best[1]}/4 success)")

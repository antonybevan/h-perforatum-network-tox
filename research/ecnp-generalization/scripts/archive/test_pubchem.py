"""Quick PubChem test"""
import requests
import time

def fetch_smiles(name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{requests.utils.quote(name)}/property/IsomericSMILES,MolecularWeight/JSON"
    try:
        r = requests.get(url, timeout=10)
        if r.ok:
            data = r.json()
            props = data['PropertyTable']['Properties'][0]
            # Try multiple SMILES keys
            smiles = props.get('IsomericSMILES') or props.get('CanonicalSMILES')
            return {'status': 'ok', 'smiles': smiles, 'mw': props.get('MolecularWeight')}
        return {'status': 'not_found', 'code': r.status_code}
    except Exception as e:
        return {'status': 'error', 'error': str(e)}

# Test
drugs = ['Acetaminophen', 'Methotrexate', 'Abacavir', 'Ibuprofen', 'Amoxicillin', 
         'Cyclosporine', 'Rifampicin', 'Isoniazid', 'Valproic acid', 'Carbamazepine']

print("Testing PubChem API:")
print("-" * 60)
success = 0
for drug in drugs:
    result = fetch_smiles(drug)
    if result['status'] == 'ok' and result['smiles']:
        print(f"OK: {drug} -> {result['smiles'][:40]}...")
        success += 1
    else:
        print(f"FAIL: {drug} -> {result}")
    time.sleep(0.3)

print(f"\nSuccess: {success}/{len(drugs)} ({success/len(drugs)*100:.0f}%)")

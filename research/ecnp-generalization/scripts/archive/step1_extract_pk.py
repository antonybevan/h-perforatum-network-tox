"""
Step 1: Extract PK/Dose Data from DrugBank
==========================================

Parse pharmacokinetic data for compounds in our validation set:
- Half-life → bins (short/medium/long)
- Dosage → daily dose estimate
- Route → oral vs parenteral
"""
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np
import re
from pathlib import Path

ROOT = Path(r'v:\new\h-perforatum-network-tox')
DRUGBANK_XML = ROOT / 'research' / 'ecnp-generalization' / 'data' / 'raw' / 'full database.xml'
OUTPUT_DIR = ROOT / 'research' / 'ecnp-generalization' / 'data' / 'curated'

# Load our compounds to match
compounds = pd.read_csv(OUTPUT_DIR / 'labeled_compounds_improved.csv')
target_drugs = set(compounds['drug_name'].str.lower())
print(f"Target drugs to extract PK for: {len(target_drugs)}")

# Parse DrugBank
ns = '{http://www.drugbank.ca}'

def parse_half_life(text):
    """Parse half-life text to hours."""
    if not text:
        return np.nan
    text = text.lower()
    
    # Extract numbers
    numbers = re.findall(r'(\d+\.?\d*)\s*(hour|hr|minute|min|day|week)', text)
    if not numbers:
        return np.nan
    
    hours = 0
    for num, unit in numbers:
        val = float(num)
        if 'min' in unit:
            hours = val / 60
        elif 'day' in unit:
            hours = val * 24
        elif 'week' in unit:
            hours = val * 24 * 7
        else:
            hours = val
    
    return hours

def half_life_bin(hours):
    """Bin half-life."""
    if pd.isna(hours):
        return 'unknown'
    if hours < 4:
        return 'short'
    elif hours < 24:
        return 'medium'
    else:
        return 'long'

def parse_dose(dosages_elem, ns):
    """Parse dosages to extract daily dose estimate."""
    if dosages_elem is None:
        return np.nan, 'unknown'
    
    max_dose = 0
    route = 'unknown'
    
    for dosage in dosages_elem:
        strength_elem = dosage.find(f'{ns}strength')
        route_elem = dosage.find(f'{ns}route')
        
        if strength_elem is not None and strength_elem.text:
            # Parse strength (e.g., "100 mg", "500mg")
            match = re.search(r'(\d+\.?\d*)\s*(mg|g|mcg|ug)', strength_elem.text.lower())
            if match:
                val = float(match.group(1))
                unit = match.group(2)
                if unit == 'g':
                    val *= 1000  # to mg
                elif unit in ['mcg', 'ug']:
                    val /= 1000  # to mg
                max_dose = max(max_dose, val)
        
        if route_elem is not None and route_elem.text:
            r = route_elem.text.lower()
            if 'oral' in r:
                route = 'oral'
            elif 'intravenous' in r or 'iv' in r:
                route = 'iv'
    
    return max_dose if max_dose > 0 else np.nan, route

print("\nParsing DrugBank XML...")
pk_data = []
matched = 0

context = ET.iterparse(DRUGBANK_XML, events=('end',))
for event, elem in context:
    if elem.tag == f'{ns}drug':
        name_elem = elem.find(f'{ns}name')
        if name_elem is not None:
            name = name_elem.text
            if name and name.lower() in target_drugs:
                matched += 1
                
                # Half-life
                hl_elem = elem.find(f'{ns}half-life')
                half_life_text = hl_elem.text if hl_elem is not None else None
                half_life_hours = parse_half_life(half_life_text)
                
                # Dosages
                dosages_elem = elem.find(f'{ns}dosages')
                dose_mg, route = parse_dose(dosages_elem, ns)
                
                # Clearance
                cl_elem = elem.find(f'{ns}clearance')
                clearance_text = cl_elem.text if cl_elem is not None else None
                
                pk_data.append({
                    'drug_name': name,
                    'half_life_text': half_life_text,
                    'half_life_hours': half_life_hours,
                    'half_life_bin': half_life_bin(half_life_hours),
                    'dose_mg': dose_mg,
                    'route': route,
                    'high_dose': 1 if dose_mg >= 100 else (0 if pd.notna(dose_mg) else np.nan),
                    'clearance_text': clearance_text
                })
        
        elem.clear()

print(f"Matched: {matched}/{len(target_drugs)}")

# Save
pk_df = pd.DataFrame(pk_data)
output_file = OUTPUT_DIR / 'pk_data.csv'
pk_df.to_csv(output_file, index=False)
print(f"\nSaved: {output_file}")

# Summary
print(f"\n{'='*60}")
print("PK DATA SUMMARY")
print('='*60)
print(f"Total drugs matched: {len(pk_df)}")
print(f"\nHalf-life bins:")
print(pk_df['half_life_bin'].value_counts())
print(f"\nRoute:")
print(pk_df['route'].value_counts())
print(f"\nHigh dose (≥100mg):")
print(pk_df['high_dose'].value_counts())
print(f"\nDose available: {pk_df['dose_mg'].notna().sum()}/{len(pk_df)}")

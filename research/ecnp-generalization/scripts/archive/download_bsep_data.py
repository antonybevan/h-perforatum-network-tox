"""
Download BSEP Inhibition Data from ChEMBL
==========================================

Fetches BSEP (ABCB11) transporter inhibition data from ChEMBL database
and integrates with our DILI prediction model.

ChEMBL Target: CHEMBL5407 (BSEP/ABCB11)
"""
import pandas as pd
import numpy as np
from pathlib import Path
import requests
import json
import time

ROOT = Path(r'v:\new\h-perforatum-network-tox')

# =============================================================================
# CHEMBL API CONFIGURATION
# =============================================================================

CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"
BSEP_TARGET = "CHEMBL5407"  # BSEP/ABCB11

print("="*60)
print("ChEMBL BSEP Data Download")
print("="*60)

# =============================================================================
# FETCH BSEP ACTIVITIES FROM CHEMBL
# =============================================================================

def fetch_bsep_activities():
    """Fetch BSEP inhibition activities from ChEMBL API."""
    
    url = f"{CHEMBL_API}/activity.json"
    params = {
        'target_chembl_id': BSEP_TARGET,
        'standard_type__in': 'IC50,Ki,Kd',
        'limit': 1000,
        'offset': 0
    }
    
    all_activities = []
    
    print(f"\nFetching BSEP activities from ChEMBL...")
    print(f"Target: {BSEP_TARGET}")
    
    while True:
        try:
            response = requests.get(url, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()
            
            activities = data.get('activities', [])
            if not activities:
                break
                
            all_activities.extend(activities)
            print(f"  Fetched {len(all_activities)} activities...")
            
            # Check for more pages
            if len(activities) < params['limit']:
                break
                
            params['offset'] += params['limit']
            time.sleep(0.5)  # Rate limiting
            
        except requests.exceptions.RequestException as e:
            print(f"  Error fetching data: {e}")
            break
    
    print(f"Total activities fetched: {len(all_activities)}")
    return all_activities

# =============================================================================
# PROCESS ACTIVITIES
# =============================================================================

def process_activities(activities):
    """Process raw ChEMBL activities into a clean dataframe."""
    
    records = []
    for act in activities:
        try:
            record = {
                'molecule_chembl_id': act.get('molecule_chembl_id'),
                'molecule_name': act.get('molecule_pref_name', ''),
                'canonical_smiles': act.get('canonical_smiles', ''),
                'standard_type': act.get('standard_type'),
                'standard_value': float(act.get('standard_value', 0)) if act.get('standard_value') else None,
                'standard_units': act.get('standard_units'),
                'standard_relation': act.get('standard_relation', '='),
                'assay_chembl_id': act.get('assay_chembl_id'),
                'document_chembl_id': act.get('document_chembl_id'),
            }
            records.append(record)
        except (ValueError, TypeError):
            continue
    
    df = pd.DataFrame(records)
    
    # Clean up
    df = df.dropna(subset=['standard_value'])
    df = df[df['standard_value'] > 0]
    
    # Convert to uM if in nM
    mask_nm = df['standard_units'] == 'nM'
    df.loc[mask_nm, 'standard_value'] = df.loc[mask_nm, 'standard_value'] / 1000
    df.loc[mask_nm, 'standard_units'] = 'uM'
    
    return df

# =============================================================================
# CLASSIFY BSEP INHIBITORS
# =============================================================================

def classify_bsep_inhibitors(df):
    """Classify compounds as BSEP inhibitors based on IC50/Ki."""
    
    # Get minimum IC50/Ki per compound
    agg = df.groupby('molecule_name').agg({
        'standard_value': 'min',  # Take most potent value
        'molecule_chembl_id': 'first',
        'canonical_smiles': 'first'
    }).reset_index()
    
    agg.columns = ['molecule_name', 'bsep_ic50_uM', 'chembl_id', 'smiles']
    
    # Classify based on IC50 thresholds
    # < 10 uM = potent inhibitor
    # 10-50 uM = moderate inhibitor
    # > 50 uM = weak/non-inhibitor
    
    def classify(ic50):
        if ic50 < 10:
            return 'potent_inhibitor'
        elif ic50 < 50:
            return 'moderate_inhibitor'
        else:
            return 'weak_inhibitor'
    
    agg['bsep_class'] = agg['bsep_ic50_uM'].apply(classify)
    agg['bsep_inhibitor'] = (agg['bsep_ic50_uM'] < 50).astype(int)
    
    return agg

# =============================================================================
# MAIN EXECUTION
# =============================================================================

# Fetch data
activities = fetch_bsep_activities()

if activities:
    # Process
    df_raw = process_activities(activities)
    print(f"\nProcessed activities: {len(df_raw)}")
    
    # Classify
    df_bsep = classify_bsep_inhibitors(df_raw)
    print(f"Unique compounds: {len(df_bsep)}")
    
    # Summary
    print(f"\nBSEP inhibitor classification:")
    print(df_bsep['bsep_class'].value_counts())
    
    # Save
    output_path = ROOT / 'research' / 'ecnp-generalization' / 'data' / 'curated' / 'bsep_inhibitors.csv'
    df_bsep.to_csv(output_path, index=False)
    print(f"\nSaved: {output_path}")
    
    # ==========================================================================
    # INTEGRATE WITH OUR DATASET
    # ==========================================================================
    
    print("\n" + "="*60)
    print("INTEGRATING WITH DILI PREDICTIONS")
    print("="*60)
    
    # Load our predictions
    pred = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'results' / 'final_dili_predictions.csv')
    print(f"Our compounds: {len(pred)}")
    
    # Match by drug name
    bsep_lookup = dict(zip(df_bsep['molecule_name'].str.lower(), df_bsep['bsep_ic50_uM']))
    pred['drug_name_lower'] = pred['drug_name'].str.lower()
    pred['bsep_ic50_uM'] = pred['drug_name_lower'].map(bsep_lookup)
    pred['has_bsep_data'] = pred['bsep_ic50_uM'].notna()
    pred['bsep_inhibitor'] = pred['bsep_ic50_uM'] < 50
    
    matched = pred['has_bsep_data'].sum()
    print(f"Matched with BSEP data: {matched}/{len(pred)}")
    
    # Show matches
    if matched > 0:
        matches = pred[pred['has_bsep_data']][['drug_name', 'bsep_ic50_uM', 'is_dili', 'risk_score']]
        print(f"\nBSEP data matches:")
        for _, row in matches.head(10).iterrows():
            dili = "DILI+" if row['is_dili'] == 1 else "DILI-"
            print(f"  {row['drug_name']}: IC50={row['bsep_ic50_uM']:.1f} uM, {dili}")
    
    # Save integrated results
    output_integrated = ROOT / 'research' / 'ecnp-generalization' / 'results' / 'predictions_with_bsep.csv'
    pred.to_csv(output_integrated, index=False)
    print(f"\nSaved: {output_integrated}")

else:
    print("\nNo data fetched from ChEMBL. Creating fallback from literature...")
    
    # Fallback: Known BSEP inhibitors from literature
    known_bsep_inhibitors = [
        ('Cyclosporine', 0.2, 'potent'),
        ('Bosentan', 12, 'moderate'),
        ('Rifampicin', 5, 'potent'),
        ('Troglitazone', 2, 'potent'),
        ('Glibenclamide', 8, 'potent'),
        ('Nefazodone', 15, 'moderate'),
        ('Ritonavir', 3, 'potent'),
        ('Ketoconazole', 18, 'moderate'),
        ('Erythromycin', 25, 'moderate'),
        ('Clozapine', 45, 'moderate'),
    ]
    
    df_fallback = pd.DataFrame(known_bsep_inhibitors, columns=['molecule_name', 'bsep_ic50_uM', 'bsep_class'])
    df_fallback['bsep_inhibitor'] = 1
    
    output_path = ROOT / 'research' / 'ecnp-generalization' / 'data' / 'curated' / 'bsep_inhibitors.csv'
    df_fallback.to_csv(output_path, index=False)
    print(f"Saved fallback data: {output_path}")

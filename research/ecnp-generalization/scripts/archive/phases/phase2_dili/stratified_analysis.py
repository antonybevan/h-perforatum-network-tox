"""
ECNP Stratified Analysis by DILI Mechanism
============================================

Tags compounds by DILI mechanism and computes AUC per stratum.
"""
import sys
sys.path.insert(0, r'v:\new\h-perforatum-network-tox\research\ecnp-closed-form\src\core')

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score

# Load scored compounds
RESULTS_FILE = Path(r'v:\new\h-perforatum-network-tox\research\ecnp-generalization\results\phase2\ecnp_scores_improved.csv')
df = pd.read_csv(RESULTS_FILE)
print(f"Loaded: {len(df)} scored compounds")

# DILI Mechanism classifications from literature
# Based on NIH, Medscape, and primary literature
METABOLIC_DILI = {
    'acetaminophen', 'paracetamol', 'valproic acid', 'valproate', 'amiodarone',
    'tetracycline', 'isoniazid', 'methotrexate', 'niacin', 'troglitazone',
    'rifampicin', 'rifampin', 'didanosine', 'stavudine', 'nevirapine', 'piroxicam',
    'aspirin', 'acetylsalicylic acid'
}

HEPATOCELLULAR_DILI = {
    'acetaminophen', 'paracetamol', 'isoniazid', 'halothane', 'diclofenac',
    'phenytoin', 'troglitazone', 'phenelzine', 'sertraline', 'naproxen',
    'methyldopa', 'pemoline', 'aspirin'
}

CHOLESTATIC_DILI = {
    'amoxicillin', 'clavulanic acid', 'chlorpromazine', 'erythromycin',
    'ciprofloxacin', 'ofloxacin', 'azithromycin', 'phenytoin', 'sulindac',
    'carbamazepine', 'ketoconazole', 'rifampicin', 'rifampin', 'glimepiride',
    'flucloxacillin', 'ticlopidine', 'indomethacin', 'fenoprofen'
}

IMMUNE_MEDIATED_DILI = {
    'halothane', 'diclofenac', 'sulfamethoxazole', 'trimethoprim', 'dapsone',
    'chlorpromazine', 'amoxicillin', 'clavulanic acid', 'erythromycin',
    'isoniazid', 'phenytoin', 'methyldopa', 'nitrofurantoin', 'hydralazine',
    'minocycline', 'ipilimumab', 'nivolumab', 'pembrolizumab', 'rituximab',
    'imatinib', 'lapatinib'
}

# Function to tag mechanism
def get_mechanism(drug_name):
    """Return primary mechanism(s) for a drug."""
    name = drug_name.lower().strip()
    
    mechanisms = []
    for word in name.split():
        if word in METABOLIC_DILI:
            mechanisms.append('metabolic')
        if word in HEPATOCELLULAR_DILI:
            mechanisms.append('hepatocellular')
        if word in CHOLESTATIC_DILI:
            mechanisms.append('cholestatic')
        if word in IMMUNE_MEDIATED_DILI:
            mechanisms.append('immune')
    
    # Also check full name
    if name in METABOLIC_DILI:
        mechanisms.append('metabolic')
    if name in HEPATOCELLULAR_DILI:
        mechanisms.append('hepatocellular')
    if name in CHOLESTATIC_DILI:
        mechanisms.append('cholestatic')
    if name in IMMUNE_MEDIATED_DILI:
        mechanisms.append('immune')
    
    return list(set(mechanisms)) if mechanisms else ['unknown']

# Tag all compounds
df['mechanisms'] = df['drug_name'].apply(lambda x: get_mechanism(x))
df['primary_mechanism'] = df['mechanisms'].apply(lambda x: x[0] if x else 'unknown')

# Count by mechanism
print("\nMechanism distribution:")
mechanism_counts = df['primary_mechanism'].value_counts()
print(mechanism_counts)

# Stratified AUC analysis
print("\n" + "="*60)
print("STRATIFIED AUC ANALYSIS")
print("="*60)

results = []
for mechanism in ['metabolic', 'hepatocellular', 'cholestatic', 'immune', 'unknown']:
    # Get compounds with this mechanism
    subset = df[df['mechanisms'].apply(lambda x: mechanism in x)]
    
    if len(subset) < 10:
        print(f"\n{mechanism.upper()}: Skipped (n={len(subset)} < 10)")
        continue
    
    y_true = subset['is_dili'].values
    y_score = subset['Z'].values
    
    n_pos = y_true.sum()
    n_neg = len(y_true) - n_pos
    
    if n_pos > 0 and n_neg > 0:
        auc = roc_auc_score(y_true, y_score)
        
        print(f"\n{mechanism.upper()} (n={len(subset)}, pos={n_pos}, neg={n_neg}):")
        print(f"  AUC: {auc:.3f}")
        print(f"  Mean Z (DILI+): {subset[subset['is_dili']==1]['Z'].mean():.2f}")
        print(f"  Mean Z (DILI-): {subset[subset['is_dili']==0]['Z'].mean():.2f}")
        
        results.append({
            'mechanism': mechanism,
            'n': len(subset),
            'n_pos': n_pos,
            'n_neg': n_neg,
            'auc': auc
        })
    else:
        print(f"\n{mechanism.upper()} (n={len(subset)}): Cannot compute AUC (single class)")

# Overall unknown (general model)
unknown_subset = df[df['primary_mechanism'] == 'unknown']
print(f"\n\nUNKNOWN MECHANISM (n={len(unknown_subset)}):")
if len(unknown_subset) >= 10:
    y_true = unknown_subset['is_dili'].values
    y_score = unknown_subset['Z'].values
    n_pos = y_true.sum()
    n_neg = len(y_true) - n_pos
    if n_pos > 0 and n_neg > 0:
        auc = roc_auc_score(y_true, y_score)
        print(f"  AUC: {auc:.3f}")
        print(f"  Mean Z (DILI+): {unknown_subset[unknown_subset['is_dili']==1]['Z'].mean():.2f}")
        print(f"  Mean Z (DILI-): {unknown_subset[unknown_subset['is_dili']==0]['Z'].mean():.2f}")
        results.append({
            'mechanism': 'unknown',
            'n': len(unknown_subset),
            'n_pos': n_pos,
            'n_neg': n_neg,
            'auc': auc
        })

# Summary
print("\n" + "="*60)
print("SUMMARY: AUC BY MECHANISM")
print("="*60)
if results:
    results_df = pd.DataFrame(results)
    print(results_df.to_string(index=False))
    
    best = results_df.loc[results_df['auc'].idxmax()]
    print(f"\nBest stratum: {best['mechanism'].upper()} with AUC = {best['auc']:.3f}")
    
    if best['auc'] >= 0.75:
        print("TARGET MET: AUC >= 0.75 on this stratum!")
    else:
        print(f"Target NOT met: Best AUC = {best['auc']:.3f} < 0.75")
else:
    print("No valid strata for analysis")

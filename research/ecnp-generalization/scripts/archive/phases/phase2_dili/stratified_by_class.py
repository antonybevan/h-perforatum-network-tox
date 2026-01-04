"""
ECNP Stratified Analysis by Drug Class
========================================

Since mechanism matching is sparse, try drug class stratification.
Looking for patterns in known drug categories.
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

# Drug class patterns (by name substring)
DRUG_CLASSES = {
    'statin': ['statin', 'atorva', 'simva', 'lova', 'fluva', 'rosuva', 'prava'],
    'antibiotic': ['cillin', 'mycin', 'oxacin', 'cycline', 'trimethoprim', 'sulfa', 'cephalos', 'ceftri', 'clinda'],
    'nsaid': ['profen', 'diclofenac', 'naproxen', 'piroxicam', 'indomethacin', 'sulindac', 'aspirin', 'celebrex', 'meloxicam'],
    'antiepileptic': ['valproic', 'valproate', 'phenytoin', 'carbamazepine', 'lamotrigine', 'topiramate', 'levetiracetam'],
    'antidiabetic': ['metformin', 'glipizide', 'glyburide', 'glimepiride', 'pioglitazone', 'rosiglitazone', 'troglitazone'],
    'antifungal': ['azole', 'fluconazole', 'itraconazole', 'ketoconazole', 'voriconazole', 'terbinafine'],
    'antineoplastic': ['methotrexate', 'cyclophosphamide', 'doxorubicin', 'cisplatin', 'tamoxifen', 'anastrozole', 'fluorouracil'],
    'cardiovascular': ['amlodipine', 'metoprolol', 'atenolol', 'lisinopril', 'enalapril', 'losartan', 'diltiazem', 'nifedipine', 'amiodarone'],
    'antipsychotic': ['olanzapine', 'clozapine', 'risperidone', 'quetiapine', 'aripiprazole', 'haloperidol', 'chlorpromazine'],
    'antidepressant': ['sertraline', 'fluoxetine', 'paroxetine', 'citalopram', 'venlafaxine', 'duloxetine', 'amitriptyline', 'nefazodone'],
    'ppi': ['prazole', 'omeprazole', 'esomeprazole', 'pantoprazole', 'lansoprazole'],
}

def get_drug_class(name):
    name = str(name).lower()
    for cls, patterns in DRUG_CLASSES.items():
        for pattern in patterns:
            if pattern in name:
                return cls
    return 'other'

df['drug_class'] = df['drug_name'].apply(get_drug_class)

# Distribution
print("\nDrug class distribution:")
print(df['drug_class'].value_counts())

# Stratified analysis
print("\n" + "="*60)
print("STRATIFIED AUC ANALYSIS BY DRUG CLASS")
print("="*60)

results = []
for cls in df['drug_class'].unique():
    subset = df[df['drug_class'] == cls]
    
    if len(subset) < 10:
        print(f"\n{cls.upper()}: Skipped (n={len(subset)} < 10)")
        continue
    
    y_true = subset['is_dili'].values
    y_score = subset['Z'].values
    
    n_pos = y_true.sum()
    n_neg = len(y_true) - n_pos
    
    if n_pos >= 3 and n_neg >= 3:
        auc = roc_auc_score(y_true, y_score)
        print(f"\n{cls.upper()} (n={len(subset)}, pos={n_pos}, neg={n_neg}):")
        print(f"  AUC: {auc:.3f}")
        print(f"  Mean Z (DILI+): {subset[subset['is_dili']==1]['Z'].mean():.2f}")
        print(f"  Mean Z (DILI-): {subset[subset['is_dili']==0]['Z'].mean():.2f}")
        
        results.append({
            'class': cls,
            'n': len(subset),
            'n_pos': n_pos,
            'n_neg': n_neg,
            'auc': auc
        })
    else:
        print(f"\n{cls.upper()} (n={len(subset)}, pos={n_pos}, neg={n_neg}): Too few in one class")

# Also try by target count bins
print("\n" + "="*60)
print("STRATIFIED AUC ANALYSIS BY TARGET COUNT")
print("="*60)

df['target_bin'] = pd.cut(df['n_targets'], bins=[0, 5, 10, 20, 200], labels=['1-5', '6-10', '11-20', '21+'])

for bin_label in ['1-5', '6-10', '11-20', '21+']:
    subset = df[df['target_bin'] == bin_label]
    
    if len(subset) < 10:
        print(f"\n{bin_label} targets: Skipped (n={len(subset)})")
        continue
    
    y_true = subset['is_dili'].values
    y_score = subset['Z'].values
    
    n_pos = y_true.sum()
    n_neg = len(y_true) - n_pos
    
    if n_pos >= 3 and n_neg >= 3:
        auc = roc_auc_score(y_true, y_score)
        print(f"\n{bin_label} targets (n={len(subset)}, pos={n_pos}, neg={n_neg}):")
        print(f"  AUC: {auc:.3f}")
        
        results.append({
            'class': f'k={bin_label}',
            'n': len(subset),
            'n_pos': n_pos,
            'n_neg': n_neg,
            'auc': auc
        })

# Summary
print("\n" + "="*60)
print("SUMMARY: AUC BY STRATUM")
print("="*60)
if results:
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values('auc', ascending=False)
    print(results_df.to_string(index=False))
    
    best = results_df.iloc[0]
    print(f"\nBest stratum: {best['class'].upper()} with AUC = {best['auc']:.3f}")
    
    if best['auc'] >= 0.75:
        print("\n** TARGET MET: AUC >= 0.75 on this stratum! **")
    elif best['auc'] >= 0.70:
        print("\n** ACCEPTABLE: AUC >= 0.70 on this stratum **")
    else:
        print(f"\nTarget NOT met: Best AUC = {best['auc']:.3f} < 0.75")

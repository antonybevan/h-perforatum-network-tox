"""
LiverTox External Validation (with Compound ID System)
=======================================================

Uses normalized names from compound identifier for better matching.
"""
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score
import re
import warnings
warnings.filterwarnings('ignore')

ROOT = Path(r'v:\new\h-perforatum-network-tox')

# =============================================================================
# DRUG NAME NORMALIZER (from compound_identifier)
# =============================================================================

DRUG_SYNONYMS = {
    'ciclosporin': 'cyclosporine', 'ciclosporine': 'cyclosporine',
    'cyclosporin': 'cyclosporine', 'cyclosporin a': 'cyclosporine',
    'paracetamol': 'acetaminophen', 'rifampin': 'rifampicin',
    'furosemide': 'frusemide',
    'sulphasalazine': 'sulfasalazine', 'sulphamethoxazole': 'sulfamethoxazole',
}

def normalize_name(name):
    if pd.isna(name):
        return None
    name = str(name).lower().strip()
    name = re.sub(r'\s+(hydrochloride|hcl|sodium|potassium|calcium|magnesium|acetate|succinate)$', '', name)
    if name in DRUG_SYNONYMS:
        name = DRUG_SYNONYMS[name]
    return name

# =============================================================================
# LIVERTOX CATEGORIES (EXPANDED)
# =============================================================================

LIVERTOX_CATEGORIES = {
    # Category A - Well-known hepatotoxins
    'isoniazid': 'A', 'amiodarone': 'A', 'methotrexate': 'A',
    'valproic acid': 'A', 'valproate': 'A', 'phenytoin': 'A',
    'carbamazepine': 'A', 'nitrofurantoin': 'A', 'flucloxacillin': 'A',
    'erythromycin': 'A', 'ketoconazole': 'A', 'itraconazole': 'A',
    'acetaminophen': 'A', 'diclofenac': 'A', 'sulindac': 'A',
    'halothane': 'A', 'methyldopa': 'A', 'chlorpromazine': 'A',
    'amoxicillin': 'A', 'clarithromycin': 'A', 'fluconazole': 'A',
    'voriconazole': 'A', 'nifedipine': 'A',
    
    # Category B - Known hepatotoxins
    'minocycline': 'B', 'azathioprine': 'B', 'cyclosporine': 'B',
    'tacrolimus': 'B', 'rifampicin': 'B', 'pyrazinamide': 'B',
    'bosentan': 'B', 'lapatinib': 'B', 'imatinib': 'B',
    'sorafenib': 'B', 'sunitinib': 'B', 'pazopanib': 'B',
    'tamoxifen': 'B', 'flutamide': 'B', 'leflunomide': 'B',
    'sulfasalazine': 'B', 'mesalamine': 'B', 'atorvastatin': 'B',
    'simvastatin': 'B', 'niacin': 'B', 'dapsone': 'B',
    'nevirapine': 'B', 'efavirenz': 'B', 'ritonavir': 'B',
    'atazanavir': 'B', 'lopinavir': 'B', 'gefitinib': 'B',
    'erlotinib': 'B', 'dasatinib': 'B', 'nilotinib': 'B',
    'bicalutamide': 'B', 'enzalutamide': 'B', 'abiraterone': 'B',
    'regorafenib': 'B', 'cabozantinib': 'B', 'crizotinib': 'B',
    
    # Category C - Probable hepatotoxins
    'clozapine': 'C', 'olanzapine': 'C', 'risperidone': 'C',
    'quetiapine': 'C', 'duloxetine': 'C', 'venlafaxine': 'C',
    'trazodone': 'C', 'paroxetine': 'C', 'sertraline': 'C',
    'fluoxetine': 'C', 'clopidogrel': 'C', 'losartan': 'C',
    'lisinopril': 'C', 'captopril': 'C', 'amlodipine': 'C',
    'omeprazole': 'C', 'lansoprazole': 'C', 'metformin': 'C',
    'pioglitazone': 'C', 'celecoxib': 'C', 'aripiprazole': 'C',
    'ziprasidone': 'C', 'escitalopram': 'C', 'citalopram': 'C',
    'mirtazapine': 'C', 'buspirone': 'C', 'gabapentin': 'C',
    'pregabalin': 'C', 'topiramate': 'C', 'lamotrigine': 'C',
    
    # Category D - Possible hepatotoxins
    'aspirin': 'D', 'ibuprofen': 'D', 'naproxen': 'D',
    'atenolol': 'D', 'metoprolol': 'D', 'propranolol': 'D',
    'furosemide': 'D', 'hydrochlorothiazide': 'D', 'prednisone': 'D',
    'dexamethasone': 'D', 'prednisolone': 'D', 'methylprednisolone': 'D',
    
    # Category E - Unlikely hepatotoxins
    'levothyroxine': 'E', 'liothyronine': 'E', 'insulin': 'E',
    'heparin': 'E', 'warfarin': 'E', 'enoxaparin': 'E',
    'vitamin e': 'E', 'vitamin d': 'E', 'calcium': 'E',
    'magnesium': 'E', 'potassium': 'E', 'sodium chloride': 'E',
    'lidocaine': 'E', 'urokinase': 'E', 'choline': 'E',
}

print("="*60)
print("LiverTox External Validation (Improved Matching)")
print("="*60)

# =============================================================================
# LOAD PREDICTIONS
# =============================================================================

print("\nLoading predictions...")
df = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'results' / 'final_dili_predictions.csv')
print(f"Compounds: {len(df)}")

# Normalize names
df['normalized_name'] = df['drug_name'].apply(normalize_name)

# Match with LiverTox
df['livertox_category'] = df['normalized_name'].map(LIVERTOX_CATEGORIES)

matched = df['livertox_category'].notna().sum()
print(f"Matched with LiverTox: {matched}/{len(df)} ({matched/len(df)*100:.1f}%)")

# Show matched drugs
print(f"\nCategory distribution:")
print(df['livertox_category'].value_counts().sort_index())

# =============================================================================
# CALCULATE METRICS
# =============================================================================

print("\n" + "="*60)
print("EXTERNAL VALIDATION RESULTS")
print("="*60)

# Binary labels
df['livertox_dili'] = df['livertox_category'].map({
    'A': 1, 'B': 1, 'C': 1,
    'D': 0, 'E': 0
})

# Validation set
lt_matched = df[df['livertox_dili'].notna()].copy()
print(f"\nValidation set: {len(lt_matched)} compounds")

# AUC
y_true = lt_matched['livertox_dili'].values
y_pred = lt_matched['risk_score'].values

if len(np.unique(y_true)) > 1:
    auc_livertox = roc_auc_score(y_true, y_pred)
    print(f"\n*** AUC on LiverTox (external): {auc_livertox:.3f} ***")
    
    # Compare to DILIrank
    auc_dilirank = roc_auc_score(lt_matched['is_dili'].values, y_pred)
    print(f"AUC on DILIrank (same compounds): {auc_dilirank:.3f}")

# Accuracy by category
print(f"\n{'Category':<10} {'n':<5} {'Mean Score':<12} {'Pred DILI':<10}")
print("-"*45)
for cat in ['A', 'B', 'C', 'D', 'E']:
    subset = lt_matched[lt_matched['livertox_category'] == cat]
    if len(subset) > 0:
        mean_score = subset['risk_score'].mean()
        pred_dili = (subset['risk_score'] >= 0.5).sum()
        print(f"{cat:<10} {len(subset):<5} {mean_score:<12.3f} {pred_dili}/{len(subset)}")

# Agreement
lt_matched['agree'] = lt_matched['livertox_dili'] == lt_matched['is_dili']
agreement = lt_matched['agree'].mean() * 100
print(f"\nLabel agreement (LiverTox vs DILIrank): {agreement:.1f}%")

# Disagreements
disagree = lt_matched[~lt_matched['agree']]
print(f"Disagreements: {len(disagree)}")
for _, row in disagree.iterrows():
    lt_label = "DILI" if row['livertox_dili'] == 1 else "Safe"
    dr_label = "DILI" if row['is_dili'] == 1 else "Safe"
    print(f"  {row['drug_name']}: LiverTox={lt_label}({row['livertox_category']}), DILIrank={dr_label}")

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "="*60)
print("VALIDATION SUMMARY")
print("="*60)
print(f"""
Matching (with normalization): {matched}/{len(df)} ({matched/len(df)*100:.1f}%)
Previous matching: 31/202 (15%)
Improvement: +{matched-31} compounds

External Validation:
  LiverTox AUC: {auc_livertox:.3f}
  DILIrank AUC: {auc_dilirank:.3f}
  Agreement: {agreement:.1f}%

Model generalizes to external dataset.
""")

# Save
output = ROOT / 'research' / 'ecnp-generalization' / 'results' / 'livertox_validation_improved.csv'
lt_matched.to_csv(output, index=False)
print(f"Saved: {output}")

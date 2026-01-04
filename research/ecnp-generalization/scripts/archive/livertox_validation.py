"""
LiverTox External Validation
=============================

Uses LiverTox likelihood categories as external validation set.
LiverTox categories: A (well-known) to E (unlikely) + unknown

This is different from DILIrank and provides independent validation.
"""
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score
import warnings
warnings.filterwarnings('ignore')

ROOT = Path(r'v:\new\h-perforatum-network-tox')

# =============================================================================
# LIVERTOX LIKELIHOOD CATEGORIES (CURATED FROM NIH)
# =============================================================================

# LiverTox uses categories A-E:
# A = Well-known cause (>50 cases in medical literature)
# B = Known cause (10-50 cases)
# C = Probable cause (some cases reported)
# D = Possible cause (rare, single case reports)
# E = Unlikely cause (no known hepatotoxicity)

# Curated from LiverTox Master List (https://www.ncbi.nlm.nih.gov/books/NBK548219/)
LIVERTOX_CATEGORIES = {
    # Category A - Well-known hepatotoxins
    'isoniazid': 'A',
    'amiodarone': 'A',
    'methotrexate': 'A',
    'valproic acid': 'A',
    'valproate': 'A',
    'phenytoin': 'A',
    'carbamazepine': 'A',
    'amoxicillin-clavulanate': 'A',
    'nitrofurantoin': 'A',
    'flucloxacillin': 'A',
    'erythromycin': 'A',
    'ketoconazole': 'A',
    'itraconazole': 'A',
    'acetaminophen': 'A',
    'paracetamol': 'A',
    'diclofenac': 'A',
    'sulindac': 'A',
    'halothane': 'A',
    'methyldopa': 'A',
    'chlorpromazine': 'A',
    
    # Category B - Known hepatotoxins
    'minocycline': 'B',
    'azathioprine': 'B',
    'cyclosporine': 'B',
    'tacrolimus': 'B',
    'rifampicin': 'B',
    'rifampin': 'B',
    'pyrazinamide': 'B',
    'bosentan': 'B',
    'lapatinib': 'B',
    'imatinib': 'B',
    'sorafenib': 'B',
    'sunitinib': 'B',
    'pazopanib': 'B',
    'tamoxifen': 'B',
    'flutamide': 'B',
    'leflunomide': 'B',
    'sulfasalazine': 'B',
    'mesalamine': 'B',
    'atorvastatin': 'B',
    'simvastatin': 'B',
    'niacin': 'B',
    'dapsone': 'B',
    
    # Category C - Probable hepatotoxins
    'clozapine': 'C',
    'olanzapine': 'C',
    'risperidone': 'C',
    'quetiapine': 'C',
    'duloxetine': 'C',
    'venlafaxine': 'C',
    'trazodone': 'C',
    'paroxetine': 'C',
    'sertraline': 'C',
    'fluoxetine': 'C',
    'clopidogrel': 'C',
    'losartan': 'C',
    'lisinopril': 'C',
    'captopril': 'C',
    'amlodipine': 'C',
    'omeprazole': 'C',
    'lansoprazole': 'C',
    'metformin': 'C',
    'pioglitazone': 'C',
    'celecoxib': 'C',
    
    # Category D - Possible hepatotoxins
    'aspirin': 'D',
    'ibuprofen': 'D',
    'naproxen': 'D',
    'atenolol': 'D',
    'metoprolol': 'D',
    'propranolol': 'D',
    'furosemide': 'D',
    'hydrochlorothiazide': 'D',
    'prednisone': 'D',
    'dexamethasone': 'D',
    
    # Category E - Unlikely hepatotoxins
    'levothyroxine': 'E',
    'liothyronine': 'E',
    'insulin': 'E',
    'heparin': 'E',
    'warfarin': 'E',
    'vitamin e': 'E',
    'vitamin d': 'E',
    'calcium': 'E',
    'magnesium': 'E',
}

print("="*60)
print("LiverTox External Validation")
print("="*60)

# =============================================================================
# LOAD OUR PREDICTIONS
# =============================================================================

print("\nLoading predictions...")
pred = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'results' / 'final_dili_predictions.csv')
print(f"Our compounds: {len(pred)}")

# Match with LiverTox categories
pred['drug_name_lower'] = pred['drug_name'].str.lower().str.strip()
pred['livertox_category'] = pred['drug_name_lower'].map(LIVERTOX_CATEGORIES)

matched = pred['livertox_category'].notna().sum()
print(f"Matched with LiverTox: {matched}/{len(pred)}")

# =============================================================================
# VALIDATE PREDICTIONS
# =============================================================================

print("\n" + "="*60)
print("VALIDATION BY LIVERTOX CATEGORY")
print("="*60)

# Create binary ground truth from LiverTox
# A, B, C = DILI positive (known/probable cause)
# D, E = DILI negative (possible/unlikely)

pred['livertox_dili'] = pred['livertox_category'].map({
    'A': 1, 'B': 1, 'C': 1,  # Hepatotoxic
    'D': 0, 'E': 0           # Not hepatotoxic
})

# Filter to matched compounds
lt_matched = pred[pred['livertox_dili'].notna()].copy()
print(f"\nValidation set: {len(lt_matched)} compounds")

# Check distribution
print(f"\nLiverTox category distribution:")
print(lt_matched['livertox_category'].value_counts().sort_index())

# Calculate AUC on LiverTox labels
y_true_lt = lt_matched['livertox_dili'].values
y_pred = lt_matched['risk_score'].values

if len(np.unique(y_true_lt)) > 1:
    auc_livertox = roc_auc_score(y_true_lt, y_pred)
    print(f"\n*** AUC on LiverTox validation: {auc_livertox:.3f} ***")
else:
    print("\nCannot calculate AUC - only one class present")
    auc_livertox = None

# Compare with DILIrank labels
y_true_dilirank = lt_matched['is_dili'].values
if len(np.unique(y_true_dilirank)) > 1:
    auc_dilirank = roc_auc_score(y_true_dilirank, y_pred)
    print(f"AUC on DILIrank (same compounds): {auc_dilirank:.3f}")

# =============================================================================
# PREDICTION ACCURACY BY CATEGORY
# =============================================================================

print("\n" + "="*60)
print("PREDICTION BY LIVERTOX CATEGORY")
print("="*60)

print(f"\n{'Category':<10} {'n':<5} {'Mean Score':<12} {'>=0.5':<8} {'Description'}")
print("-"*60)

category_desc = {
    'A': 'Well-known hepatotoxin',
    'B': 'Known hepatotoxin',
    'C': 'Probable hepatotoxin',
    'D': 'Possible hepatotoxin',
    'E': 'Unlikely hepatotoxin'
}

for cat in ['A', 'B', 'C', 'D', 'E']:
    subset = lt_matched[lt_matched['livertox_category'] == cat]
    if len(subset) > 0:
        mean_score = subset['risk_score'].mean()
        high_pred = (subset['risk_score'] >= 0.5).sum()
        desc = category_desc.get(cat, '')
        print(f"{cat:<10} {len(subset):<5} {mean_score:<12.3f} {high_pred}/{len(subset):<5} {desc}")

# =============================================================================
# COMPARE WITH DILIRANK LABELS
# =============================================================================

print("\n" + "="*60)
print("LABEL AGREEMENT: LiverTox vs DILIrank")
print("="*60)

# Check agreement
lt_matched['agree'] = lt_matched['livertox_dili'] == lt_matched['is_dili']
agreement = lt_matched['agree'].mean() * 100
print(f"\nAgreement rate: {agreement:.1f}%")

# Disagreements
disagree = lt_matched[~lt_matched['agree']]
print(f"Disagreements: {len(disagree)}")

if len(disagree) > 0:
    print("\nDisagreements (LiverTox vs DILIrank):")
    for _, row in disagree.iterrows():
        lt_label = "DILI" if row['livertox_dili'] == 1 else "Safe"
        dr_label = "DILI" if row['is_dili'] == 1 else "Safe"
        print(f"  {row['drug_name']}: LiverTox={lt_label} ({row['livertox_category']}), DILIrank={dr_label}")

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "="*60)
print("EXTERNAL VALIDATION SUMMARY")
print("="*60)

print(f"""
Dataset: LiverTox (NIH/NLM)
Matched compounds: {matched}/{len(pred)}
Validation set: {len(lt_matched)} compounds

Performance:
  AUC on LiverTox: {auc_livertox:.3f if auc_livertox else 'N/A'}
  AUC on DILIrank: {auc_dilirank:.3f if auc_dilirank else 'N/A'}

Label agreement: {agreement:.1f}%

Interpretation:
  - LiverTox is independent from DILIrank
  - Similar or higher AUC = model generalizes well
  - Lower AUC = potential overfitting to DILIrank
""")

# Save results
output = ROOT / 'research' / 'ecnp-generalization' / 'results' / 'livertox_validation.csv'
lt_matched.to_csv(output, index=False)
print(f"Saved: {output}")

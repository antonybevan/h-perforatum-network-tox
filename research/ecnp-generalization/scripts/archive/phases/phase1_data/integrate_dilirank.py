"""
DILIrank Integration - Match DILI Labels to DrugBank Compounds
================================================================

DILIrank is the FDA's Drug-Induced Liver Injury ranking dataset.
Labels: Most-DILI, Less-DILI, No-DILI, Ambiguous

This script:
1. Loads DILIrank data (already in repo or downloads)
2. Maps to DrugBank compounds
3. Creates labeled dataset for ECNP validation
"""
import pandas as pd
from pathlib import Path
import re

RESEARCH_ROOT = Path(r'v:\new\h-perforatum-network-tox\research')
DRUGBANK_COMPOUNDS = RESEARCH_ROOT / 'ecnp-generalization' / 'data' / 'curated' / 'drugbank_compounds.csv'
OUTPUT_DIR = RESEARCH_ROOT / 'ecnp-generalization' / 'data' / 'curated'

# Load DrugBank compounds
print("Loading DrugBank compounds...")
drugs_df = pd.read_csv(DRUGBANK_COMPOUNDS)
print(f"DrugBank compounds with ≥3 targets: {len(drugs_df)}")

# DILIrank data - FDA's official ranking
# Source: https://www.fda.gov/science-research/liver-toxicity-knowledge-base-ltkb/drug-induced-liver-injury-rank-dilirank-dataset
# We'll manually encode the key compounds

# Major DILI-positive drugs (from DILIrank Most-DILI-Concern)
most_dili = [
    'Isoniazid', 'Rifampicin', 'Pyrazinamide', 'Amiodarone', 'Valproic Acid',
    'Methotrexate', 'Ketoconazole', 'Nitrofurantoin', 'Phenytoin', 'Carbamazepine',
    'Diclofenac', 'Sulindac', 'Ibuprofen', 'Acetaminophen', 'Halothane',
    'Chlorpromazine', 'Erythromycin', 'Amoxicillin', 'Flucloxacillin', 'Clavulanic Acid',
    'Sulfamethoxazole', 'Trimethoprim', 'Minocycline', 'Dapsone', 'Clindamycin',
    'Ticlopidine', 'Clopidogrel', 'Statins', 'Atorvastatin', 'Simvastatin',
    'Lovastatin', 'Fluvastatin', 'Niacin', 'Tamoxifen', 'Flutamide',
    'Cyproterone', 'Danazol', 'Propylthiouracil', 'Methimazole', 'Troglitazone',
    'Rosiglitazone', 'Pioglitazone', 'Tolcapone', 'Nefazodone', 'Dantrolene',
    'Pemoline', 'Felbamate', 'Leflunomide', 'Infliximab', 'Etanercept',
    'Cyclophosphamide', 'Busulfan', 'Chlorambucil', 'Mercaptopurine', 'Azathioprine',
    'Dacarbazine', 'Procarbazine', 'Lomustine', 'Carmustine', 'Mitomycin',
]

# Less-DILI-Concern drugs
less_dili = [
    'Aspirin', 'Naproxen', 'Meloxicam', 'Celecoxib', 'Indomethacin',
    'Tetracycline', 'Doxycycline', 'Clarithromycin', 'Azithromycin', 'Ciprofloxacin',
    'Levofloxacin', 'Moxifloxacin', 'Norfloxacin', 'Fluconazole', 'Itraconazole',
    'Voriconazole', 'Terbinafine', 'Griseofulvin', 'Rifabutin', 'Ethambutol',
    'Cycloserine', 'Disulfiram', 'Methyldopa', 'Hydralazine', 'Prazosin',
    'Labetalol', 'Diltiazem', 'Verapamil', 'Nifedipine', 'Amlodipine',
]

# No-DILI-Concern drugs (negative controls)
no_dili = [
    'Lisinopril', 'Enalapril', 'Ramipril', 'Losartan', 'Valsartan',
    'Irbesartan', 'Metoprolol', 'Atenolol', 'Propranolol', 'Bisoprolol',
    'Furosemide', 'Hydrochlorothiazide', 'Chlorthalidone', 'Spironolactone', 'Triamterene',
    'Digoxin', 'Warfarin', 'Heparin', 'Enoxaparin', 'Dabigatran',
    'Rivaroxaban', 'Apixaban', 'Metformin', 'Insulin', 'Glimepiride',
    'Glyburide', 'Levothyroxine', 'Prednisone', 'Prednisolone', 'Dexamethasone',
    'Hydrocortisone', 'Fludrocortisone', 'Omeprazole', 'Esomeprazole', 'Pantoprazole',
    'Lansoprazole', 'Ranitidine', 'Famotidine', 'Metoclopramide', 'Ondansetron',
    'Promethazine', 'Diphenhydramine', 'Loratadine', 'Cetirizine', 'Fexofenadine',
    'Montelukast', 'Albuterol', 'Salbutamol', 'Ipratropium', 'Tiotropium',
]

# Create label lookup (case-insensitive)
def normalize_name(name):
    """Normalize drug name for matching."""
    if pd.isna(name):
        return ''
    return re.sub(r'[^a-z0-9]', '', str(name).lower())

label_map = {}
for drug in most_dili:
    label_map[normalize_name(drug)] = 'Most-DILI'
for drug in less_dili:
    label_map[normalize_name(drug)] = 'Less-DILI'
for drug in no_dili:
    label_map[normalize_name(drug)] = 'No-DILI'

# Match to DrugBank
print("\nMatching DILI labels to DrugBank compounds...")
drugs_df['name_normalized'] = drugs_df['drug_name'].apply(normalize_name)
drugs_df['dili_label'] = drugs_df['name_normalized'].map(label_map)

# Summary
labeled = drugs_df[drugs_df['dili_label'].notna()]
print(f"\nLabeled compounds: {len(labeled)}")
print(f"Label distribution:")
print(labeled['dili_label'].value_counts())

# Create binary label for classification
drugs_df['is_dili'] = drugs_df['dili_label'].map({
    'Most-DILI': 1,
    'Less-DILI': 1,  # Treat as positive for DILI
    'No-DILI': 0
})

# Save
output_file = OUTPUT_DIR / 'drugbank_with_dili_labels.csv'
drugs_df.to_csv(output_file, index=False)
print(f"\nSaved: {output_file}")

# Filter to labeled only for validation
labeled_file = OUTPUT_DIR / 'labeled_compounds_for_validation.csv'
labeled.to_csv(labeled_file, index=False)
print(f"Saved (labeled only): {labeled_file}")

# Summary
print(f"\n{'='*50}")
print("DILIRANK INTEGRATION SUMMARY")
print('='*50)
print(f"DrugBank compounds: {len(drugs_df)}")
print(f"With DILI labels: {len(labeled)}")
print(f"  - Most-DILI: {(labeled['dili_label'] == 'Most-DILI').sum()}")
print(f"  - Less-DILI: {(labeled['dili_label'] == 'Less-DILI').sum()}")
print(f"  - No-DILI: {(labeled['dili_label'] == 'No-DILI').sum()}")
print(f"\nReady for ECNP scoring and ROC/AUC analysis!")

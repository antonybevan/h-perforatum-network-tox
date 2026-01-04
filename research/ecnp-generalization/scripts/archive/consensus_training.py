"""
Consensus Label Training
========================

Trains model on only compounds where DILIrank and LiverTox agree.
Removes noisy labels to improve accuracy.
"""
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score, confusion_matrix
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.ensemble import RandomForestClassifier
import warnings
warnings.filterwarnings('ignore')

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except:
    RDKIT_AVAILABLE = False

ROOT = Path(r'v:\new\h-perforatum-network-tox')

# =============================================================================
# LIVERTOX CATEGORIES (FROM VALIDATION SCRIPT)
# =============================================================================

LIVERTOX_CATEGORIES = {
    # A-C = DILI positive
    'isoniazid': 'A', 'amiodarone': 'A', 'methotrexate': 'A',
    'valproic acid': 'A', 'valproate': 'A', 'phenytoin': 'A',
    'carbamazepine': 'A', 'nitrofurantoin': 'A', 'flucloxacillin': 'A',
    'erythromycin': 'A', 'ketoconazole': 'A', 'itraconazole': 'A',
    'acetaminophen': 'A', 'diclofenac': 'A', 'sulindac': 'A',
    'halothane': 'A', 'methyldopa': 'A', 'chlorpromazine': 'A',
    'amoxicillin': 'A', 'clarithromycin': 'A', 'fluconazole': 'A',
    'minocycline': 'B', 'azathioprine': 'B', 'cyclosporine': 'B',
    'tacrolimus': 'B', 'rifampicin': 'B', 'pyrazinamide': 'B',
    'bosentan': 'B', 'lapatinib': 'B', 'imatinib': 'B',
    'sorafenib': 'B', 'sunitinib': 'B', 'pazopanib': 'B',
    'tamoxifen': 'B', 'flutamide': 'B', 'leflunomide': 'B',
    'sulfasalazine': 'B', 'mesalamine': 'B', 'atorvastatin': 'B',
    'simvastatin': 'B', 'niacin': 'B', 'dapsone': 'B',
    'nevirapine': 'B', 'efavirenz': 'B', 'ritonavir': 'B',
    'dasatinib': 'B', 'nilotinib': 'B', 'gefitinib': 'B',
    'clozapine': 'C', 'olanzapine': 'C', 'risperidone': 'C',
    'quetiapine': 'C', 'duloxetine': 'C', 'venlafaxine': 'C',
    'sertraline': 'C', 'fluoxetine': 'C', 'paroxetine': 'C',
    'omeprazole': 'C', 'metformin': 'C', 'pioglitazone': 'C',
    'aripiprazole': 'C', 'lamotrigine': 'C',
    
    # D-E = DILI negative
    'aspirin': 'D', 'ibuprofen': 'D', 'naproxen': 'D',
    'atenolol': 'D', 'metoprolol': 'D', 'propranolol': 'D',
    'furosemide': 'D', 'prednisone': 'D', 'dexamethasone': 'D',
    'levothyroxine': 'E', 'liothyronine': 'E', 'insulin': 'E',
    'heparin': 'E', 'warfarin': 'E', 'enoxaparin': 'E',
    'lidocaine': 'E', 'urokinase': 'E', 'choline': 'E',
}

print("="*60)
print("Consensus Label Training")
print("="*60)

# =============================================================================
# LOAD DATA
# =============================================================================

print("\nLoading data...")
df = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'results' / 'ecfp_model_results.csv')
print(f"Total compounds: {len(df)}")

# Add LiverTox labels
df['drug_name_lower'] = df['drug_name'].str.lower().str.strip()
df['livertox_cat'] = df['drug_name_lower'].map(LIVERTOX_CATEGORIES)
df['livertox_dili'] = df['livertox_cat'].map({
    'A': 1, 'B': 1, 'C': 1, 'D': 0, 'E': 0
})

# Find consensus
df['has_livertox'] = df['livertox_dili'].notna()
df['consensus'] = np.nan
df.loc[df['has_livertox'], 'consensus'] = (
    df.loc[df['has_livertox'], 'is_dili'] == df.loc[df['has_livertox'], 'livertox_dili']
)

# Disagreements
disagree = df[df['consensus'] == False]
print(f"\nDisagreements (to remove): {len(disagree)}")
for _, row in disagree.iterrows():
    lt = "DILI" if row['livertox_dili'] == 1 else "Safe"
    dr = "DILI" if row['is_dili'] == 1 else "Safe"
    print(f"  {row['drug_name']}: DILIrank={dr}, LiverTox={lt}({row['livertox_cat']})")

# =============================================================================
# CREATE CONSENSUS DATASET
# =============================================================================

print("\n" + "="*60)
print("CONSENSUS DATASET")
print("="*60)

# Option 1: Remove disagreements only
df_consensus = df[df['consensus'] != False].copy()
print(f"After removing disagreements: {len(df_consensus)}")

# Option 2: Use consensus labels where available
df_consensus['final_label'] = df_consensus['is_dili'].copy()
# Keep original DILIrank label (LiverTox is only for validation)

print(f"DILI+: {(df_consensus['final_label']==1).sum()}")
print(f"DILI-: {(df_consensus['final_label']==0).sum()}")

# =============================================================================
# TRAIN ON CONSENSUS DATA
# =============================================================================

print("\n" + "="*60)
print("TRAINING ON CONSENSUS LABELS")
print("="*60)

# Generate ECFP fingerprints
def smiles_to_ecfp(smiles, radius=2, nBits=1024):
    if not RDKIT_AVAILABLE or pd.isna(smiles):
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits)
            return np.array(fp)
    except:
        pass
    return None

# Generate fingerprints
print("Generating fingerprints...")
fingerprints = []
valid_idx = []

for idx, row in df_consensus.iterrows():
    if pd.notna(row.get('smiles')):
        fp = smiles_to_ecfp(row['smiles'])
        if fp is not None:
            fingerprints.append(fp)
            valid_idx.append(idx)

X_fp = np.array(fingerprints)
df_valid = df_consensus.loc[valid_idx].copy()
y = df_valid['final_label'].values

print(f"Valid compounds: {len(df_valid)}")
print(f"Fingerprint shape: {X_fp.shape}")

# Add tabular features
tabular_cols = ['logp', 'mw', 'hbd', 'hba', 'tpsa', 'k', 'ecnp_z']
tabular_cols = [c for c in tabular_cols if c in df_valid.columns]

X_tab = df_valid[tabular_cols].values
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
X_tab_scaled = scaler.fit_transform(X_tab)

# Combined features
X_combined = np.hstack([X_fp, X_tab_scaled])
print(f"Combined features: {X_combined.shape[1]}")

# Train model
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
model = RandomForestClassifier(n_estimators=100, max_depth=10, random_state=42, n_jobs=-1)

y_pred_consensus = cross_val_predict(model, X_combined, y, cv=cv, method='predict_proba')[:, 1]
auc_consensus = roc_auc_score(y, y_pred_consensus)

print(f"\n*** AUC on consensus data: {auc_consensus:.3f} ***")

# Compare to original
print(f"Previous AUC (all data): 0.884")
print(f"Change: {auc_consensus - 0.884:+.3f}")

# =============================================================================
# COMPARE ERROR RATES
# =============================================================================

print("\n" + "="*60)
print("ERROR ANALYSIS")
print("="*60)

df_valid['y_pred'] = y_pred_consensus
df_valid['pred_binary'] = (df_valid['y_pred'] > 0.5).astype(int)
df_valid['correct'] = df_valid['pred_binary'] == df_valid['final_label']
df_valid['fn'] = (df_valid['pred_binary'] == 0) & (df_valid['final_label'] == 1)
df_valid['fp'] = (df_valid['pred_binary'] == 1) & (df_valid['final_label'] == 0)

accuracy = df_valid['correct'].mean()
fn_count = df_valid['fn'].sum()
fp_count = df_valid['fp'].sum()

print(f"Accuracy: {accuracy:.1%}")
print(f"False negatives: {fn_count} (was 9)")
print(f"False positives: {fp_count} (was 29)")

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "="*60)
print("CONSENSUS TRAINING SUMMARY")
print("="*60)

print(f"""
Original dataset: 202 compounds
Removed (disagreements): {len(disagree)} compounds
  - Heparin, Levothyroxine, Ibuprofen, Enoxaparin
Consensus dataset: {len(df_valid)} compounds

Results:
  Original AUC: 0.884
  Consensus AUC: {auc_consensus:.3f}
  Change: {auc_consensus - 0.884:+.3f}

Error rates:
  False negatives: {fn_count} (was 9)
  False positives: {fp_count} (was 29)
  Accuracy: {accuracy:.1%}
""")

# Save
output = ROOT / 'research' / 'ecnp-generalization' / 'results' / 'consensus_model_results.csv'
df_valid.to_csv(output, index=False)
print(f"Saved: {output}")

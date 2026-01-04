"""
Data Integrity Revalidation
============================

Check for mismatches between the 901-drug (PubChem) and 202-drug (DrugBank) datasets.

1. Label consistency - same drugs have same labels?
2. SMILES validity - all SMILES valid and consistent?
3. Overlap analysis - how many drugs in common?
4. Class balance comparison
5. Feature distributions
"""
import pandas as pd
import numpy as np
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except:
    RDKIT_AVAILABLE = False

ROOT = Path(r'v:\new\h-perforatum-network-tox')

print("="*70)
print("DATA INTEGRITY REVALIDATION")
print("="*70)

# =============================================================================
# LOAD BOTH DATASETS
# =============================================================================

print("\n--- Loading Datasets ---")

# Tier 1: Full DILIrank from PubChem
tier1_path = ROOT / 'research' / 'ecnp-generalization' / 'results' / 'dilirank_full_smiles.csv'
df_tier1 = pd.read_csv(tier1_path)
print(f"Tier 1 (PubChem): {len(df_tier1)} drugs")

# Tier 2: Our original subset
tier2_path = ROOT / 'research' / 'ecnp-generalization' / 'results' / 'ecfp_model_results.csv'
df_tier2 = pd.read_csv(tier2_path)
print(f"Tier 2 (DrugBank): {len(df_tier2)} drugs")

# Original DILIrank
dilirank_path = ROOT / 'data' / 'external' / 'DILIrank_2.0.xlsx'
dilirank = pd.read_excel(dilirank_path, header=1)
print(f"Original DILIrank: {len(dilirank)} drugs")

# =============================================================================
# CHECK 1: COLUMN CONSISTENCY
# =============================================================================

print("\n" + "="*70)
print("CHECK 1: COLUMN STRUCTURE")
print("="*70)

print(f"\nTier 1 columns: {list(df_tier1.columns)}")
print(f"\nTier 2 columns (first 15): {list(df_tier2.columns)[:15]}...")

# =============================================================================
# CHECK 2: OVERLAP ANALYSIS
# =============================================================================

print("\n" + "="*70)
print("CHECK 2: OVERLAP ANALYSIS")
print("="*70)

# Normalize names for matching
def normalize(name):
    if pd.isna(name):
        return ''
    return str(name).lower().strip()

tier1_names = set(df_tier1['dilirank_name'].apply(normalize))
tier2_names = set(df_tier2['drug_name'].apply(normalize)) if 'drug_name' in df_tier2.columns else set()

# Also try dilirank_name in tier2
if 'dilirank_name' in df_tier2.columns:
    tier2_names_alt = set(df_tier2['dilirank_name'].apply(normalize))
else:
    tier2_names_alt = set()

overlap = tier1_names & tier2_names
overlap_alt = tier1_names & tier2_names_alt

print(f"Tier 1 unique names: {len(tier1_names)}")
print(f"Tier 2 unique names (drug_name): {len(tier2_names)}")
print(f"Tier 2 unique names (dilirank_name): {len(tier2_names_alt)}")
print(f"Overlap (drug_name): {len(overlap)}")
print(f"Overlap (dilirank_name): {len(overlap_alt)}")

# Find the overlapping drugs
overlap_list = []
for idx, row in df_tier2.iterrows():
    name_norm = normalize(row.get('dilirank_name') or row.get('drug_name'))
    match = df_tier1[df_tier1['dilirank_name'].apply(normalize) == name_norm]
    if len(match) > 0:
        overlap_list.append({
            'name': name_norm,
            'tier1_label': match['is_dili'].values[0],
            'tier2_label': row['is_dili'],
            'tier1_smiles': match['smiles'].values[0],
            'tier2_smiles': row['smiles']
        })

df_overlap = pd.DataFrame(overlap_list)
print(f"\nMatched drugs with both labels: {len(df_overlap)}")

# =============================================================================
# CHECK 3: LABEL CONSISTENCY
# =============================================================================

print("\n" + "="*70)
print("CHECK 3: LABEL CONSISTENCY")
print("="*70)

if len(df_overlap) > 0:
    label_match = (df_overlap['tier1_label'] == df_overlap['tier2_label']).sum()
    label_mismatch = (df_overlap['tier1_label'] != df_overlap['tier2_label']).sum()
    
    print(f"Labels match: {label_match}")
    print(f"Labels mismatch: {label_mismatch}")
    
    if label_mismatch > 0:
        print("\n[!] LABEL MISMATCHES:")
        mismatches = df_overlap[df_overlap['tier1_label'] != df_overlap['tier2_label']]
        for _, row in mismatches.iterrows():
            print(f"  {row['name']}: Tier1={row['tier1_label']}, Tier2={row['tier2_label']}")
else:
    print("[!] No overlapping drugs found for label comparison")

# =============================================================================
# CHECK 4: SMILES CONSISTENCY
# =============================================================================

print("\n" + "="*70)
print("CHECK 4: SMILES CONSISTENCY")
print("="*70)

if len(df_overlap) > 0 and RDKIT_AVAILABLE:
    # Compare canonical SMILES
    smiles_match = 0
    smiles_mismatch = 0
    
    for _, row in df_overlap.iterrows():
        try:
            mol1 = Chem.MolFromSmiles(row['tier1_smiles'])
            mol2 = Chem.MolFromSmiles(row['tier2_smiles'])
            if mol1 and mol2:
                can1 = Chem.MolToSmiles(mol1, canonical=True)
                can2 = Chem.MolToSmiles(mol2, canonical=True)
                if can1 == can2:
                    smiles_match += 1
                else:
                    smiles_mismatch += 1
            else:
                smiles_mismatch += 1
        except:
            smiles_mismatch += 1
    
    print(f"SMILES match (canonical): {smiles_match}")
    print(f"SMILES mismatch: {smiles_mismatch}")
else:
    print("Cannot check SMILES (no overlap or RDKit unavailable)")

# =============================================================================
# CHECK 5: CLASS BALANCE
# =============================================================================

print("\n" + "="*70)
print("CHECK 5: CLASS BALANCE COMPARISON")
print("="*70)

tier1_pos = df_tier1['is_dili'].sum()
tier1_neg = (df_tier1['is_dili'] == 0).sum()
tier1_ratio = tier1_pos / len(df_tier1) if len(df_tier1) > 0 else 0

tier2_pos = df_tier2['is_dili'].sum()
tier2_neg = (df_tier2['is_dili'] == 0).sum()
tier2_ratio = tier2_pos / len(df_tier2) if len(df_tier2) > 0 else 0

print(f"""
Tier 1 (901 drugs):
  DILI+: {tier1_pos} ({tier1_ratio:.1%})
  DILI-: {tier1_neg} ({1-tier1_ratio:.1%})

Tier 2 (202 drugs):
  DILI+: {tier2_pos} ({tier2_ratio:.1%})
  DILI-: {tier2_neg} ({1-tier2_ratio:.1%})
""")

# =============================================================================
# CHECK 6: ORIGINAL DILIRANK LABELS
# =============================================================================

print("\n" + "="*70)
print("CHECK 6: ORIGINAL DILIRANK LABEL DISTRIBUTION")
print("="*70)

print(dilirank['vDILI-Concern'].value_counts())

# Compare our binary conversion
print("\n[Verifying binary conversion logic:]")
print("DILI+ = vMost-DILI-concern OR vLess-DILI-concern")
print("DILI- = vNo-DILI-concern")
print("Excluded = Ambiguous-DILI-concern")

# =============================================================================
# CHECK 7: TIER 2 ORIGINAL SOURCE
# =============================================================================

print("\n" + "="*70)
print("CHECK 7: TIER 2 DATA SOURCE")
print("="*70)

# Check what columns indicate source
if 'dilirank' in df_tier2.columns:
    print("Tier 2 has 'dilirank' column with values:")
    print(df_tier2['dilirank'].value_counts())
    
    # Check is_dili vs dilirank correspondence
    print("\nChecking is_dili vs dilirank consistency:")
    for label in df_tier2['dilirank'].unique():
        subset = df_tier2[df_tier2['dilirank'] == label]
        dili_vals = subset['is_dili'].unique()
        print(f"  dilirank='{label}': is_dili values = {dili_vals}")

# =============================================================================
# CHECK 8: INVESTIGATE AUC DISCREPANCY
# =============================================================================

print("\n" + "="*70)
print("CHECK 8: POTENTIAL CAUSES OF AUC DISCREPANCY")
print("="*70)

print("""
Tier 1 AUC: 0.758 (901 drugs)
Tier 2 AUC: 0.870 (202 drugs)

Potential causes:
1. Label noise in larger dataset
2. Different drug populations (different mechanisms)
3. Selection bias in 202-drug subset
4. SMILES quality differences
5. Class balance differences
""")

# Check if Tier 2 is a curated subset
tier2_only = tier2_names_alt - tier1_names
print(f"\nDrugs in Tier 2 but NOT in Tier 1: {len(tier2_only)}")
if len(tier2_only) > 0:
    print("Examples:", list(tier2_only)[:10])

tier1_only = tier1_names - tier2_names_alt
print(f"\nDrugs in Tier 1 but NOT in Tier 2: {len(tier1_only)}")
if len(tier1_only) > 0:
    print("Examples:", list(tier1_only)[:10])

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "="*70)
print("INTEGRITY SUMMARY")
print("="*70)

print(f"""
Overlap: {len(df_overlap)} drugs in both datasets
Label matches: {label_match if len(df_overlap) > 0 else 'N/A'}
Label mismatches: {label_mismatch if len(df_overlap) > 0 else 'N/A'}
SMILES matches: {smiles_match if RDKIT_AVAILABLE and len(df_overlap) > 0 else 'N/A'}

Tier 1 class balance: {tier1_ratio:.1%} DILI+
Tier 2 class balance: {tier2_ratio:.1%} DILI+
""")

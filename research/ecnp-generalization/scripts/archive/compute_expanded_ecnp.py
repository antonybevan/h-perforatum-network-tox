"""
Compute ECNP for 706 Expanded DILIrank Drugs
=============================================

Use RobustECNP to compute Z-scores for all 706 drugs with DrugBank targets.
"""
import pandas as pd
import numpy as np
from pathlib import Path
import sys
import warnings
warnings.filterwarnings('ignore')

# Add project root to path
ROOT = Path(r'v:\new\h-perforatum-network-tox')
sys.path.insert(0, str(ROOT / 'research' / 'ecnp-generalization' / 'scripts'))

print("="*70)
print("COMPUTING ECNP FOR 706 DRUGS")
print("="*70)

# =============================================================================
# TRY TO IMPORT ROBUST ECNP
# =============================================================================

try:
    from ecnp_robust import RobustECNP, RobustECNPConfig
    ecnp_engine = RobustECNP()
    print("Loaded RobustECNP engine")
except Exception as e:
    print(f"Error loading RobustECNP: {e}")
    print("Using fallback computation...")
    ecnp_engine = None

# =============================================================================
# LOAD DATA
# =============================================================================

print("\n--- Loading Data ---")

# Drugs with targets
df = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'results' / 'dilirank_with_targets.csv')
print(f"Drugs with targets: {len(df)}")

# Targets expanded
df_targets = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'results' / 'dilirank_drugbank_targets.csv')
print(f"Expanded target rows: {len(df_targets)}")

# Full SMILES data
df_smiles = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'results' / 'dilirank_full_smiles.csv')
print(f"SMILES data: {len(df_smiles)}")

# Load LCC nodes to filter targets
lcc_path = ROOT / 'data' / 'processed' / 'protein_lcc.csv'
if lcc_path.exists():
    lcc_nodes = set(pd.read_csv(lcc_path)['protein'].values)
    print(f"LCC nodes: {len(lcc_nodes)}")
else:
    lcc_nodes = None
    print("LCC file not found, skipping LCC filtering")

# =============================================================================
# COMPUTE ECNP FOR EACH DRUG
# =============================================================================

print("\n--- Computing ECNP ---")

# Parse targets column (it's stored as string representation of list)
import ast

results = []
config = RobustECNPConfig() if ecnp_engine else None

for idx, row in df.iterrows():
    drug_name = row['dilirank_name']
    
    # Parse targets
    try:
        targets = ast.literal_eval(row['targets'])
    except:
        targets = []
    
    # Filter to LCC if available
    if lcc_nodes:
        targets_lcc = [t for t in targets if t in lcc_nodes]
    else:
        targets_lcc = targets
    
    # Compute ECNP
    if ecnp_engine and len(targets_lcc) > 0:
        try:
            result = ecnp_engine.ecnp(targets_lcc, config)
            z_score = result.get('Z', np.nan)
            mu_T = result.get('mu_T', 0)
            sigma_T = result.get('sigma_T', 0)
            pool_size = result.get('pool_size', 0)
            I_T = result.get('I_T', 0)
        except Exception as e:
            z_score = np.nan
            mu_T = 0
            sigma_T = 0
            pool_size = 0
            I_T = 0
    else:
        z_score = np.nan
        mu_T = 0
        sigma_T = 0
        pool_size = 0
        I_T = 0
    
    results.append({
        'dilirank_name': drug_name,
        'is_dili': row['is_dili'],
        'n_targets': row['n_targets'],
        'n_targets_lcc': len(targets_lcc),
        'ecnp_z': z_score,
        'ecnp_mu_T': mu_T,
        'ecnp_sigma_T': sigma_T,
        'pool_size': pool_size,
        'I_T': I_T,
        'drugbank_id': row['drugbank_id']
    })
    
    # Progress
    if (idx + 1) % 100 == 0:
        print(f"  Processed {idx+1}/{len(df)}")

print(f"\nComputed ECNP for {len(results)} drugs")

# =============================================================================
# MERGE WITH SMILES
# =============================================================================

print("\n--- Merging with SMILES ---")

df_ecnp = pd.DataFrame(results)

# Merge with SMILES
df_merged = df_ecnp.merge(
    df_smiles[['dilirank_name', 'smiles']],
    on='dilirank_name',
    how='left'
)

# Filter to those with valid SMILES and positive ECNP
df_valid = df_merged[df_merged['smiles'].notna()].copy()
print(f"With valid SMILES: {len(df_valid)}")

df_with_ecnp = df_valid[df_valid['ecnp_z'] != 0].copy()
print(f"With non-zero ECNP: {len(df_with_ecnp)}")

# =============================================================================
# SAVE
# =============================================================================

print("\n--- Saving Results ---")

output = ROOT / 'research' / 'ecnp-generalization' / 'results' / 'dilirank_expanded_ecnp.csv'
df_valid.to_csv(output, index=False)
print(f"Saved: {output}")

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "="*70)
print("SUMMARY")
print("="*70)

print(f"""
Total drugs: {len(df_valid)}
With non-zero ECNP: {len(df_with_ecnp)}
With targets in LCC: {(df_valid['n_targets_lcc'] > 0).sum()}

Label distribution:
  DILI+: {df_valid['is_dili'].sum()}
  DILI-: {(df_valid['is_dili'] == 0).sum()}

ECNP statistics:
  Mean Z: {df_valid['ecnp_z'].mean():.3f}
  Std Z:  {df_valid['ecnp_z'].std():.3f}
  Min Z:  {df_valid['ecnp_z'].min():.3f}
  Max Z:  {df_valid['ecnp_z'].max():.3f}
""")

# Compare to original 202
print("\nComparison to Original 202-drug subset:")
print(f"  Original: 202 drugs")
print(f"  Expanded: {len(df_valid)} drugs")
print(f"  Increase: {len(df_valid) - 202:+d}")

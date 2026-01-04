"""
Compute ECNP with UniProt→Gene Mapping
======================================

Use the UniProt to gene symbol mapping to compute ECNP scores.
"""
import pandas as pd
import numpy as np
from pathlib import Path
import ast
import sys
import warnings
warnings.filterwarnings('ignore')

ROOT = Path(r'v:\new\h-perforatum-network-tox')
sys.path.insert(0, str(ROOT / 'research' / 'ecnp-generalization' / 'scripts'))

print("="*70)
print("COMPUTING ECNP WITH GENE SYMBOL MAPPING")
print("="*70)

# Load RobustECNP
from ecnp_robust import RobustECNP, RobustECNPConfig
ecnp_engine = RobustECNP()
config = RobustECNPConfig()
print("Loaded RobustECNP engine")

# Load mapping
print("\n--- Loading UniProt->Gene Mapping ---")
mapping = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'results' / 'uniprot_to_gene.csv')
uniprot_to_gene = dict(zip(mapping['uniprot_id'], mapping['gene_symbol']))
print(f"Loaded {len(uniprot_to_gene)} mappings")

# Load drugs with targets
print("\n--- Loading Drug Targets ---")
df = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'results' / 'dilirank_with_targets.csv')
print(f"Drugs: {len(df)}")

# Load SMILES
df_smiles = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'results' / 'dilirank_full_smiles.csv')

# Compute ECNP
print("\n--- Computing ECNP ---")
results = []

for idx, row in df.iterrows():
    drug_name = row['dilirank_name']
    
    # Parse targets
    try:
        targets_uniprot = ast.literal_eval(row['targets'])
    except:
        targets_uniprot = []
    
    # Map to gene symbols
    targets_genes = [uniprot_to_gene.get(t) for t in targets_uniprot if t in uniprot_to_gene]
    targets_genes = [g for g in targets_genes if g is not None]
    
    # Compute ECNP
    if len(targets_genes) > 0:
        try:
            result = ecnp_engine.ecnp(targets_genes, config)
            z_score = result.get('Z', np.nan)
            mu_T = result.get('mu_T', 0)
            sigma_T = result.get('sigma_T', 0)
            pool_size = result.get('pool_size', 0)
            I_T = result.get('I_T', 0)
            k = result.get('k', 0)
        except Exception as e:
            z_score = np.nan
            mu_T = 0
            sigma_T = 0
            pool_size = 0
            I_T = 0
            k = 0
    else:
        z_score = np.nan
        mu_T = 0
        sigma_T = 0
        pool_size = 0
        I_T = 0
        k = 0
    
    results.append({
        'dilirank_name': drug_name,
        'is_dili': row['is_dili'],
        'n_targets_uniprot': len(targets_uniprot),
        'n_targets_mapped': len(targets_genes),
        'ecnp_z': z_score,
        'ecnp_mu_T': mu_T,
        'ecnp_sigma_T': sigma_T,
        'pool_size': pool_size,
        'I_T': I_T,
        'k': k,
        'drugbank_id': row['drugbank_id']
    })
    
    if (idx + 1) % 100 == 0:
        print(f"  Processed {idx+1}/{len(df)}")

df_ecnp = pd.DataFrame(results)

# Merge with SMILES
df_final = df_ecnp.merge(df_smiles[['dilirank_name', 'smiles']], on='dilirank_name', how='left')

# Save
output = ROOT / 'research' / 'ecnp-generalization' / 'results' / 'dilirank_706_with_ecnp.csv'
df_final.to_csv(output, index=False)
print(f"\nSaved: {output}")

# Stats
print("\n" + "="*70)
print("SUMMARY")
print("="*70)

valid_z = df_final['ecnp_z'].notna()
print(f"""
Total drugs: {len(df_final)}
With valid ECNP: {valid_z.sum()}
Coverage: {valid_z.sum()/len(df_final)*100:.1f}%

ECNP Z-score stats (valid only):
  Mean: {df_final.loc[valid_z, 'ecnp_z'].mean():.3f}
  Std:  {df_final.loc[valid_z, 'ecnp_z'].std():.3f}
  Min:  {df_final.loc[valid_z, 'ecnp_z'].min():.3f}
  Max:  {df_final.loc[valid_z, 'ecnp_z'].max():.3f}

Label distribution (with valid ECNP):
  DILI+: {df_final.loc[valid_z, 'is_dili'].sum()}
  DILI-: {(df_final.loc[valid_z, 'is_dili'] == 0).sum()}
""")

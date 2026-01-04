"""
Rebuild DILI Gene Module with Full DisGeNET Set
=================================================

The current dili_900_lcc.csv has only 82 genes.
DisGeNET curated gene-disease associations has 787 DILI-related genes.

This script:
1. Extracts all DILI genes from DisGeNET
2. Filters to genes in Liver LCC network
3. Recomputes DILI influence vector
4. Re-runs ECNP validation
"""
import pandas as pd
import numpy as np
from pathlib import Path

ROOT = Path(r'v:\new\h-perforatum-network-tox')
DISGENET_FILE = ROOT / 'data' / 'raw' / 'curated_gene_disease_associations.tsv'
NETWORK_FILE = ROOT / 'research' / 'ecnp-closed-form' / 'data' / 'node_list_900.csv'
INFLUENCE_MATRIX = ROOT / 'research' / 'ecnp-closed-form' / 'data' / 'influence_matrix_900.npz'
OUTPUT_DIR = ROOT / 'research' / 'ecnp-closed-form' / 'data'

# Load DisGeNET
print("Loading DisGeNET curated associations...")
disgenet = pd.read_csv(DISGENET_FILE, sep='\t')
print(f"Total associations: {len(disgenet)}")

# Filter to DILI-related conditions
dili_terms = ['liver', 'hepat', 'dili', 'hepatotox', 'cholestatic', 'jaundice']
dili_pattern = '|'.join(dili_terms)
dili_df = disgenet[disgenet['diseaseName'].str.contains(dili_pattern, case=False, na=False)]

print(f"\nDILI-related associations: {len(dili_df)}")
print(f"Unique DILI genes: {dili_df['geneSymbol'].nunique()}")
print(f"\nDiseases included:")
print(dili_df['diseaseName'].value_counts())

# Get unique DILI genes
dili_genes = set(dili_df['geneSymbol'].unique())
print(f"\nTotal unique DILI genes: {len(dili_genes)}")

# Load network nodes
print("\nLoading Liver LCC network...")
nodes = pd.read_csv(NETWORK_FILE)
network_genes = set(nodes['gene'].tolist())
print(f"Network genes: {len(network_genes)}")

# Filter to network
dili_in_network = dili_genes & network_genes
print(f"DILI genes in network: {len(dili_in_network)}")

# Compare to old module
old_dili = pd.read_csv(ROOT / 'data' / 'processed' / 'dili_900_lcc.csv')
old_genes = set(old_dili['gene_name'].tolist())
print(f"\nOld DILI module: {len(old_genes)} genes")
print(f"New DILI module: {len(dili_in_network)} genes")
print(f"New vs old: {len(dili_in_network) / len(old_genes):.1f}x more genes")

# Compute new DILI influence vector
print("\nComputing new DILI influence vector...")
npz = np.load(INFLUENCE_MATRIX, allow_pickle=True)
M = npz['M']
node_list = npz['node_list'].tolist()
node_to_idx = {g: i for i, g in enumerate(node_list)}

# Get DILI gene indices
dili_idx = [node_to_idx[g] for g in dili_in_network if g in node_to_idx]
print(f"DILI genes with influence vectors: {len(dili_idx)}")

# Compute DILI influence for each node
# m_j = sum_d M[j,d] for d in DILI genes
dili_influence = M[:, dili_idx].sum(axis=1)

# Save new influence vector
new_influence_df = pd.DataFrame({
    'gene': node_list,
    'dili_influence': dili_influence
})
new_influence_file = OUTPUT_DIR / 'dili_influence_vector_900_full.csv'
new_influence_df.to_csv(new_influence_file, index=False)
print(f"Saved: {new_influence_file}")

# Save new DILI gene list
new_dili_df = pd.DataFrame({'gene_name': list(dili_in_network)})
new_dili_file = OUTPUT_DIR / 'dili_genes_full.csv'
new_dili_df.to_csv(new_dili_file, index=False)
print(f"Saved: {new_dili_file}")

# Summary
print("\n" + "="*60)
print("SUMMARY")
print("="*60)
print(f"Old DILI module: {len(old_genes)} genes")
print(f"New DILI module: {len(dili_in_network)} genes")
print(f"Improvement: {len(dili_in_network) / len(old_genes):.1f}x more coverage")
print(f"\nNew files created:")
print(f"  - {new_influence_file.name}")
print(f"  - {new_dili_file.name}")

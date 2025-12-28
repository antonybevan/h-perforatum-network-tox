#!/usr/bin/env python3
"""
Create Unified Pre-Filtered Data Source

This script creates a single CSV with all pre-filtered data for consistency:
1. Liver expression (TPM ≥ 1)
2. Targets validated against both 700 and 900 networks
3. Network membership flags
4. DILI gene status

Output: data/processed/unified_analysis_data.csv
"""

import sys
from pathlib import Path
import pandas as pd
import networkx as nx

project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root / 'src'))

DATA_DIR = project_root / 'data'
GTEX_FILE = DATA_DIR / 'raw' / 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct'
MIN_TPM = 1.0

print("=" * 80)
print(" CREATING UNIFIED PRE-FILTERED DATA SOURCE")
print("=" * 80)

# 1. Load GTEx liver expression
print("\n[1/5] Loading GTEx liver expression...")
if not GTEX_FILE.exists():
    print(f"ERROR: GTEx file not found: {GTEX_FILE}")
    sys.exit(1)

df_gtex = pd.read_csv(GTEX_FILE, sep='\t', skiprows=2)
expression = {}
for _, row in df_gtex.iterrows():
    gene = row.get('Description', row.get('Name', ''))
    tpm = row.get('Liver', 0)
    if gene and pd.notna(tpm) and tpm >= MIN_TPM:
        expression[gene] = float(tpm)

print(f"  Genes with TPM ≥ {MIN_TPM}: {len(expression)}")

# 2. Load networks
print("\n[2/5] Loading networks...")
df700 = pd.read_parquet(DATA_DIR / 'processed' / 'network_700.parquet')
G700 = nx.from_pandas_edgelist(df700, 'gene1', 'gene2')
print(f"  Network 700: {G700.number_of_nodes()} nodes")

df900 = pd.read_parquet(DATA_DIR / 'processed' / 'network_900.parquet')
G900 = nx.from_pandas_edgelist(df900, 'protein1', 'protein2')
print(f"  Network 900: {G900.number_of_nodes()} nodes")

# Set of genes in both networks
genes_in_both = set(G700.nodes()) & set(G900.nodes())
print(f"  Genes in BOTH networks: {len(genes_in_both)}")

# 3. Load DILI genes
print("\n[3/5] Loading DILI genes...")
dili_700 = set(pd.read_csv(DATA_DIR / 'processed' / 'dili_700_lcc.csv')['protein_id'])
dili_900 = set(pd.read_csv(DATA_DIR / 'processed' / 'dili_900_lcc.csv')['gene_name'])
dili_all = dili_700 | dili_900
print(f"  DILI genes (700): {len(dili_700)}")
print(f"  DILI genes (900): {len(dili_900)}")

# 4. Load and validate targets
print("\n[4/5] Loading and validating targets...")
targets_df = pd.read_csv(DATA_DIR / 'processed' / 'targets.csv')

# Get unique targets per compound
hyp_targets = set(targets_df[targets_df['compound'] == 'Hyperforin']['gene_name'])
quer_targets = set(targets_df[targets_df['compound'] == 'Quercetin']['gene_name'])

# Filter to targets in both networks
hyp_valid = hyp_targets & genes_in_both
quer_valid = quer_targets & genes_in_both

print(f"  Hyperforin: {len(hyp_valid)}/{len(hyp_targets)} in both networks")
print(f"  Quercetin: {len(quer_valid)}/{len(quer_targets)} in both networks")

# 5. Create unified dataframe
print("\n[5/5] Creating unified data file...")

# Build gene-level data
data = []
all_genes = genes_in_both | hyp_targets | quer_targets

for gene in sorted(all_genes):
    row = {
        'gene_symbol': gene,
        'liver_tpm': expression.get(gene, 0.0),
        'expressed_liver': gene in expression,  # TPM >= 1
        'in_network_700': gene in G700,
        'in_network_900': gene in G900,
        'in_both_networks': gene in genes_in_both,
        'is_dili_gene': gene in dili_all,
        'hyperforin_target': gene in hyp_targets,
        'quercetin_target': gene in quer_targets,
        'hyperforin_target_valid': gene in hyp_valid,
        'quercetin_target_valid': gene in quer_valid,
    }
    data.append(row)

unified_df = pd.DataFrame(data)

# Save
output_file = DATA_DIR / 'processed' / 'unified_analysis_data.csv'
unified_df.to_csv(output_file, index=False)
print(f"\n✓ Saved: {output_file}")

# Summary
print("\n" + "=" * 80)
print(" SUMMARY")
print("=" * 80)
print(f"\nTotal genes in file: {len(unified_df)}")
print(f"Expressed in liver (TPM ≥ 1): {unified_df['expressed_liver'].sum()}")
print(f"In both networks: {unified_df['in_both_networks'].sum()}")
print(f"DILI genes: {unified_df['is_dili_gene'].sum()}")
print(f"\nVALID TARGETS (in both networks):")
print(f"  Hyperforin: {unified_df['hyperforin_target_valid'].sum()}")
print(f"  Quercetin: {unified_df['quercetin_target_valid'].sum()}")

# Create filtered targets CSV
print("\n" + "-" * 80)
print("Creating validated targets file...")

validated_targets = []
for gene in hyp_valid:
    validated_targets.append({
        'gene_symbol': gene,
        'compound': 'Hyperforin',
        'liver_tpm': expression.get(gene, 0.0),
        'in_network_700': True,
        'in_network_900': True
    })
for gene in quer_valid:
    validated_targets.append({
        'gene_symbol': gene,
        'compound': 'Quercetin', 
        'liver_tpm': expression.get(gene, 0.0),
        'in_network_700': True,
        'in_network_900': True
    })

validated_df = pd.DataFrame(validated_targets)
validated_file = DATA_DIR / 'processed' / 'targets_validated.csv'
validated_df.to_csv(validated_file, index=False)
print(f"✓ Saved: {validated_file}")

# Create filtered expression CSV for liver
print("\n" + "-" * 80)
print("Creating liver expression file (TPM ≥ 1)...")

liver_expr = [{'gene_symbol': g, 'liver_tpm': t} for g, t in expression.items()]
liver_df = pd.DataFrame(liver_expr).sort_values('liver_tpm', ascending=False)
liver_file = DATA_DIR / 'processed' / 'liver_expression_filtered.csv'
liver_df.to_csv(liver_file, index=False)
print(f"✓ Saved: {liver_file}")
print(f"  {len(liver_df)} genes with TPM ≥ {MIN_TPM}")

print("\n" + "=" * 80)
print(" DONE - USE THESE FILES FOR ALL ANALYSES")
print("=" * 80)
print(f"\n1. {output_file}")
print(f"   → Complete gene-level data with all flags")
print(f"\n2. {validated_file}")
print(f"   → Validated targets in both networks")
print(f"\n3. {liver_file}")
print(f"   → Pre-filtered liver expression (TPM ≥ 1)")

#!/usr/bin/env python3
"""
Create Consistent LCC-Filtered Data Files

Updates:
1. liver_proteome.csv - Liver-expressed genes (TPM >= 1)
2. targets_lcc.csv - Targets filtered to liver LCC (9 Hyperforin, 62 Quercetin)
3. network_liver_lcc.parquet - Liver LCC subnetwork

This ensures consistency across all analyses.
"""

import sys
from pathlib import Path
import pandas as pd
import networkx as nx

project_root = Path(__file__).resolve().parent.parent
DATA_DIR = project_root / 'data'
GTEX_FILE = DATA_DIR / 'raw' / 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct'
MIN_TPM = 1.0

print("=" * 80)
print(" CREATING CONSISTENT LCC-FILTERED DATA FILES")
print("=" * 80)

# 1. Load and filter GTEx liver expression
print("\n[1/6] Loading GTEx liver expression (TPM ≥ 1)...")
df_gtex = pd.read_csv(GTEX_FILE, sep='\t', skiprows=2)

liver_genes = {}
for _, row in df_gtex.iterrows():
    gene = row.get('Description', row.get('Name', ''))
    tpm = row.get('Liver', 0)
    if gene and pd.notna(tpm) and tpm >= MIN_TPM:
        liver_genes[gene] = float(tpm)

print(f"  Genes with liver TPM ≥ {MIN_TPM}: {len(liver_genes)}")

# Save liver_proteome.csv
liver_df = pd.DataFrame([
    {'gene_symbol': g, 'liver_tpm': t} 
    for g, t in liver_genes.items()
]).sort_values('liver_tpm', ascending=False)

liver_file = DATA_DIR / 'processed' / 'liver_proteome.csv'
liver_df.to_csv(liver_file, index=False)
print(f"  ✓ Saved: {liver_file}")

# 2. Load STRING networks
print("\n[2/6] Loading STRING networks...")
df700 = pd.read_parquet(DATA_DIR / 'processed' / 'network_700.parquet')
G700_full = nx.from_pandas_edgelist(df700, 'gene1', 'gene2')
print(f"  Network 700: {G700_full.number_of_nodes()} nodes, {G700_full.number_of_edges()} edges")

df900 = pd.read_parquet(DATA_DIR / 'processed' / 'network_900.parquet')
G900_full = nx.from_pandas_edgelist(df900, 'protein1', 'protein2')
print(f"  Network 900: {G900_full.number_of_nodes()} nodes, {G900_full.number_of_edges()} edges")

# 3. Create liver-filtered LCC for each network
print("\n[3/6] Creating liver-filtered LCC networks...")

def get_liver_lcc(G, liver_genes):
    """Filter network to liver genes and extract LCC."""
    liver_nodes = [n for n in G.nodes() if n in liver_genes]
    G_liver = G.subgraph(liver_nodes).copy()
    if len(G_liver) == 0:
        return set()
    lcc = max(nx.connected_components(G_liver), key=len)
    return lcc

lcc_700 = get_liver_lcc(G700_full, liver_genes)
lcc_900 = get_liver_lcc(G900_full, liver_genes)

print(f"  Liver LCC (700): {len(lcc_700)} nodes")
print(f"  Liver LCC (900): {len(lcc_900)} nodes")

# Get intersection (nodes in both LCCs)
lcc_both = lcc_700 & lcc_900
print(f"  In BOTH LCCs: {len(lcc_both)} nodes")

# 4. Save liver LCC networks
print("\n[4/6] Saving liver LCC networks...")

# For network 700
edges_700 = [(u, v) for u, v in G700_full.edges() if u in lcc_700 and v in lcc_700]
lcc_700_df = pd.DataFrame(edges_700, columns=['gene1', 'gene2'])
lcc_700_file = DATA_DIR / 'processed' / 'network_700_liver_lcc.parquet'
lcc_700_df.to_parquet(lcc_700_file)
print(f"  ✓ Saved: {lcc_700_file} ({len(edges_700)} edges)")

# For network 900
edges_900 = [(u, v) for u, v in G900_full.edges() if u in lcc_900 and v in lcc_900]
lcc_900_df = pd.DataFrame(edges_900, columns=['gene1', 'gene2'])
lcc_900_file = DATA_DIR / 'processed' / 'network_900_liver_lcc.parquet'
lcc_900_df.to_parquet(lcc_900_file)
print(f"  ✓ Saved: {lcc_900_file} ({len(edges_900)} edges)")

# 5. Load and filter targets to LCC
print("\n[5/6] Filtering targets to liver LCC...")
targets_df = pd.read_csv(DATA_DIR / 'processed' / 'targets.csv')

hyp_all = set(targets_df[targets_df['compound'] == 'Hyperforin']['gene_name'])
quer_all = set(targets_df[targets_df['compound'] == 'Quercetin']['gene_name'])

# Filter to liver LCC (use intersection of both networks for consistency)
hyp_lcc = hyp_all & lcc_both
quer_lcc = quer_all & lcc_both

print(f"  Hyperforin: {len(hyp_lcc)}/{len(hyp_all)} in liver LCC")
print(f"  Quercetin: {len(quer_lcc)}/{len(quer_all)} in liver LCC")

# Create targets_lcc.csv
lcc_targets = []
for gene in sorted(hyp_lcc):
    lcc_targets.append({
        'gene_symbol': gene,
        'compound': 'Hyperforin',
        'liver_tpm': liver_genes.get(gene, 0.0),
        'in_lcc_700': gene in lcc_700,
        'in_lcc_900': gene in lcc_900
    })

for gene in sorted(quer_lcc):
    lcc_targets.append({
        'gene_symbol': gene,
        'compound': 'Quercetin',
        'liver_tpm': liver_genes.get(gene, 0.0),
        'in_lcc_700': gene in lcc_700,
        'in_lcc_900': gene in lcc_900
    })

targets_lcc_df = pd.DataFrame(lcc_targets)
targets_lcc_file = DATA_DIR / 'processed' / 'targets_lcc.csv'
targets_lcc_df.to_csv(targets_lcc_file, index=False)
print(f"  ✓ Saved: {targets_lcc_file}")

# 6. Summary
print("\n" + "=" * 80)
print(" SUMMARY - CONSISTENT DATA FILES")
print("=" * 80)

print(f"""
Files Created/Updated:
  1. {liver_file}
     → {len(liver_df)} genes with liver TPM ≥ {MIN_TPM}

  2. {lcc_700_file}
     → Liver LCC network (700): {len(lcc_700)} nodes, {len(edges_700)} edges

  3. {lcc_900_file}
     → Liver LCC network (900): {len(lcc_900)} nodes, {len(edges_900)} edges

  4. {targets_lcc_file}
     → Hyperforin: {len(hyp_lcc)} targets in liver LCC
     → Quercetin: {len(quer_lcc)} targets in liver LCC

FINAL TARGET COUNTS FOR ALL ANALYSES:
  ┌─────────────┬─────────┬───────────┐
  │ Compound    │ Total   │ In LCC    │
  ├─────────────┼─────────┼───────────┤
  │ Hyperforin  │ {len(hyp_all):>7} │ {len(hyp_lcc):>9} │
  │ Quercetin   │ {len(quer_all):>7} │ {len(quer_lcc):>9} │
  └─────────────┴─────────┴───────────┘

USE THESE FILES FOR ALL FUTURE ANALYSES!
""")

# Show which Hyperforin targets were excluded
excluded_hyp = hyp_all - hyp_lcc
if excluded_hyp:
    print(f"Hyperforin targets excluded (not in liver LCC): {sorted(excluded_hyp)}")

excluded_quer = quer_all - quer_lcc
if excluded_quer:
    print(f"\nQuercetin targets excluded: {len(excluded_quer)} targets")
    print(f"  (run with --verbose to see full list)")

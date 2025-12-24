"""
Script to regenerate DILI files from correct source data.
Extracts Drug-Induced Liver Injury genes from curated_gene_disease_associations.tsv
"""

import pandas as pd
import networkx as nx
import sys
from pathlib import Path

sys.path.append('src')
from network_tox.utils.data_loader import load_liver_genes
from network_tox.core.network import filter_to_tissue

DATA_DIR = Path('data')

print("="*80)
print("REGENERATING CORRECT DILI FILES")
print("="*80)

# Step 1: Extract DILI genes from curated source
print("\n1. Extracting DILI genes from curated_gene_disease_associations.tsv...")
curated = pd.read_csv(DATA_DIR / 'raw/curated_gene_disease_associations.tsv', sep='\t')

# Filter for Drug-Induced Liver Injury
dili = curated[curated['diseaseName'] == 'Drug-Induced Liver Injury'].copy()
print(f"   Found {len(dili)} DILI gene associations")

# Get unique genes
dili_genes_df = dili[['geneSymbol', 'geneId', 'score']].drop_duplicates(subset=['geneSymbol'])
dili_genes_df = dili_genes_df.rename(columns={'geneSymbol': 'gene_name'})
dili_genes_df = dili_genes_df.sort_values('score', ascending=False)
print(f"   Unique DILI genes: {len(dili_genes_df)}")
print(f"   Top genes: {list(dili_genes_df['gene_name'].head(10))}")

# Save raw DILI genes
raw_path = DATA_DIR / 'raw/dili_genes_raw.csv'
dili_genes_df.to_csv(raw_path, index=False)
print(f"   Saved to {raw_path}")

# Step 2: Load GTEx liver genes
print("\n2. Loading GTEx liver-expressed genes...")
gtex_path = DATA_DIR / 'raw/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct'
liver_genes = load_liver_genes(gtex_path)
print(f"   Liver-expressed genes: {len(liver_genes)}")

# Step 3: Filter for each network threshold
for threshold in ['900', '700']:
    print(f"\n3. Processing {threshold} threshold network...")
    
    # Load network
    net_path = DATA_DIR / f'processed/network_{threshold}.parquet'
    df_net = pd.read_parquet(net_path)
    
    if 'protein1' in df_net.columns:
        G = nx.from_pandas_edgelist(df_net, 'protein1', 'protein2')
    else:
        G = nx.from_pandas_edgelist(df_net, 'gene1', 'gene2')
    
    print(f"   Network nodes: {G.number_of_nodes()}")
    
    # Filter to liver
    G_liver = filter_to_tissue(G, liver_genes)
    print(f"   Liver network nodes: {G_liver.number_of_nodes()}")
    
    # Filter DILI genes to LCC
    dili_in_lcc = [g for g in dili_genes_df['gene_name'] if g in G_liver]
    print(f"   DILI genes in liver LCC: {len(dili_in_lcc)}")
    
    # Create LCC-filtered DILI file
    dili_lcc_df = dili_genes_df[dili_genes_df['gene_name'].isin(dili_in_lcc)].copy()
    
    # Add columns to match expected format
    dili_lcc_df['protein_id'] = dili_lcc_df['gene_name']  # For compatibility with validation scripts
    dili_lcc_df['diseaseName'] = 'Drug-Induced Liver Injury'
    dili_lcc_df['diseaseId'] = 'umls:C0860207'
    
    # Save
    lcc_path = DATA_DIR / f'processed/dili_{threshold}_lcc.csv'
    dili_lcc_df.to_csv(lcc_path, index=False)
    print(f"   Saved to {lcc_path}")

print("\n" + "="*80)
print("DILI FILES REGENERATED SUCCESSFULLY")
print("="*80)

# Verify
print("\nVerification:")
for threshold in ['900', '700']:
    df = pd.read_csv(DATA_DIR / f'processed/dili_{threshold}_lcc.csv')
    print(f"  {threshold}: {len(df)} genes, disease={df['diseaseName'].unique()[0]}")

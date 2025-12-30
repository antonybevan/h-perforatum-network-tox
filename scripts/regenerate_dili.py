#!/usr/bin/env python3
"""
Regenerate DILI genes from DisGeNET source data.

Extracts Drug-Induced Liver Injury genes from curated_gene_disease_associations.tsv
and creates raw and LCC-filtered versions.

Usage:
    python scripts/regenerate_dili.py
"""

import pandas as pd
import networkx as nx
from pathlib import Path

DATA_DIR = Path('data')


def main():
    print("=" * 70)
    print(" REGENERATE DILI GENES FROM DISGENET SOURCE")
    print("=" * 70)
    print()
    
    # Step 1: Load curated associations from DisGeNET
    source_file = DATA_DIR / 'raw' / 'curated_gene_disease_associations.tsv'
    
    if not source_file.exists():
        print(f"ERROR: Source file not found: {source_file}")
        print("Download from: https://www.disgenet.org/downloads")
        return
    
    print(f"[1] Loading DisGeNET curated associations...")
    print(f"    Source: {source_file}")
    
    curated = pd.read_csv(source_file, sep='\t')
    print(f"    Total associations: {len(curated)}")
    
    # Step 2: Filter for Drug-Induced Liver Injury
    print()
    print("[2] Filtering for Drug-Induced Liver Injury (DILI)...")
    
    dili = curated[curated['diseaseName'] == 'Drug-Induced Liver Injury'].copy()
    print(f"    DILI associations: {len(dili)}")
    
    # Step 3: Create unique gene list
    print()
    print("[3] Creating unique gene list...")
    
    dili_genes_df = dili[['geneSymbol', 'geneId', 'score', 'diseaseName', 'diseaseId']].copy()
    dili_genes_df = dili_genes_df.rename(columns={'geneSymbol': 'gene_name'})
    dili_genes_df = dili_genes_df.drop_duplicates(subset=['gene_name'])
    dili_genes_df = dili_genes_df.sort_values('score', ascending=False)
    
    print(f"    Unique DILI genes: {len(dili_genes_df)}")
    print(f"    Top 10 genes: {list(dili_genes_df['gene_name'].head(10))}")
    
    # Step 4: Save raw DILI genes
    raw_path = DATA_DIR / 'raw' / 'dili_genes_raw.csv'
    dili_genes_df.to_csv(raw_path, index=False)
    print()
    print(f"[4] Saved raw DILI genes: {raw_path}")
    print(f"    Genes: {len(dili_genes_df)}")
    
    # Step 5: Create LCC-filtered versions for each network threshold
    print()
    print("[5] Creating LCC-filtered versions...")
    
    # Load liver proteome
    liver_df = pd.read_csv(DATA_DIR / 'processed' / 'liver_proteome.csv')
    liver_genes = set(liver_df['gene_symbol'])
    print(f"    Liver proteome genes: {len(liver_genes)}")
    
    for threshold in [700, 900]:
        print()
        print(f"    Processing ≥{threshold} network...")
        
        # Load LCC network
        net_file = DATA_DIR / 'processed' / f'network_{threshold}_liver_lcc.parquet'
        
        if not net_file.exists():
            print(f"    WARNING: {net_file} not found, skipping")
            continue
        
        net_df = pd.read_parquet(net_file)
        
        if 'gene1' in net_df.columns:
            G = nx.from_pandas_edgelist(net_df, 'gene1', 'gene2')
        else:
            G = nx.from_pandas_edgelist(net_df, 'protein1', 'protein2')
        
        lcc_nodes = set(G.nodes())
        print(f"    LCC nodes: {len(lcc_nodes)}")
        
        # Filter DILI genes to LCC
        dili_lcc = dili_genes_df[dili_genes_df['gene_name'].isin(lcc_nodes)].copy()
        
        # Add protein_id column for compatibility
        dili_lcc['protein_id'] = dili_lcc['gene_name']
        
        lcc_path = DATA_DIR / 'processed' / f'dili_{threshold}_lcc.csv'
        dili_lcc.to_csv(lcc_path, index=False)
        
        print(f"    Saved: {lcc_path}")
        print(f"    DILI genes in LCC: {len(dili_lcc)}/{len(dili_genes_df)}")
    
    # Summary
    print()
    print("=" * 70)
    print(" SUMMARY")
    print("=" * 70)
    print()
    print(f"Source: curated_gene_disease_associations.tsv (DisGeNET)")
    print(f"Disease: Drug-Induced Liver Injury")
    print()
    print("Generated files:")
    print(f"  - data/raw/dili_genes_raw.csv ({len(dili_genes_df)} genes)")
    print(f"  - data/processed/dili_700_lcc.csv")
    print(f"  - data/processed/dili_900_lcc.csv")
    print()
    print("Done!")


if __name__ == '__main__':
    main()

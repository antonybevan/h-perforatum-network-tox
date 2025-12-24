"""
Comprehensive verification of network and target files.
"""

import pandas as pd
import networkx as nx
import sys
sys.path.append('src')
from pathlib import Path
from network_tox.core.network import filter_to_tissue

DATA_DIR = Path('data')

print("="*80)
print("NETWORK AND TARGET FILES VERIFICATION")
print("="*80)

# 1. NETWORK PARQUET FILES
print("\n1. NETWORK PARQUET FILES")
print("-" * 80)

for threshold in ['900', '700']:
    net_path = DATA_DIR / f'processed/network_{threshold}.parquet'
    df = pd.read_parquet(net_path)
    
    # Get node columns
    col1, col2 = df.columns[0], df.columns[1]
    nodes = set(df[col1]) | set(df[col2])
    
    # Create graph
    G = nx.from_pandas_edgelist(df, col1, col2)
    
    print(f"\nNetwork {threshold}:")
    print(f"  File: {net_path.name}")
    print(f"  Columns: {list(df.columns)}")
    print(f"  Edges: {len(df):,}")
    print(f"  Unique nodes (from df): {len(nodes):,}")
    print(f"  Nodes in graph: {G.number_of_nodes():,}")
    print(f"  Graph edges: {G.number_of_edges():,}")
    
    # Check liver filtering status
    liver_df = pd.read_csv(DATA_DIR / 'processed/liver_proteome.csv')
    liver_genes = set(liver_df['gene_symbol'])
    G_liver = filter_to_tissue(G, liver_genes)
    
    print(f"  After liver filter: {G_liver.number_of_nodes():,} nodes")
    print(f"  Expected from validation: {7677 if threshold=='900' else 9773}")
    print(f"  MATCH: {G_liver.number_of_nodes() == (7677 if threshold=='900' else 9773)}")

# 2. TARGET CSV FILES
print("\n2. TARGET CSV FILES")
print("-" * 80)

for threshold in ['900', '700']:
    tgt_path = DATA_DIR / f'processed/targets_{threshold}.csv'
    df = pd.read_csv(tgt_path)
    
    print(f"\nTargets {threshold}:")
    print(f"  File: {tgt_path.name}")
    print(f"  Columns: {list(df.columns)}")
    print(f"  Total rows: {len(df)}")
    
    # Hyperforin
    hf = df[df['compound'] == 'Hyperforin']
    print(f"\n  Hyperforin targets: {len(hf)}")
    print(f"    NR1I2 present: {'NR1I2' in hf['gene_name'].values} [{'OK' if 'NR1I2' in hf['gene_name'].values else 'FAIL'}]")
    if 'NR1I2' in hf['gene_name'].values:
        nr1i2_row = hf[hf['gene_name'] == 'NR1I2']
        print(f"    NR1I2 protein ID: {nr1i2_row['protein_id'].values[0]}")
        print(f"    NR1I2 source: {nr1i2_row['source'].values[0]}")
    
    print(f"    Gene names: {list(hf['gene_name'].head(5))}...")
    
    # Quercetin
    qu = df[df['compound'] == 'Quercetin']
    print(f"\n  Quercetin targets: {len(qu)}")
    print(f"    Gene names: {list(qu['gene_name'].head(5))}...")
    
    # Liver network validation
    net_df = pd.read_parquet(DATA_DIR / f'processed/network_{threshold}.parquet')
    col1, col2 = net_df.columns[0], net_df.columns[1]
    G = nx.from_pandas_edgelist(net_df, col1, col2)
    liver_df = pd.read_csv(DATA_DIR / 'processed/liver_proteome.csv')
    liver_genes = set(liver_df['gene_symbol'])
    G_liver = filter_to_tissue(G, liver_genes)
    
    hf_in_liver = [g for g in hf['gene_name'] if g in G_liver]
    qu_in_liver = [g for g in qu['gene_name'] if g in G_liver]
    
    print(f"\n  In liver network:")
    print(f"    Hyperforin: {len(hf)} total -> {len(hf_in_liver)} in liver")
    print(f"    Quercetin: {len(qu)} total -> {len(qu_in_liver)} in liver")
    print(f"    Expected from validation: Hyperforin=9, Quercetin=63")
    print(f"    MATCH: {len(hf_in_liver)==9 and len(qu_in_liver)==63}")

print("\n" + "="*80)
print("VERIFICATION COMPLETE")
print("="*80)

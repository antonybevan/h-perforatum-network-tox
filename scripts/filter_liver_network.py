"""
Filter network to liver-expressed genes and extract LCC.

Usage:
    python scripts/filter_liver_network.py --input data/processed/network_900_raw.parquet --output data/processed/network_900.parquet
"""

import pandas as pd
import networkx as nx
import argparse
import sys

sys.path.append('src')
from network_tox.core.network import filter_to_tissue

def filter_network(input_file, liver_genes_file, output_file):
    """Filter network to liver genes and extract LCC."""
    
    # Load network
    print(f"Loading network from {input_file}...")
    df = pd.read_parquet(input_file)
    print(f"  Edges: {len(df):,}")
    
    # Get unique genes
    col1, col2 = df.columns[0], df.columns[1]
    all_genes = set(df[col1]) | set(df[col2])
    print(f"  Genes: {len(all_genes):,}")
    
    # Load liver genes
    print(f"\nLoading liver genes from {liver_genes_file}...")
    liver_df = pd.read_csv(liver_genes_file)
    liver_genes = set(liver_df['gene_symbol'])
    print(f"  Liver genes: {len(liver_genes):,}")
    
    # Create graph
    print("\nFiltering to liver-expressed genes...")
    G = nx.from_pandas_edgelist(df, col1, col2)
    G_liver = filter_to_tissue(G, liver_genes)
    
    print(f"  Nodes after liver filter: {G_liver.number_of_nodes():,}")
    print(f"  Edges after liver filter: {G_liver.number_of_edges():,}")
    
    # Convert back to dataframe
    edges = []
    for u, v in G_liver.edges():
        edges.append({col1: u, col2: v})
    
    df_filtered = pd.DataFrame(edges)
    
    # Save
    df_filtered.to_parquet(output_file, index=False)
    print(f"\n✓ Saved filtered network to {output_file}")
    
    return df_filtered

def main():
    parser = argparse.ArgumentParser(description='Filter network to liver genes')
    parser.add_argument('--input', type=str, required=True,
                       help='Input network parquet file')
    parser.add_argument('--output', type=str, required=True,
                       help='Output filtered network parquet file')
    parser.add_argument('--liver', type=str, default=str(Path('data') / 'processed' / 'liver_proteome.csv'),
                       help='Liver genes CSV file')
    
    args = parser.parse_args()
    
    filter_network(args.input, args.liver, args.output)
    print("\n✓ Liver filtering complete!")

if __name__ == '__main__':
    main()

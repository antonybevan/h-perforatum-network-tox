import sys
import networkx as nx
from pathlib import Path

# Add src to path
sys.path.append('src')

from network_tox.utils.data_loader import load_liver_genes
from network_tox.core.network import construct_liver_network

DATA_DIR = Path('data')

def main():
    print("Generating network with threshold 700...")
    gtex_path = DATA_DIR / 'raw/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct'

    if not gtex_path.exists():
        print(f"Error: {gtex_path} not found.")
        return

    liver_genes = load_liver_genes(gtex_path)
    print(f"Loaded {len(liver_genes)} liver genes.")

    # Check string file existence
    string_path = DATA_DIR / 'raw/protein.links.v11.5.txt.gz'
    if not string_path.exists():
        print(f"Error: {string_path} not found.")
        return

    G = construct_liver_network(liver_genes, score_threshold=700)
    print(f"Network generated: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    out_path = DATA_DIR / 'processed/network_700.parquet'
    nx.to_pandas_edgelist(G).to_parquet(out_path)
    print(f"Saved to {out_path}")

if __name__ == "__main__":
    main()

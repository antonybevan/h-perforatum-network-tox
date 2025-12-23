"""
Final Liver Analysis Script for H. perforatum Network Toxicology.
"""

import pandas as pd
import networkx as nx
from pathlib import Path
import sys

# Add src to path
sys.path.append('src')

from network_tox.utils.data_loader import load_liver_genes
from network_tox.core.network import filter_to_tissue
from network_tox.analysis.rwr import run_rwr

DATA_DIR = Path('data')
RESULTS_DIR = Path('results')

def get_hyperforin_targets():
    """Load and map Hyperforin targets."""
    print("Loading Hyperforin targets...")
    targets_df = pd.read_csv(DATA_DIR / 'raw/targets_raw.csv')
    mapping_df = pd.read_csv(DATA_DIR / 'external/uniprot_mapping.csv', header=None, comment='#', names=['protein_id', 'gene_name'])

    # Filter for Hyperforin
    hyp_targets = targets_df[targets_df['compound'] == 'Hyperforin']

    # Map to gene names
    # First, creating a dictionary for mapping
    mapping_dict = dict(zip(mapping_df['protein_id'], mapping_df['gene_name']))

    genes = []
    for pid in hyp_targets['protein_id']:
        if pid in mapping_dict:
            genes.append(mapping_dict[pid])
        else:
            # Fallback if the protein_id looks like a gene name (unlikely given the data, but good practice)
            # In this dataset, P08684 is mapped to CYP3A4, etc.
            pass

    print(f"Found {len(genes)} mapped targets for Hyperforin.")
    return list(set(genes))

def load_network(liver_genes):
    """Load network or generate synthetic if missing."""
    net_path = DATA_DIR / 'processed/network_900.parquet'

    if net_path.exists():
        print(f"Loading network from {net_path}...")
        # Assuming parquet contains edgelist
        df = pd.read_parquet(net_path)
        G = nx.from_pandas_edgelist(df, 'gene1', 'gene2')
    else:
        print(f"WARNING: {net_path} not found.")
        print("Generating SYNTHETIC network for demonstration purposes based on liver proteome.")

        # Create a synthetic scale-free graph
        # Use liver genes as nodes
        nodes = list(liver_genes)
        # Taking a subset to ensure reasonable execution time if list is huge,
        # but liver_proteome.csv is ~13k genes.
        # Barabasi-Albert generates edges based on attachment.
        # We'll map integer nodes to gene names.

        n_nodes = len(nodes)
        if n_nodes < 10:
            raise ValueError("Not enough liver genes to build network.")

        print(f"Building synthetic graph with {n_nodes} nodes...")
        # m=5 edges per new node
        G_synthetic = nx.barabasi_albert_graph(n_nodes, 5, seed=42)

        # Relabel nodes
        mapping = {i: nodes[i] for i in range(n_nodes)}
        G = nx.relabel_nodes(G_synthetic, mapping)

    return G

def main():
    print("--- Starting H. perforatum Liver Analysis ---")

    # 1. Load Liver Genes
    print("Loading liver genes...")
    # Use the raw GCT file as source of truth
    gtex_path = DATA_DIR / 'raw/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct'
    liver_genes = load_liver_genes(gtex_path)
    print(f"Loaded {len(liver_genes)} liver-expressed genes.")

    # 2. Load and Filter Network
    G = load_network(liver_genes)

    print("Filtering network to liver context...")
    G_liver = filter_to_tissue(G, liver_genes)
    print(f"Liver Network: {G_liver.number_of_nodes()} nodes, {G_liver.number_of_edges()} edges")

    # 3. Get Targets
    targets = get_hyperforin_targets()
    valid_targets = [t for t in targets if t in G_liver]
    print(f"Targets in Liver Network: {len(valid_targets)}/{len(targets)}")

    if not valid_targets:
        print("Error: No targets found in the network. Cannot run RWR.")
        return

    # 4. Run RWR
    print("Running Random Walk with Restart...")
    scores = run_rwr(G_liver, valid_targets, restart_prob=0.7)

    # 5. Save Results
    print("Saving results...")
    results_df = pd.DataFrame(list(scores.items()), columns=['gene', 'rwr_score'])
    results_df = results_df.sort_values('rwr_score', ascending=False)

    output_path = RESULTS_DIR / 'liver_analysis_output.csv'
    results_df.to_csv(output_path, index=False)
    print(f"Results saved to {output_path}")

    # 6. Show Findings
    print("\n--- TOP 10 INFLUENTIAL GENES ---")
    print(results_df.head(10).to_string(index=False))

    # Verify Hyperforin Targets ranking
    print("\n--- HYPERFORIN TARGETS RANKING ---")
    target_ranks = results_df[results_df['gene'].isin(valid_targets)]
    print(target_ranks.to_string(index=False))

if __name__ == "__main__":
    main()

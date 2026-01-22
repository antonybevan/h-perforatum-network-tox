"""Core network operations."""

import networkx as nx
import pandas as pd
import gzip


def load_string_network(threshold, links_file, info_file):
    """
    Load STRING network at specified confidence threshold.
    
    Args:
        threshold: Minimum combined score
        links_file: Path to string_links.txt.gz
        info_file: Path to string_info.txt.gz
        
    Returns:
        NetworkX graph (LCC)
    """
    with gzip.open(info_file, 'rt') as f:
        df_info = pd.read_csv(f, sep='\t')
    id_map = dict(zip(df_info['#string_protein_id'], df_info['preferred_name']))
    
    with gzip.open(links_file, 'rt') as f:
        df_links = pd.read_csv(f, sep=' ')
    
    df = df_links[df_links['combined_score'] >= threshold].copy()
    df['gene1'] = df['protein1'].map(id_map)
    df['gene2'] = df['protein2'].map(id_map)
    df = df.dropna(subset=['gene1', 'gene2'])
    
    G = nx.Graph()
    G.add_edges_from(zip(df['gene1'], df['gene2']))
    
    # Extract LCC
    lcc = max(nx.connected_components(G), key=len)
    return G.subgraph(lcc).copy()


def filter_to_tissue(G, tissue_genes):
    """
    Filter network to tissue-expressed genes.
    
    Args:
        G: NetworkX graph
        tissue_genes: Set of gene symbols
        
    Returns:
        Filtered graph (LCC)
    """
    nodes = [n for n in G.nodes() if n in tissue_genes]
    G_tissue = G.subgraph(nodes).copy()
    
    if len(G_tissue) > 0:
        components = list(nx.connected_components(G_tissue))
        if len(components) > 1:
            lcc = max(components, key=len)
            return G_tissue.subgraph(lcc).copy()
    
    return G_tissue

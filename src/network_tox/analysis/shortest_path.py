"""Shortest path analysis."""

import numpy as np
import networkx as nx

def calculate_shortest_path(G, drug_targets, disease_genes):
    """
    Calculate shortest-path proximity (d_c).

    Formula: dc = 1/|T| * sum(min(dist(t, d)) for t in T) where d in D

    Args:
        G: NetworkX graph
        drug_targets: List of target genes
        disease_genes: List of disease genes

    Returns:
        Mean minimum distance
    """
    targets_in = [t for t in drug_targets if t in G]
    disease_in = [d for d in disease_genes if d in G]

    if not targets_in or not disease_in:
        return np.nan

    distances = []
    for target in targets_in:
        min_dist = float('inf')
        for disease in disease_in:
            try:
                dist = nx.shortest_path_length(G, target, disease)
                min_dist = min(min_dist, dist)
            except nx.NetworkXNoPath:
                continue
        if min_dist != float('inf'):
            distances.append(min_dist)

    return np.mean(distances) if distances else np.nan

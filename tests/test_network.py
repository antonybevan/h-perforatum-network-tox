"""
Unit tests for network filtering functions.
"""

import pytest
import networkx as nx
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from network_tox.core.network import filter_to_tissue


def test_filter_to_tissue_basic():
    """Test basic tissue filtering."""
    # Create simple test network
    G = nx.Graph()
    G.add_edges_from([('A', 'B'), ('B', 'C'), ('C', 'D')])
    
    # Filter to subset
    tissue_genes = {'A', 'B', 'C'}
    
    G_filtered = filter_to_tissue(G, tissue_genes)
    
    assert G_filtered.number_of_nodes() == 3
    assert 'D' not in G_filtered.nodes()
    assert ('A', 'B') in G_filtered.edges()


def test_filter_to_tissue_empty():
    """Test filtering with empty tissue set."""
    G = nx.Graph()
    G.add_edge('A', 'B')
    
    G_filtered = filter_to_tissue(G, set())
    
    assert G_filtered.number_of_nodes() == 0


def test_filter_to_tissue_all():
    """Test filtering where all nodes are in tissue."""
    G = nx.Graph()
    G.add_edges_from([('A', 'B'), ('B', 'C')])
    
    tissue_genes = {'A', 'B', 'C'}
    
    G_filtered = filter_to_tissue(G, tissue_genes)
    
    assert G_filtered.number_of_nodes() == 3
    assert G_filtered.number_of_edges() == 2


def test_network_lcc_extraction():
    """Test largest connected component extraction."""
    # Create network with multiple components
    G = nx.Graph()
    G.add_edges_from([('A', 'B'), ('B', 'C')])  # Component 1
    G.add_edges_from([('X', 'Y')])  # Component 2
    
    # Get LCC
    lcc = max(nx.connected_components(G), key=len)
    G_lcc = G.subgraph(lcc).copy()
    
    assert G_lcc.number_of_nodes() == 3
    assert 'X' not in G_lcc.nodes()

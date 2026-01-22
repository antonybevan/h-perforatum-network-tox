"""
Unit tests for proximity calculations.
"""

import numpy as np
import networkx as nx
from network_tox.core import proximity


def test_calculate_shortest_path():
    """Test shortest-path proximity calculation."""
    # Create simple test graph
    G = nx.Graph()
    edges = [('A', 'B'), ('B', 'C'), ('C', 'D'), ('D', 'E')]
    G.add_edges_from(edges)
    
    targets = ['A', 'B']
    disease = ['D', 'E']
    
    d_c = proximity.calculate_shortest_path(G, targets, disease)
    
    # A->D = 3, A->E = 4, min = 3
    # B->D = 2, B->E = 3, min = 2
    # Mean = (3 + 2) / 2 = 2.5
    assert d_c == 2.5


def test_calculate_shortest_path_no_overlap():
    """Test when no targets in network."""
    G = nx.Graph()
    G.add_edges_from([('A', 'B')])
    
    d_c = proximity.calculate_shortest_path(G, ['X'], ['A'])
    assert np.isnan(d_c)


def test_calculate_shortest_path_disconnected():
    """Test with disconnected components."""
    G = nx.Graph()
    G.add_edges_from([('A', 'B'), ('C', 'D')])
    
    d_c = proximity.calculate_shortest_path(G, ['A'], ['C'])
    assert np.isnan(d_c)

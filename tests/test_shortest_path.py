"""Tests for shortest path analysis."""

import numpy as np
import networkx as nx
from network_tox.analysis import shortest_path

def test_calculate_shortest_path_simple():
    """Test simple shortest path calculation."""
    # A - B - C - D
    G = nx.Graph()
    G.add_edges_from([('A', 'B'), ('B', 'C'), ('C', 'D')])

    # Target: A, Disease: D
    # Dist: 3
    d_c = shortest_path.calculate_shortest_path(G, ['A'], ['D'])
    assert d_c == 3.0

    # Target: A, Disease: C, D
    # A->C = 2, A->D = 3. Min = 2.
    d_c = shortest_path.calculate_shortest_path(G, ['A'], ['C', 'D'])
    assert d_c == 2.0

    # Target: A, B. Disease: D
    # A->D = 3. B->D = 2.
    # Mean(3, 2) = 2.5
    d_c = shortest_path.calculate_shortest_path(G, ['A', 'B'], ['D'])
    assert d_c == 2.5

def test_calculate_shortest_path_disconnected():
    """Test disconnected graph."""
    G = nx.Graph()
    G.add_edges_from([('A', 'B'), ('C', 'D')])

    # A->C impossible
    d_c = shortest_path.calculate_shortest_path(G, ['A'], ['C'])
    assert np.isnan(d_c)

def test_calculate_shortest_path_no_overlap():
    """Test no overlap with network."""
    G = nx.Graph()
    G.add_edges_from([('A', 'B')])

    d_c = shortest_path.calculate_shortest_path(G, ['X'], ['A'])
    assert np.isnan(d_c)

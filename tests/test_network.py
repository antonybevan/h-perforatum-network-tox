"""
Unit tests for core network operations.

Run with: pytest tests/
"""

import pytest
import networkx as nx
from network_tox.core import network


def test_load_string_network():
    """Test STRING network loading creates single component."""
    # This would require mock data in practice
    pass  # Placeholder for integration test


def test_filter_to_tissue():
    """Test tissue filtering maintains connectivity."""
    # Create test graph
    G = nx.Graph()
    G.add_edges_from([('A', 'B'), ('B', 'C'), ('C', 'D'), ('E', 'F')])
    
    # Filter to subset
    tissue_genes = {'A', 'B', 'C', 'D'}
    G_filtered = network.filter_to_tissue(G, tissue_genes)
    
    # Should return LCC
    assert G_filtered.number_of_nodes() == 4
    assert G_filtered.number_of_edges() == 3
    assert nx.is_connected(G_filtered)


def test_filter_to_tissue_preserves_lcc():
    """Test that only LCC is returned when multiple components."""
    G = nx.Graph()
    G.add_edges_from([('A', 'B'), ('B', 'C'), ('D', 'E')])
    
    tissue_genes = {'A', 'B', 'C', 'D', 'E'}
    G_filtered = network.filter_to_tissue(G, tissue_genes)
    
    # Should return larger component only
    assert G_filtered.number_of_nodes() == 3
    assert nx.is_connected(G_filtered)

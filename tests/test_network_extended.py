"""
Extended tests for Network module.
Targets: src/network_tox/core/network.py
Goal: Increase coverage from 19% to 80%+
"""

import pytest
import networkx as nx
import sys
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.network_tox.core.network import filter_to_tissue


class TestFilterToTissue:
    """Test suite for filter_to_tissue function."""
    
    def test_filter_basic(self):
        """Test basic filtering to tissue genes."""
        G = nx.Graph()
        G.add_edges_from([('A', 'B'), ('B', 'C'), ('C', 'D'), ('D', 'E')])
        tissue_genes = {'A', 'B', 'C'}
        
        result = filter_to_tissue(G, tissue_genes)
        
        assert set(result.nodes()) == {'A', 'B', 'C'}
        assert result.number_of_edges() == 2  # A-B, B-C
    
    def test_filter_empty_graph(self):
        """Test filtering empty graph."""
        G = nx.Graph()
        tissue_genes = {'A', 'B'}
        
        result = filter_to_tissue(G, tissue_genes)
        
        assert len(result) == 0
    
    def test_filter_no_matching_genes(self):
        """Test when no genes match tissue."""
        G = nx.Graph()
        G.add_edges_from([('A', 'B'), ('B', 'C')])
        tissue_genes = {'X', 'Y', 'Z'}
        
        result = filter_to_tissue(G, tissue_genes)
        
        assert len(result) == 0
    
    def test_filter_all_genes_match(self):
        """Test when all genes match tissue."""
        G = nx.Graph()
        G.add_edges_from([('A', 'B'), ('B', 'C'), ('C', 'A')])
        tissue_genes = {'A', 'B', 'C', 'D', 'E'}  # Superset
        
        result = filter_to_tissue(G, tissue_genes)
        
        assert set(result.nodes()) == {'A', 'B', 'C'}
    
    def test_filter_extracts_lcc(self):
        """Test that filtering extracts largest connected component."""
        G = nx.Graph()
        # Large component: A-B-C-D-E
        G.add_edges_from([('A', 'B'), ('B', 'C'), ('C', 'D'), ('D', 'E')])
        # Small component: X-Y
        G.add_edges_from([('X', 'Y')])
        
        tissue_genes = {'A', 'B', 'C', 'X', 'Y'}
        
        result = filter_to_tissue(G, tissue_genes)
        
        # Should keep only LCC (A-B-C)
        assert 'A' in result
        assert 'B' in result
        assert 'C' in result
        # X-Y is smaller component, should be excluded
        assert 'X' not in result
        assert 'Y' not in result
    
    def test_filter_equal_components(self):
        """Test with equally sized components (keeps one)."""
        G = nx.Graph()
        G.add_edges_from([('A', 'B')])
        G.add_edges_from([('X', 'Y')])
        
        tissue_genes = {'A', 'B', 'X', 'Y'}
        
        result = filter_to_tissue(G, tissue_genes)
        
        # Should have 2 nodes (one component)
        assert len(result) == 2
    
    def test_filter_single_node(self):
        """Test filtering results in single isolated node."""
        G = nx.Graph()
        G.add_edges_from([('A', 'B'), ('B', 'C')])
        tissue_genes = {'A'}  # Only matches one node
        
        result = filter_to_tissue(G, tissue_genes)
        
        # Single isolated node has no edges
        assert len(result) == 1
        assert 'A' in result
    
    def test_filter_preserves_graph_structure(self):
        """Test that filtering preserves edge structure."""
        G = nx.Graph()
        G.add_edges_from([('A', 'B'), ('B', 'C'), ('A', 'C')])  # Triangle
        tissue_genes = {'A', 'B', 'C'}
        
        result = filter_to_tissue(G, tissue_genes)
        
        # All edges should be preserved
        assert result.has_edge('A', 'B')
        assert result.has_edge('B', 'C')
        assert result.has_edge('A', 'C')
    
    def test_filter_breaks_edges(self):
        """Test that filtering removes edges when nodes are filtered."""
        G = nx.Graph()
        G.add_edges_from([('A', 'B'), ('B', 'C'), ('C', 'D')])
        tissue_genes = {'A', 'B', 'D'}  # Missing C
        
        result = filter_to_tissue(G, tissue_genes)
        
        # A-B connected, D isolated
        # LCC is A-B
        assert 'A' in result
        assert 'B' in result
        assert 'D' not in result  # Isolated after C removed
    
    def test_filter_with_set_input(self):
        """Test that function works with set input."""
        G = nx.Graph()
        G.add_edges_from([('A', 'B'), ('B', 'C')])
        tissue_genes = set(['A', 'B', 'C'])
        
        result = filter_to_tissue(G, tissue_genes)
        
        assert len(result) == 3
    
    def test_filter_with_list_input(self):
        """Test that function works with list input."""
        G = nx.Graph()
        G.add_edges_from([('A', 'B'), ('B', 'C')])
        tissue_genes = ['A', 'B', 'C']
        
        result = filter_to_tissue(G, tissue_genes)
        
        assert len(result) == 3


class TestFilterToTissueEdgeCases:
    """Edge cases for filter_to_tissue."""
    
    def test_filter_large_graph(self):
        """Test filtering on larger graph."""
        G = nx.barabasi_albert_graph(100, 3, seed=42)
        # Keep about half the nodes
        tissue_genes = set(range(0, 50))
        
        result = filter_to_tissue(G, tissue_genes)
        
        # Should have some nodes
        assert len(result) > 0
        # All nodes should be in tissue genes
        assert all(n in tissue_genes for n in result.nodes())
    
    def test_filter_returns_copy(self):
        """Test that filtering returns a copy, not a view."""
        G = nx.Graph()
        G.add_edges_from([('A', 'B'), ('B', 'C')])
        tissue_genes = {'A', 'B', 'C'}
        
        result = filter_to_tissue(G, tissue_genes)
        
        # Modifying result shouldn't affect original
        result.add_node('NEW')
        assert 'NEW' not in G

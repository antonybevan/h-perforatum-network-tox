"""
Extended tests for RWR (Random Walk with Restart) module.
Targets: src/network_tox/analysis/rwr.py
Goal: Increase coverage from 14% to 80%+
"""

import pytest
import networkx as nx
import numpy as np
import sys
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.network_tox.analysis.rwr import run_rwr


class TestRunRWR:
    """Test suite for run_rwr function."""
    
    def test_rwr_empty_graph(self):
        """Test RWR on empty graph returns empty dict."""
        G = nx.Graph()
        result = run_rwr(G, ['A'])
        assert result == {}
    
    def test_rwr_single_node_graph(self):
        """Test RWR on single node graph."""
        G = nx.Graph()
        G.add_node('A')
        result = run_rwr(G, ['A'])
        assert 'A' in result
        assert result['A'] > 0.5  # Single isolated node with restart_prob=0.7
    
    def test_rwr_linear_graph(self):
        """Test RWR on linear graph (A-B-C-D)."""
        G = nx.path_graph(4)  # 0-1-2-3
        result = run_rwr(G, [0])
        
        # Seed should have highest score
        assert result[0] > result[3]
        # Scores should decrease with distance
        assert result[0] > result[1] > result[2] > result[3]
    
    def test_rwr_star_graph(self):
        """Test RWR on star graph (hub connected to leaves)."""
        G = nx.star_graph(5)  # 0 is hub, 1-5 are leaves
        result = run_rwr(G, [0])
        
        # Hub should have high score
        assert result[0] > 0.5
        # All leaves should have similar scores
        leaf_scores = [result[i] for i in range(1, 6)]
        assert max(leaf_scores) - min(leaf_scores) < 0.1
    
    def test_rwr_multiple_seeds(self):
        """Test RWR with multiple seed nodes."""
        G = nx.cycle_graph(6)  # 0-1-2-3-4-5-0
        result = run_rwr(G, [0, 3])  # Opposite sides
        
        # Both seeds should have high scores
        assert result[0] > 0.1
        assert result[3] > 0.1
        # Middle nodes should have lower scores
        assert result[0] >= result[1]
        assert result[3] >= result[2]
    
    def test_rwr_no_valid_seeds(self):
        """Test RWR when seeds don't exist in graph."""
        G = nx.path_graph(3)  # 0-1-2
        result = run_rwr(G, ['X', 'Y', 'Z'])  # Non-existent nodes
        
        # Should return zero scores for all nodes
        assert all(v == 0.0 for v in result.values())
    
    def test_rwr_partial_valid_seeds(self):
        """Test RWR when only some seeds exist in graph."""
        G = nx.path_graph(3)  # 0-1-2
        result = run_rwr(G, [0, 'nonexistent'])
        
        # Should work with valid seed
        assert result[0] > 0
        assert 0 in result
    
    def test_rwr_scores_sum_approximately_one(self):
        """Test that RWR scores sum to approximately 1."""
        G = nx.barabasi_albert_graph(50, 3, seed=42)
        result = run_rwr(G, [0, 1, 2])
        
        total = sum(result.values())
        assert 0.99 < total < 1.01  # Should be close to 1
    
    def test_rwr_restart_prob_effect(self):
        """Test that higher restart_prob keeps walker closer to seeds."""
        G = nx.path_graph(10)  # Long path
        
        high_restart = run_rwr(G, [0], restart_prob=0.9)
        low_restart = run_rwr(G, [0], restart_prob=0.3)
        
        # High restart should concentrate more on seed
        assert high_restart[0] > low_restart[0]
        # Low restart should spread more to far nodes
        assert high_restart[9] < low_restart[9]
    
    def test_rwr_convergence(self):
        """Test that RWR converges within iteration limit."""
        G = nx.barabasi_albert_graph(100, 3, seed=42)
        result = run_rwr(G, [0], max_iter=100)
        
        # Should produce valid results
        assert len(result) == 100
        assert all(0 <= v <= 1 for v in result.values())
    
    def test_rwr_disconnected_graph(self):
        """Test RWR on graph with disconnected components."""
        G = nx.Graph()
        G.add_edges_from([('A', 'B'), ('B', 'C')])  # Component 1
        G.add_edges_from([('X', 'Y'), ('Y', 'Z')])  # Component 2
        
        result = run_rwr(G, ['A'])
        
        # Seed component should have positive scores
        assert result['A'] > 0
        assert result['B'] > 0
        # Disconnected component should have zero (or very low from restart)
        # Actually in RWR with restart, disconnected nodes get 0
        assert result['X'] == 0 or result['X'] < 0.01
    
    def test_rwr_weighted_vs_unweighted(self):
        """Test that function works on unweighted graphs."""
        G = nx.Graph()
        G.add_weighted_edges_from([('A', 'B', 1), ('B', 'C', 10)])
        
        # Current implementation ignores weights
        result = run_rwr(G, ['A'])
        assert 'A' in result
        assert 'B' in result
        assert 'C' in result


class TestRWREdgeCases:
    """Edge case tests for RWR."""
    
    def test_rwr_self_loop(self):
        """Test RWR with self-loops in graph."""
        G = nx.Graph()
        G.add_edges_from([('A', 'B'), ('B', 'C')])
        G.add_edge('B', 'B')  # Self-loop
        
        result = run_rwr(G, ['A'])
        assert 'A' in result
    
    def test_rwr_complete_graph(self):
        """Test RWR on complete graph."""
        G = nx.complete_graph(5)
        result = run_rwr(G, [0])
        
        # All nodes connected equally
        assert len(result) == 5
        # Seed should still have highest score
        assert result[0] >= max(result[1], result[2], result[3], result[4])

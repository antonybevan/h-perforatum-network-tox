"""
Unit tests for Random Walk with Restart - CI compatible.
"""

import pytest
import networkx as nx
import numpy as np
import sys
from pathlib import Path

# Add src to path for imports
sys.path.append(str(Path(__file__).resolve().parent.parent / 'src'))

from network_tox.analysis.rwr import run_rwr


def test_rwr_returns_scores():
    """Test that RWR returns proper scores."""
    G = nx.Graph()
    G.add_edges_from([('A', 'B'), ('B', 'C'), ('C', 'D')])
    
    seeds = ['A']
    scores = run_rwr(G, seeds, restart_prob=0.15)
    
    assert isinstance(scores, dict)
    assert len(scores) > 0
    assert all(0 <= v <= 1 for v in scores.values())


def test_rwr_seed_has_high_score():
    """Test that seed node has high score."""
    G = nx.Graph()
    G.add_edges_from([('A', 'B'), ('B', 'C'), ('C', 'D')])
    
    seeds = ['B']
    scores = run_rwr(G, seeds, restart_prob=0.15)
    
    # Seed should have highest or near-highest score
    assert scores['B'] >= max(scores['A'], scores['D'])


def test_rwr_multiple_seeds():
    """Test RWR with multiple seeds."""
    G = nx.Graph()
    G.add_edges_from([('A', 'B'), ('B', 'C'), ('C', 'D')])
    
    seeds = ['A', 'D']
    scores = run_rwr(G, seeds, restart_prob=0.15)
    
    assert isinstance(scores, dict)
    assert len(scores) == 4


def test_rwr_scores_sum_to_one():
    """Test that RWR scores sum to approximately 1."""
    G = nx.Graph()
    G.add_edges_from([('A', 'B'), ('B', 'C')])
    
    seeds = ['A']
    scores = run_rwr(G, seeds, restart_prob=0.15)
    
    total = sum(scores.values())
    assert abs(total - 1.0) < 1e-5

def test_rwr_convergence():
    """Test that RWR converges."""
    G = nx.path_graph(10)
    seeds = [0]

    # Run with standard tolerance
    scores = run_rwr(G, seeds, restart_prob=0.15, tol=1e-6)

    # Verify values are stable (simple check: sum is 1)
    assert abs(sum(scores.values()) - 1.0) < 1e-6

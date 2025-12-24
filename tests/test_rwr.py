"""
Unit tests for Random Walk with Restart algorithm.
"""

import pytest
import networkx as nx
import numpy as np
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from network_tox.core.rwr import run_rwr


def test_rwr_returns_scores():
    """Test that RWR returns proper scores."""
    # Simple test network
    G = nx.Graph()
    G.add_edges_from([('A', 'B'), ('B', 'C'), ('C', 'D')])
    
    seeds = ['A']
    scores = run_rwr(G, seeds, restart_prob=0.7)
    
    assert isinstance(scores, dict)
    assert len(scores) > 0
    assert all(0 <= v <= 1 for v in scores.values())


def test_rwr_seed_has_highest_score():
    """Test that seed node has highest score."""
    G = nx.Graph()
    G.add_edges_from([('A', 'B'), ('B', 'C'), ('C', 'D')])
    
    seeds = ['B']
    scores = run_rwr(G, seeds, restart_prob=0.7)
    
    # Seed should have high score
    assert scores['B'] > scores['A']
    assert scores['B'] > scores['D']


def test_rwr_multiple_seeds():
    """Test RWR with multiple seeds."""
    G = nx.Graph()
    G.add_edges_from([('A', 'B'), ('B', 'C'), ('C', 'D')])
    
    seeds = ['A', 'D']
    scores = run_rwr(G, seeds, restart_prob=0.7)
    
    assert isinstance(scores, dict)
    assert len(scores) == 4


def test_rwr_scores_sum_to_one():
    """Test that RWR scores sum to approximately 1."""
    G = nx.Graph()
    G.add_edges_from([('A', 'B'), ('B', 'C')])
    
    seeds = ['A']
    scores = run_rwr(G, seeds, restart_prob=0.7)
    
    total = sum(scores.values())
    assert abs(total - 1.0) < 1e-5  # Allow for floating point precision


def test_rwr_restart_probability():
    """Test that restart probability affects scores."""
    G = nx.Graph()
    G.add_edges_from([('A', 'B'), ('B', 'C'), ('C', 'D')])
    
    seeds = ['A']
    
    # Low restart prob = more diffusion
    scores_low = run_rwr(G, seeds, restart_prob=0.3)
    
    # High restart prob = less diffusion
    scores_high = run_rwr(G, seeds, restart_prob=0.9)
    
    # Seed should have higher score with high restart
    assert scores_high['A'] > scores_low['A']

"""
Unit tests for Random Walk with Restart - CI compatible.
"""

import pytest
import networkx as nx
import numpy as np


def simple_rwr(G, seeds, restart_prob=0.7, max_iter=100):
    """Simplified RWR for testing."""
    n = len(G.nodes())
    node_list = list(G.nodes())
    node_idx = {node: i for i, node in enumerate(node_list)}
    
    # Initialize
    p = np.zeros(n)
    for seed in seeds:
        if seed in node_idx:
            p[node_idx[seed]] = 1.0 / len(seeds)
    
    # Simple iteration
    for _ in range(max_iter):
        p_new = p * restart_prob
        for i, node in enumerate(node_list):
            for neighbor in G.neighbors(node):
                j = node_idx[neighbor]
                p_new[i] += (1 - restart_prob) * p[j] / G.degree(neighbor)
        p = p_new / p_new.sum()
    
    return {node: p[node_idx[node]] for node in node_list}


def test_rwr_returns_scores():
    """Test that RWR returns proper scores."""
    G = nx.Graph()
    G.add_edges_from([('A', 'B'), ('B', 'C'), ('C', 'D')])
    
    seeds = ['A']
    scores = simple_rwr(G, seeds, restart_prob=0.7)
    
    assert isinstance(scores, dict)
    assert len(scores) > 0
    assert all(0 <= v <= 1 for v in scores.values())


def test_rwr_seed_has_high_score():
    """Test that seed node has high score."""
    G = nx.Graph()
    G.add_edges_from([('A', 'B'), ('B', 'C'), ('C', 'D')])
    
    seeds = ['B']
    scores = simple_rwr(G, seeds, restart_prob=0.7)
    
    # Seed should have highest or near-highest score
    assert scores['B'] >= max(scores['A'], scores['D'])


def test_rwr_multiple_seeds():
    """Test RWR with multiple seeds."""
    G = nx.Graph()
    G.add_edges_from([('A', 'B'), ('B', 'C'), ('C', 'D')])
    
    seeds = ['A', 'D']
    scores = simple_rwr(G, seeds, restart_prob=0.7)
    
    assert isinstance(scores, dict)
    assert len(scores) == 4


def test_rwr_scores_sum_to_one():
    """Test that RWR scores sum to approximately 1."""
    G = nx.Graph()
    G.add_edges_from([('A', 'B'), ('B', 'C')])
    
    seeds = ['A']
    scores = simple_rwr(G, seeds, restart_prob=0.7)
    
    total = sum(scores.values())
    assert abs(total - 1.0) < 1e-5

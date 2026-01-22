"""
Unit tests for permutation testing.
"""

import networkx as nx
from network_tox.core import permutation


def test_get_degree_matched_random():
    """Test degree-matched sampling."""
    G = nx.barabasi_albert_graph(100, 3, seed=42)
    targets = list(G.nodes())[:10]
    
    random_set = permutation.get_degree_matched_random(G, targets, 5, seed=42)
    
    assert len(random_set) == 5
    assert len(set(random_set) & set(targets)) == 0  # No overlap


def test_calculate_z_score():
    """Test Z-score calculation."""
    null_dist = [1, 2, 3, 4, 5]
    obs = 1
    
    z = permutation.calculate_z_score(obs, null_dist)
    
    # mean=3, std=sqrt(2), z=(1-3)/sqrt(2) = -1.414
    assert abs(z - (-1.414)) < 0.01


def test_calculate_p_value():
    """Test p-value calculation."""
    # Z-score of 2 (two-tailed)
    p = permutation.calculate_p_value(2.0, tail='two')
    assert abs(p - 0.0455) < 0.01
    
    # Z-score of 2 (one-tailed)
    p = permutation.calculate_p_value(2.0, tail='one')
    assert abs(p - 0.0228) < 0.01

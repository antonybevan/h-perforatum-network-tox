"""Tests for RWR analysis."""

import networkx as nx
from network_tox.analysis import rwr

def test_run_rwr_simple():
    """Test RWR on a simple graph."""
    # A - B - C
    G = nx.Graph()
    G.add_edges_from([('A', 'B'), ('B', 'C')])

    # Restart at A
    scores = rwr.run_rwr(G, ['A'], restart_prob=0.5)

    # A should be highest
    assert scores['A'] > scores['B']
    assert scores['B'] > scores['C']

    # Sum should be approx 1
    total = sum(scores.values())
    assert abs(total - 1.0) < 1e-6

def test_run_rwr_convergence():
    """Test convergence."""
    G = nx.cycle_graph(10)
    nodes = list(G.nodes())
    scores = rwr.run_rwr(G, [nodes[0]], restart_prob=0.1)

    assert len(scores) == 10
    assert abs(sum(scores.values()) - 1.0) < 1e-6

def test_run_rwr_no_seeds():
    """Test with no valid seeds."""
    G = nx.path_graph(5)
    scores = rwr.run_rwr(G, ['INVALID'])

    assert all(s == 0.0 for s in scores.values())

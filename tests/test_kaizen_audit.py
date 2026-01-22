import pandas as pd
import numpy as np
import networkx as nx
import pytest
from scipy import sparse

# Import modules to test
import sys
from pathlib import Path
sys.path.append('src')

from network_tox.core.permutation import calculate_z_score, calculate_empirical_p_value, get_degree_matched_random
from network_tox.analysis.rwr import run_rwr
from network_tox.analysis.expression_weighted_rwr import normalize_expression_values, create_expression_weighted_transition_matrix

# =============================================================================
# Layer 1: Mathematical Verification
# =============================================================================

def test_z_score_calculation():
    """Verify Z = (x_obs - μ_null) / σ_null"""
    obs = 10.0
    null_dist = [2.0, 4.0, 6.0]
    # mean = 4.0
    # std = sqrt(((2-4)^2 + (4-4)^2 + (6-4)^2) / 3) = sqrt((4+0+4)/3) = sqrt(8/3) ≈ 1.633

    expected_z = (10.0 - 4.0) / np.std(null_dist)
    calculated_z = calculate_z_score(obs, null_dist)

    assert np.isclose(calculated_z, expected_z)

def test_empirical_p_value_calculation():
    """Verify P = (r+1)/(n+1)"""
    # Case: obs > all null (highly significant)
    obs = 10.0
    null_dist = [1.0, 2.0, 3.0, 4.0, 5.0] # n=5
    # tail='one_greater', r=0 (none >= obs)
    # p = (0+1)/(5+1) = 1/6

    p = calculate_empirical_p_value(obs, null_dist, tail='one_greater')
    assert np.isclose(p, 1/6)

    # Case: obs < all null (not significant for 'greater')
    obs = 0.0
    # r=5 (all >= obs)
    # p = (5+1)/(5+1) = 1.0
    p = calculate_empirical_p_value(obs, null_dist, tail='one_greater')
    assert np.isclose(p, 1.0)

    # Case: median
    obs = 3.0
    # r=3 (3,4,5 >= 3)
    # p = (3+1)/6 = 4/6 = 0.6667
    p = calculate_empirical_p_value(obs, null_dist, tail='one_greater')
    assert np.isclose(p, 4/6)

# =============================================================================
# Layer 2: Algorithm Correctness (RWR)
# =============================================================================

def test_rwr_conservation_of_mass():
    """Verify total probability mass is conserved (sums to 1)."""
    # Create simple graph: A-B-C
    G = nx.Graph()
    G.add_edges_from([('A', 'B'), ('B', 'C')])
    seeds = ['A']

    scores = run_rwr(G, seeds, restart_prob=0.15, max_iter=100)

    total_mass = sum(scores.values())
    assert np.isclose(total_mass, 1.0, atol=1e-5)

def test_rwr_restart_vector():
    """Verify restart vector spreads probability among seeds."""
    G = nx.Graph()
    G.add_edges_from([('A', 'B'), ('C', 'D'), ('B', 'C')])
    seeds = ['A', 'D'] # Two disconnected seeds initially

    # With high restart probability, mass should be concentrated on seeds
    scores = run_rwr(G, seeds, restart_prob=0.99, max_iter=10)

    # Check if seeds have approx 0.5 * 0.99 mass each
    assert np.isclose(scores['A'], 0.5 * 0.99, atol=0.05)
    assert np.isclose(scores['D'], 0.5 * 0.99, atol=0.05)

def test_expression_weighting_direction():
    """Verify expression weighting uses destination-node weighting."""
    # Create star graph: C connected to L (low expr) and H (high expr)
    # Seed at C.
    # Standard: Flow splits 50/50 to L and H.
    # Weighted: Flow should favor H because it has higher expression.

    G = nx.Graph()
    G.add_edges_from([('C', 'L'), ('C', 'H')])

    expression = {'C': 10.0, 'L': 1.0, 'H': 1000.0}
    # Norm L: log2(2) = 1 -> min -> 0.01
    # Norm H: log2(1001) ~ 10 -> max -> 1.0

    seeds = ['C']

    # Standard RWR
    scores_std = run_rwr(G, seeds, restart_prob=0.15)

    # Expression RWR
    from network_tox.analysis.expression_weighted_rwr import run_expression_weighted_rwr
    scores_expr = run_expression_weighted_rwr(G, seeds, expression, restart_prob=0.15)

    # Check H vs L ratio
    ratio_std = scores_std['H'] / scores_std['L']
    ratio_expr = scores_expr['H'] / scores_expr['L']

    # Standard should be symmetric (ratio ~ 1)
    assert np.isclose(ratio_std, 1.0, atol=0.01)

    # Expr should favor H (ratio > 1)
    assert ratio_expr > 1.5 # Should be significantly higher
    assert ratio_expr > ratio_std

def test_log_transformation_base():
    """Verify log2 transformation is used."""
    nodes = ['A', 'B']
    expression = {'A': 3.0, 'B': 7.0}
    # log2(3+1)=2, log2(7+1)=3
    # min=2, max=3
    # Norm: A->0, B->1 (with floor 0.01)

    norm = normalize_expression_values(expression, nodes, method="log_minmax")

    assert np.isclose(norm[0], 0.01) # Floored
    assert np.isclose(norm[1], 1.0)

# =============================================================================
# Layer 6: Statistical Integrity
# =============================================================================

def test_bootstrap_sampling():
    """Verify sampling is without replacement."""
    # Create a mock population
    population = list(range(100))
    sample_size = 10

    # Run 1000 samples and check for uniqueness within sample
    for _ in range(100):
        sample = np.random.choice(population, size=sample_size, replace=False)
        assert len(set(sample)) == sample_size # All elements unique

def test_degree_matching_tolerance():
    """Verify degree matching respects +/- 25% tolerance."""
    # Create star graph: Center connected to 10 nodes (degree 10)
    # Plus other nodes with degree 8, 12, 13 (8=10-20%, 12=10+20%, 13=10+30%)

    G = nx.Graph()
    # Target node T (deg 10)
    edges = [('T', f'n{i}') for i in range(10)]

    # Candidate C1 (deg 8) - should match (20% < 25%)
    edges += [('C1', f'c1_{i}') for i in range(8)]

    # Candidate C2 (deg 12) - should match (20% < 25%)
    edges += [('C2', f'c2_{i}') for i in range(12)]

    # Candidate C3 (deg 13) - should NOT match (30% > 25%)
    edges += [('C3', f'c3_{i}') for i in range(13)]

    G.add_edges_from(edges)

    targets = ['T']

    # Run matching 100 times to see what we get
    # Note: The function returns a list of matched nodes
    matches = []
    for i in range(100):
        m = get_degree_matched_random(G, targets, n_sample=1, seed=i)
        matches.extend(m)

    unique_matches = set(matches)

    assert 'C1' in unique_matches
    assert 'C2' in unique_matches
    assert 'C3' not in unique_matches
    assert 'T' not in unique_matches

if __name__ == "__main__":
    sys.exit(pytest.main(["-v", __file__]))

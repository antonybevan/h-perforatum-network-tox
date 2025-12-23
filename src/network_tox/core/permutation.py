"""Permutation testing."""

import numpy as np
from scipy import stats


def get_degree_matched_random(G, targets, n_sample, seed=None):
    """
    Get degree-matched random nodes.
    
    Args:
        G: NetworkX graph
        targets: Target nodes
        n_sample: Number to sample
        seed: Random seed
        
    Returns:
        List of random nodes
    """
    if seed is not None:
        np.random.seed(seed)
    
    degrees = dict(G.degree())
    target_degrees = [degrees[t] for t in targets if t in degrees]
    
    random_set = []
    nodes = list(G.nodes())
    
    for deg in target_degrees[:n_sample]:
        tol = max(1, int(deg * 0.15))
        candidates = [n for n in nodes 
                     if abs(degrees[n] - deg) <= tol 
                     and n not in targets 
                     and n not in random_set]
        if not candidates:
            candidates = [n for n in nodes if n not in targets and n not in random_set]
        if candidates:
            random_set.append(np.random.choice(candidates))
    
    return random_set


def calculate_z_score(obs, null_dist):
    """Calculate Z-score from null distribution."""
    null_mean = np.mean(null_dist)
    null_std = np.std(null_dist)
    
    if null_std > 0:
        return (obs - null_mean) / null_std
    return 0.0


def calculate_p_value(z_score, tail='two'):
    """Calculate p-value from Z-score."""
    if tail == 'two':
        return 2 * (1 - stats.norm.cdf(abs(z_score)))
    else:  # one-tailed
        return 1 - stats.norm.cdf(z_score)

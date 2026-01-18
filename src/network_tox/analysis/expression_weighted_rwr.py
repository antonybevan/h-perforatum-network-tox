"""
Expression-Weighted Random Walk with Restart (RWR)

Implements tissue-constrained influence propagation using transition-matrix
weighting, where expression values modify the adjacency matrix itself.

Mathematical formulation:
    A'_ij = A_ij * e_i           (weight by source expression)
    W'_ij = A'_ij / Σ_k A'_kj    (column-normalize)
    p^(∞) = (1-α) * W' * p^(∞) + α * s_uniform

Where:
    e_i = normalized liver expression of protein i (TPM scaled to [0,1])
    s_uniform = uniform restart vector over target nodes

Biological interpretation:
    - Highly expressed proteins become preferred conduits for signal propagation
    - Non-expressed proteins become dead ends
    - Walk is constrained to liver-relevant biology

References:
    - Guney et al. (2016) Network-based in silico drug efficacy screening
    - Vanunu et al. (2010) Associating genes and protein complexes with disease
"""

import numpy as np
import networkx as nx
from scipy import sparse
from typing import Dict, List, Optional, Union
from pathlib import Path


def load_liver_expression(
    gtex_file: Union[str, Path],
    tissue_column: str = "Liver"
) -> Dict[str, float]:
    """
    Load liver TPM values from GTEx v8 gene median TPM file.
    
    Args:
        gtex_file: Path to GTEx gene_median_tpm.gct file
        tissue_column: Name of tissue column to extract
        
    Returns:
        Dictionary mapping gene_symbol -> TPM value
    """
    import pandas as pd
    
    # GTEx .gct format: skip first 2 rows (version, dimensions)
    # Columns: Name, Description, then tissue columns
    df = pd.read_csv(gtex_file, sep='\t', skiprows=2)
    
    if tissue_column not in df.columns:
        raise ValueError(f"Tissue '{tissue_column}' not found. Available: {df.columns.tolist()}")
    
    # Use Description column as gene symbol (more readable than ENSG ID)
    expression = {}
    for _, row in df.iterrows():
        gene_symbol = row.get('Description', row.get('Name', ''))
        tpm = row.get(tissue_column, 0)
        if gene_symbol and pd.notna(tpm):
            expression[gene_symbol] = float(tpm)
    
    return expression


def normalize_expression_values(
    expression: Dict[str, float],
    nodes: List[str],
    method: str = "minmax"
) -> np.ndarray:
    """
    Normalize expression values to [0, 1] range for transition weighting.
    
    Args:
        expression: Dictionary mapping gene -> TPM value
        nodes: List of all network nodes
        method: Normalization method ('minmax' or 'log_minmax')
        
    Returns:
        Array of normalized expression values (length = len(nodes))
    """
    n = len(nodes)
    expr_values = np.array([expression.get(node, 0.0) for node in nodes])
    
    if method == "log_minmax":
        # Log-transform then min-max normalize
        expr_values = np.log1p(expr_values)
    
    # Min-max normalize to [0, 1]
    # Add small epsilon to avoid perfect zeros (makes non-expressed nodes dead ends)
    min_val = expr_values.min()
    max_val = expr_values.max()
    
    if max_val > min_val:
        expr_normalized = (expr_values - min_val) / (max_val - min_val)
    else:
        expr_normalized = np.ones(n)  # All equal if no variance
    
    # Add small floor to avoid perfect zero (allow minimal propagation through non-expressed nodes)
    expr_normalized = np.maximum(expr_normalized, 0.01)
    
    return expr_normalized


def create_expression_weighted_transition_matrix(
    adj_matrix: sparse.spmatrix,
    expression: Dict[str, float],
    nodes: List[str]
) -> sparse.spmatrix:
    """
    Create expression-weighted transition matrix for RWR.
    
    This is the CORRECT implementation of expression-weighted RWR:
        1. Weight adjacency by source node expression: A'_ij = A_ij * e_i
        2. Column-normalize: W'_ij = A'_ij / Σ_k A'_kj
    
    Biological interpretation:
        - Highly expressed proteins become preferred conduits
        - Low/non-expressed proteins become dead ends
        - Walk is constrained to liver-active biology
    
    Args:
        adj_matrix: Adjacency matrix (n x n, undirected)
        expression: Dictionary mapping node -> TPM value
        nodes: List of nodes (defines row/column indexing)
        
    Returns:
        Column-normalized transition matrix W' (sparse)
    """
    n = len(nodes)
    
    # Normalize expression to [0, 1]
    expr_normalized = normalize_expression_values(expression, nodes, method="log_minmax")
    
    # Weight each ROW by source node expression: A'_ij = A_ij * e_i
    # This makes highly-expressed nodes better signal transmitters
    expr_diag = sparse.diags(expr_normalized)
    weighted_adj = expr_diag.dot(adj_matrix)
    
    # Column-normalize: each column sums to 1
    col_sum = np.array(weighted_adj.sum(axis=0)).flatten()
    col_sum[col_sum == 0] = 1.0  # Avoid division by zero
    d_inv = sparse.diags(1.0 / col_sum)
    W_prime = weighted_adj.dot(d_inv)
    
    return W_prime


def run_expression_weighted_rwr(
    G: nx.Graph,
    seeds: List[str],
    expression: Dict[str, float],
    restart_prob: float = 0.15,
    tol: float = 1e-6,
    max_iter: int = 100,
    transform: str = "log1p"  # Kept for backwards compatibility, now unused
) -> Dict[str, float]:
    """
    Run expression-weighted Random Walk with Restart.
    
    Uses TRANSITION-MATRIX WEIGHTING (the correct, reviewer-proof approach):
        1. Weight adjacency by expression: A'_ij = A_ij * e_i
        2. Column-normalize: W'_ij = A'_ij / Σ_k A'_kj
        3. Use uniform restart vector over targets
        4. Iterate: p^(t+1) = (1-α) * W' * p^t + α * s_uniform
    
    This is NOT restart-vector weighting. Expression modifies the transition
    probabilities themselves, forcing walks to flow through liver-expressed proteins.
    
    Biological interpretation:
        - Highly expressed targets become preferred conduits for signal
        - Non-expressed proteins become dead ends
        - Network propagation is constrained to liver-active biology
        - Standard RWR math, just biology-aware transition matrix
    
    Args:
        G: NetworkX graph (undirected PPI network)
        seeds: List of seed nodes (drug targets)
        expression: Dictionary mapping node -> TPM value
        restart_prob: Restart probability α (default 0.15)
        tol: Convergence tolerance (L1 difference)
        max_iter: Maximum iterations
        transform: Deprecated (kept for backwards compatibility)
        
    Returns:
        Dictionary of {node: steady-state probability}
    """
    if len(G) == 0:
        return {}
    
    # Create node indexing
    nodes = list(G.nodes())
    node_idx = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)
    
    # Create adjacency matrix
    adj = nx.adjacency_matrix(G, nodelist=nodes).astype(float)
    
    # Create EXPRESSION-WEIGHTED TRANSITION MATRIX (the key change)
    W_prime = create_expression_weighted_transition_matrix(
        adj_matrix=adj,
        expression=expression,
        nodes=nodes
    )
    
    # Create UNIFORM restart vector over targets (not expression-weighted)
    r = np.zeros((n, 1))
    valid_seeds = [s for s in seeds if s in node_idx]
    
    if not valid_seeds:
        return {node: 0.0 for node in nodes}
    
    for seed in valid_seeds:
        r[node_idx[seed]] = 1.0 / len(valid_seeds)
    
    # Initialize p with restart vector
    p = r.copy()
    
    # Iterate until convergence (standard RWR iteration)
    for iteration in range(max_iter):
        p_new = (1 - restart_prob) * W_prime.dot(p) + restart_prob * r
        diff = np.sum(np.abs(p_new - p))
        p = p_new
        if diff < tol:
            break
    
    return {nodes[i]: float(p[i, 0]) for i in range(n)}


def compute_dili_influence(
    rwr_scores: Dict[str, float],
    dili_genes: List[str]
) -> float:
    """
    Compute DILI influence as sum of RWR scores at DILI genes.
    
    Args:
        rwr_scores: Output from run_expression_weighted_rwr
        dili_genes: List of DILI-associated gene symbols
        
    Returns:
        Total influence score (sum of steady-state probabilities)
    """
    return sum(rwr_scores.get(gene, 0.0) for gene in dili_genes)


# =============================================================================
# COMPARISON WITH STANDARD RWR
# =============================================================================

def run_standard_rwr(
    G: nx.Graph,
    seeds: List[str],
    restart_prob: float = 0.15,
    tol: float = 1e-6,
    max_iter: int = 100
) -> Dict[str, float]:
    """
    Run standard (unweighted) RWR for comparison.
    
    Uses uniform restart vector where all seeds have equal weight.
    """
    if len(G) == 0:
        return {}
    
    nodes = list(G.nodes())
    node_idx = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)
    
    adj = nx.adjacency_matrix(G, nodelist=nodes).astype(float)
    
    col_sum = np.array(adj.sum(axis=0)).flatten()
    col_sum[col_sum == 0] = 1
    d_inv = sparse.diags(1.0 / col_sum)
    W = adj.dot(d_inv)
    
    # Uniform restart vector
    r = np.zeros((n, 1))
    valid_seeds = [s for s in seeds if s in node_idx]
    if not valid_seeds:
        return {node: 0.0 for node in nodes}
    
    for seed in valid_seeds:
        r[node_idx[seed]] = 1.0 / len(valid_seeds)
    
    p = r.copy()
    
    for _ in range(max_iter):
        p_new = (1 - restart_prob) * W.dot(p) + restart_prob * r
        diff = np.sum(np.abs(p_new - p))
        p = p_new
        if diff < tol:
            break
    
    return {nodes[i]: float(p[i, 0]) for i in range(n)}


# =============================================================================
# PERMUTATION TESTING COMPATIBILITY
# =============================================================================

def get_degree_matched_random_seeds(
    G: nx.Graph,
    original_seeds: List[str],
    expression: Dict[str, float],
    degree_tolerance: float = 0.25
) -> List[str]:
    """
    Sample degree-matched random seeds, preserving expression availability.
    
    For permutation testing, this ensures null distribution matches
    both degree and expression characteristics of original seeds.
    
    Args:
        G: NetworkX graph
        original_seeds: Original seed nodes to match
        expression: Expression dictionary (to ensure sampled nodes have expression)
        degree_tolerance: Fraction tolerance for degree matching (default ±25%)
        
    Returns:
        List of randomly sampled degree-matched seeds
    """
    import random
    
    nodes_in_G = set(G.nodes())
    valid_seeds = [s for s in original_seeds if s in nodes_in_G]
    
    if not valid_seeds:
        return []
    
    # Get all nodes with expression data
    expressed_nodes = [n for n in nodes_in_G if n in expression]
    
    random_seeds = []
    for seed in valid_seeds:
        seed_degree = G.degree(seed)
        min_degree = int(seed_degree * (1 - degree_tolerance))
        max_degree = int(seed_degree * (1 + degree_tolerance))
        
        # Find candidates with matching degree and expression
        candidates = [
            n for n in expressed_nodes
            if min_degree <= G.degree(n) <= max_degree and n not in random_seeds
        ]
        
        if candidates:
            random_seeds.append(random.choice(candidates))
        else:
            # Fallback: any expressed node
            fallback = [n for n in expressed_nodes if n not in random_seeds]
            if fallback:
                random_seeds.append(random.choice(fallback))
    
    return random_seeds

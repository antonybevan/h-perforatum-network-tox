#!/usr/bin/env python3
"""
OPTIMIZED EWI Bootstrap Sensitivity Analysis.
Creates the weighted transition matrix once and runs 100 RWR iterations.
"""
import pandas as pd
import numpy as np
import networkx as nx
from pathlib import Path
import sys
from tqdm import tqdm
from scipy.sparse import issparse

sys.path.append('src')
from network_tox.analysis.expression_weighted_rwr import (
    load_liver_expression,
    create_expression_weighted_transition_matrix
)

DATA_DIR = Path('data')
RESULTS_DIR = Path('results')
GTEX_FILE = DATA_DIR / 'raw' / 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct'

# Configuration
N_BOOTSTRAP = 100
RANDOM_SEED = 42

def run_rwr_on_matrix(W, r, restart_prob=0.15, tol=1e-6, max_iter=100):
    """Internal RWR iteration on pre-built matrix."""
    r = r.reshape(-1, 1)
    p = r.copy()
    for _ in range(max_iter):
        p_new = (1 - restart_prob) * W.dot(p) + restart_prob * r
        diff = np.sum(np.abs(p_new - p))
        p = p_new
        if diff < tol:
            break
    return p.flatten()

def main():
    print("Optimized EWI Bootstrap Sensitivity Analysis")
    np.random.seed(RANDOM_SEED)
    
    # Load expression
    expression = load_liver_expression(GTEX_FILE, tissue_column="Liver")
    
    # Load LCC-filtered targets
    targets_df = pd.read_csv(DATA_DIR / 'processed' / 'targets_lcc.csv')
    hyp_targets = list(targets_df[targets_df['compound'] == 'Hyperforin']['gene_symbol'])
    quer_targets = list(targets_df[targets_df['compound'] == 'Quercetin']['gene_symbol'])
    SAMPLE_SIZE = len(hyp_targets)
    
    # Load network
    net_df = pd.read_parquet(DATA_DIR / 'processed' / 'network_900_liver_lcc.parquet')
    G = nx.from_pandas_edgelist(net_df, 'gene1', 'gene2')
    nodes = list(G.nodes())
    node_idx = {n: i for i, n in enumerate(nodes)}
    adj = nx.adjacency_matrix(G, nodelist=nodes).astype(float)
    
    # Pre-calculate EWI Transition Matrix
    print("Building EWI transition matrix...")
    W_prime = create_expression_weighted_transition_matrix(
        adj_matrix=adj,
        expression=expression,
        nodes=nodes
    )
    
    # Load DILI genes
    dili_df = pd.read_csv(DATA_DIR / 'processed' / 'dili_900_lcc.csv')
    dili_indices = [node_idx[g] for g in dili_df['gene_name'] if g in node_idx]
    
    # Filter targets
    hyp_indices = [node_idx[t] for t in hyp_targets if t in node_idx]
    quer_indices = [node_idx[t] for t in quer_targets if t in node_idx]
    
    # Observed EWI
    p0_hyp = np.zeros(len(nodes))
    p0_hyp[hyp_indices] = 1.0 / len(hyp_indices)
    hyp_scores = run_rwr_on_matrix(W_prime, p0_hyp, restart_prob=0.15)
    hyp_ewi = hyp_scores[dili_indices].sum()
    print(f"Hyperforin Observed EWI: {hyp_ewi:.6f}")
    
    # Bootstrap
    bootstrap_results = []
    print(f"Running {N_BOOTSTRAP} optimized iterations...")
    for i in tqdm(range(N_BOOTSTRAP), desc="EWI Bootstrap"):
        # Sample
        sample_indices = np.random.choice(quer_indices, size=SAMPLE_SIZE, replace=False)
        p0 = np.zeros(len(nodes))
        p0[sample_indices] = 1.0 / SAMPLE_SIZE
        
        # Run RWR on pre-built matrix
        scores = run_rwr_on_matrix(W_prime, p0, restart_prob=0.15)
        influence = scores[dili_indices].sum()
        bootstrap_results.append(influence)
        
    mean_ewi = np.mean(bootstrap_results)
    ratio = hyp_ewi / mean_ewi
    
    print(f"\nResults:")
    print(f"  Quercetin Bootstrap Mean (EWI): {mean_ewi:.6f}")
    print(f"  Hyperforin Robust Ratio (EWI):  {ratio:.2f}x")
    
    # Save results
    summary = pd.DataFrame([{
        'observed_ewi': hyp_ewi,
        'bootstrap_mean_ewi': mean_ewi,
        'robust_ratio': ratio,
        'n_bootstrap': N_BOOTSTRAP
    }])
    summary.to_csv(RESULTS_DIR / 'tables' / 'ewi_bootstrap_summary.csv', index=False)
    print(f"\nSaved to results/tables/ewi_bootstrap_summary.csv")

if __name__ == '__main__':
    main()

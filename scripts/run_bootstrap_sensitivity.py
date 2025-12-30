#!/usr/bin/env python3
"""
Bootstrap Sensitivity Analysis.

Samples N random targets from Quercetin (matching Hyperforin count)
to test if Hyperforin's influence advantage is due to target count.
"""

import pandas as pd
import numpy as np
import networkx as nx
from pathlib import Path
import sys
from tqdm import tqdm

sys.path.append('src')
from network_tox.analysis.rwr import run_rwr

DATA_DIR = Path('data')
RESULTS_DIR = Path('results')

# Configuration
N_BOOTSTRAP = 100
RANDOM_SEED = 42

def main():
    print("=" * 70)
    print(" BOOTSTRAP SENSITIVITY ANALYSIS")
    print("=" * 70)
    print()
    
    np.random.seed(RANDOM_SEED)
    
    # Load LCC-filtered targets
    targets_df = pd.read_csv(DATA_DIR / 'processed' / 'targets_lcc.csv')
    hyp_targets = list(targets_df[targets_df['compound'] == 'Hyperforin']['gene_symbol'])
    quer_targets = list(targets_df[targets_df['compound'] == 'Quercetin']['gene_symbol'])
    
    SAMPLE_SIZE = len(hyp_targets)  # Dynamic - now 10
    
    print(f"Configuration:")
    print(f"  Bootstrap iterations: {N_BOOTSTRAP}")
    print(f"  Sample size: {SAMPLE_SIZE} (matching Hyperforin)")
    print(f"  Random seed: {RANDOM_SEED}")
    print()
    
    # Load network
    net_df = pd.read_parquet(DATA_DIR / 'processed' / 'network_900_liver_lcc.parquet')
    if 'gene1' in net_df.columns:
        G = nx.from_pandas_edgelist(net_df, 'gene1', 'gene2')
    else:
        G = nx.from_pandas_edgelist(net_df, 'protein1', 'protein2')
    
    print(f"Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    
    # Load DILI genes
    dili_df = pd.read_csv(DATA_DIR / 'processed' / 'dili_900_lcc.csv')
    dili_genes = [g for g in dili_df['gene_name'] if g in G]
    print(f"DILI genes in network: {len(dili_genes)}")
    
    # Filter targets to network
    hyp_in_net = [t for t in hyp_targets if t in G]
    quer_in_net = [t for t in quer_targets if t in G]
    
    print(f"Hyperforin targets in network: {len(hyp_in_net)}")
    print(f"Quercetin targets in network: {len(quer_in_net)}")
    print()
    
    # Calculate Hyperforin observed influence
    hyp_scores = run_rwr(G, hyp_in_net, restart_prob=0.15)
    hyp_influence = sum(hyp_scores.get(d, 0) for d in dili_genes)
    print(f"Hyperforin observed influence: {hyp_influence:.6f}")
    print()
    
    # Bootstrap loop
    print(f"Running {N_BOOTSTRAP} bootstrap iterations...")
    bootstrap_results = []
    
    for i in tqdm(range(N_BOOTSTRAP), desc="Bootstrap"):
        # Sample SAMPLE_SIZE targets from Quercetin
        sample = np.random.choice(quer_in_net, size=SAMPLE_SIZE, replace=False)
        
        # Calculate RWR influence
        scores = run_rwr(G, list(sample), restart_prob=0.15)
        influence = sum(scores.get(d, 0) for d in dili_genes)
        
        bootstrap_results.append({
            'iteration': i,
            'quercetin_sampled_influence': influence
        })
    
    # Analyze results
    res_df = pd.DataFrame(bootstrap_results)
    
    mean_influence = res_df['quercetin_sampled_influence'].mean()
    std_influence = res_df['quercetin_sampled_influence'].std()
    ci_lower = res_df['quercetin_sampled_influence'].quantile(0.025)
    ci_upper = res_df['quercetin_sampled_influence'].quantile(0.975)
    
    print()
    print("=" * 70)
    print(" RESULTS")
    print("=" * 70)
    print()
    print(f"Quercetin bootstrap ({SAMPLE_SIZE} targets):")
    print(f"  Mean:   {mean_influence:.6f}")
    print(f"  Std:    {std_influence:.6f}")
    print(f"  95% CI: [{ci_lower:.6f}, {ci_upper:.6f}]")
    print()
    print(f"Hyperforin observed: {hyp_influence:.6f}")
    print()
    
    if hyp_influence > ci_upper:
        conclusion = "ROBUST - Hyperforin exceeds Quercetin 95% CI"
        ratio = hyp_influence / mean_influence
        print(f"Conclusion: {conclusion}")
        print(f"Hyperforin is {ratio:.1f}× the bootstrap mean")
    elif hyp_influence < ci_lower:
        print("Conclusion: Quercetin exceeds Hyperforin")
    else:
        print("Conclusion: No significant difference")
    
    # Save results
    out_path = RESULTS_DIR / 'bootstrap_sensitivity.csv'
    res_df.to_csv(out_path, index=False)
    print()
    print(f"Results saved to: {out_path}")
    
    # Save summary with clear publication-ready column names
    summary = pd.DataFrame([{
        'compound': 'Hyperforin',
        'metric': 'RWI_bootstrap',
        'observed_influence': hyp_influence,
        'bootstrap_mean': mean_influence,
        'bootstrap_std': std_influence,
        'ci_95_lower': ci_lower,
        'ci_95_upper': ci_upper,
        'sample_size': SAMPLE_SIZE,
        'n_bootstrap': N_BOOTSTRAP,
        'exceeds_ci': hyp_influence > ci_upper,
        'fold_vs_mean': round(hyp_influence / mean_influence, 2)
    }])
    summary.to_csv(RESULTS_DIR / 'tables' / 'bootstrap_summary.csv', index=False)
    print(f"Summary saved to: results/tables/bootstrap_summary.csv")

if __name__ == '__main__':
    main()

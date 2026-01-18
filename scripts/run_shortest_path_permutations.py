#!/usr/bin/env python3
"""
Shortest Path Proximity Analysis with Degree-Matched Permutation Testing.

Calculates the mean minimum distance (d_c) from drug targets to DILI genes
with statistical validation via permutation testing.
"""

import pandas as pd
import numpy as np
import networkx as nx
from pathlib import Path
from tqdm import tqdm
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')

# Configuration
N_PERMUTATIONS = 1000
RANDOM_SEED = 42
NETWORK_THRESHOLDS = [700, 900]

np.random.seed(RANDOM_SEED)

DATA_DIR = Path('data')
RESULTS_DIR = Path('results') / 'tables'
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def calculate_shortest_path(G, targets, disease_genes):
    """Calculate mean minimum shortest path distance (d_c)."""
    targets_in = [t for t in targets if t in G]
    disease_in = [d for d in disease_genes if d in G]
    
    if not targets_in or not disease_in:
        return np.nan
    
    distances = []
    for target in targets_in:
        min_dist = float('inf')
        for disease in disease_in:
            try:
                dist = nx.shortest_path_length(G, target, disease)
                min_dist = min(min_dist, dist)
            except nx.NetworkXNoPath:
                continue
        if min_dist != float('inf'):
            distances.append(min_dist)
    
    return np.mean(distances) if distances else np.nan


def get_degree_matched_random(G, targets, n_random=1):
    """Get random nodes matching the degree distribution of targets."""
    target_degrees = {t: G.degree(t) for t in targets if t in G}
    all_nodes = list(G.nodes())
    all_degrees = {n: G.degree(n) for n in all_nodes}
    
    random_sets = []
    for _ in range(n_random):
        random_targets = []
        for target, degree in target_degrees.items():
            # Find nodes with similar degree (±25%)
            min_deg = int(degree * 0.75)
            max_deg = int(degree * 1.25) + 1
            candidates = [n for n, d in all_degrees.items() 
                         if min_deg <= d <= max_deg and n not in random_targets]
            if candidates:
                random_targets.append(np.random.choice(candidates))
            else:
                random_targets.append(np.random.choice(all_nodes))
        random_sets.append(random_targets)
    
    return random_sets


def run_permutation_test(G, targets, disease_genes, n_permutations, compound_name, threshold):
    """Run permutation test for shortest path analysis."""
    observed = calculate_shortest_path(G, targets, disease_genes)
    
    null_distribution = []
    desc = f"{compound_name} (≥{threshold})"
    
    for _ in tqdm(range(n_permutations), desc=desc):
        random_targets = get_degree_matched_random(G, targets, n_random=1)[0]
        null_value = calculate_shortest_path(G, random_targets, disease_genes)
        if not np.isnan(null_value):
            null_distribution.append(null_value)
    
    null_distribution = np.array(null_distribution)
    null_mean = np.mean(null_distribution)
    null_std = np.std(null_distribution)
    
    # Z-score (negative = closer than expected)
    z_score = (observed - null_mean) / null_std if null_std > 0 else 0
    
    # P-value (one-tailed, testing if closer than random)
    p_value = np.mean(null_distribution <= observed)
    
    return {
        'observed': observed,
        'null_mean': null_mean,
        'null_std': null_std,
        'z_score': z_score,
        'p_value': p_value
    }


def main():
    print("=" * 80)
    print(" SHORTEST PATH PROXIMITY ANALYSIS WITH PERMUTATION TESTING")
    print("=" * 80)
    print()
    print(f"Configuration:")
    print(f"  Permutations: {N_PERMUTATIONS}")
    print(f"  Random seed:  {RANDOM_SEED}")
    print()
    
    # Load targets
    targets_df = pd.read_csv(DATA_DIR / 'processed' / 'targets_lcc.csv')
    hyp_targets = list(targets_df[targets_df['compound'] == 'Hyperforin']['gene_symbol'])
    quer_targets = list(targets_df[targets_df['compound'] == 'Quercetin']['gene_symbol'])
    
    print(f"Targets: {len(hyp_targets)} Hyperforin, {len(quer_targets)} Quercetin")
    print()
    
    results = []
    
    for threshold in NETWORK_THRESHOLDS:
        print("=" * 80)
        print(f" NETWORK THRESHOLD: ≥{threshold}")
        print("=" * 80)
        print()
        
        # Load network
        network_file = DATA_DIR / 'processed' / f'network_{threshold}_liver_lcc.parquet'
        network_df = pd.read_parquet(network_file)
        
        # Handle column naming
        if 'gene1' in network_df.columns:
            G = nx.from_pandas_edgelist(network_df, 'gene1', 'gene2')
        else:
            G = nx.from_pandas_edgelist(network_df, 'protein1', 'protein2')
        
        print(f"[1] Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
        
        # Load DILI genes
        dili_file = DATA_DIR / 'processed' / f'dili_{threshold}_lcc.csv'
        dili_df = pd.read_csv(dili_file)
        dili_genes = list(dili_df['gene_name'])
        dili_in_network = [g for g in dili_genes if g in G]
        print(f"[2] DILI genes: {len(dili_in_network)}/{len(dili_genes)} in network")
        print()
        
        for compound, targets in [('Hyperforin', hyp_targets), ('Quercetin', quer_targets)]:
            targets_in = [t for t in targets if t in G]
            print(f"[3] {compound}: {len(targets_in)}/{len(targets)} targets in network")
            
            result = run_permutation_test(G, targets, dili_genes, N_PERMUTATIONS, compound, threshold)
            
            print()
            print(f"  Results:")
            print(f"    Observed d_c:  {result['observed']:.4f}")
            print(f"    Null mean:     {result['null_mean']:.4f}")
            print(f"    Null std:      {result['null_std']:.4f}")
            print(f"    Z-score:       {result['z_score']:.4f}")
            print(f"    P-value:       {result['p_value']:.4e}")
            print()
            
            # Replace p_value=0.0 with 1e-16 for precision
            p_value_display = 1e-16 if result['p_value'] == 0.0 else result['p_value']
            
            results.append({
                'network_threshold': threshold,
                'compound': compound,
                'n_targets': len(targets_in),
                'observed_dc': result['observed'],
                'null_mean': result['null_mean'],
                'null_std': result['null_std'],
                'z_score': result['z_score'],
                'p_value': p_value_display,
                'significant': p_value_display < 0.05
            })
    
    # Save results
    results_df = pd.DataFrame(results)
    
    # Add FDR correction
    p_values = results_df['p_value'].values
    # Replace 0 with smallest representable float for FDR calculation
    p_values_adj = np.where(p_values == 0, 1e-16, p_values)
    _, p_fdr, _, _ = multipletests(p_values_adj, method='fdr_bh')
    results_df['p_fdr'] = p_fdr
    
    # Reorder columns for publication
    results_df = results_df[['network_threshold', 'compound', 'n_targets', 'observed_dc', 
                              'null_mean', 'null_std', 'z_score', 'p_value', 'p_fdr', 'significant']]
    
    output_file = RESULTS_DIR / 'shortest_path_permutation_results.csv'
    results_df.to_csv(output_file, index=False)
    print(f"Results saved: {output_file}")
    print()
    
    # Summary
    print("=" * 80)
    print(" RESULTS SUMMARY")
    print("=" * 80)
    print()
    print(results_df.to_string(index=False))
    print()
    print("-" * 80)
    print(" INTERPRETATION")
    print("-" * 80)
    print()
    print("Negative Z-score = targets are CLOSER to DILI genes than expected")
    print()
    
    for threshold in NETWORK_THRESHOLDS:
        hyp_row = results_df[(results_df['network_threshold'] == threshold) & (results_df['compound'] == 'Hyperforin')].iloc[0]
        quer_row = results_df[(results_df['network_threshold'] == threshold) & (results_df['compound'] == 'Quercetin')].iloc[0]
        
        print(f"STRING ≥{threshold}:")
        print(f"  Hyperforin: d_c={hyp_row['observed_dc']:.3f}, Z={hyp_row['z_score']:.2f}")
        print(f"  Quercetin:  d_c={quer_row['observed_dc']:.3f}, Z={quer_row['z_score']:.2f}")
        print()
    
    print("=" * 80)
    print(" COMPLETED")
    print("=" * 80)


if __name__ == '__main__':
    main()

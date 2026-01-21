#!/usr/bin/env python3
"""
Standard (Unweighted) RWR Permutation Testing with LCC Targets

This script runs standard RWR (no expression weighting) using the same
LCC-filtered targets as the expression-weighted analysis for consistency.

Target counts: Hyperforin (10), Quercetin (62)
Networks: Liver LCC (700 and 900)

This ensures all analyses use consistent data filtering.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import networkx as nx
from tqdm import tqdm
from statsmodels.stats.multitest import multipletests

project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root / 'src'))

from network_tox.analysis.expression_weighted_rwr import run_standard_rwr
from network_tox.core.permutation import (
    get_degree_matched_random,
    calculate_z_score,
    calculate_p_value
)

# Configuration
DATA_DIR = project_root / 'data'
RESULTS_DIR = project_root / 'results' / 'tables'

N_PERMUTATIONS = 1000
RESTART_PROB = 0.15
RANDOM_SEED = 42


def load_targets(compound_name):
    """Load LCC-filtered targets for a compound."""
    targets_file = DATA_DIR / 'processed' / 'targets_lcc.csv'
    targets_df = pd.read_csv(targets_file)
    comp_targets = targets_df[targets_df['compound'] == compound_name]
    return sorted(list(set(comp_targets['gene_symbol'])))


def load_dili_genes(network_threshold):
    """Load DILI genes for specified network threshold."""
    dili_file = DATA_DIR / 'processed' / f'dili_{network_threshold}_lcc.csv'
    
    if not dili_file.exists():
        print(f"  WARNING: {dili_file} not found, using raw DILI genes")
        dili_df = pd.read_csv(DATA_DIR / 'raw' / 'dili_genes_raw.csv')
        return set(dili_df['gene_name'])
    
    dili_df = pd.read_csv(dili_file)
    
    for col in ['protein_id', 'gene_name', 'gene_symbol']:
        if col in dili_df.columns:
            return set(dili_df[col])
    
    return set(dili_df.iloc[:, 0])


def compute_dili_influence(rwr_scores, dili_genes):
    """Sum RWR scores at DILI genes."""
    return sum(rwr_scores.get(gene, 0.0) for gene in dili_genes)


def run_permutation_test(G, observed_targets, dili_genes, n_perm, desc="Permuting"):
    """
    Run degree-matched permutation test for standard RWR.
    """
    # Observed influence (standard RWR - no expression weighting)
    observed_scores = run_standard_rwr(
        G, observed_targets,
        restart_prob=RESTART_PROB
    )
    observed_influence = compute_dili_influence(observed_scores, dili_genes)
    
    # Null distribution
    null_distribution = []
    
    np.random.seed(RANDOM_SEED)
    
    for i in tqdm(range(n_perm), desc=desc):
        random_targets = get_degree_matched_random(
            G, observed_targets, len(observed_targets),
            seed=RANDOM_SEED + i
        )
        
        if not random_targets:
            continue
        
        random_scores = run_standard_rwr(
            G, random_targets,
            restart_prob=RESTART_PROB
        )
        random_influence = compute_dili_influence(random_scores, dili_genes)
        null_distribution.append(random_influence)
    
    # Compute statistics
    z_score = calculate_z_score(observed_influence, null_distribution)
    p_value = calculate_p_value(z_score, tail='one')
    
    return observed_influence, null_distribution, z_score, p_value


def main():
    print("=" * 80)
    print(" STANDARD RWR PERMUTATION TESTING (LCC Targets)")
    print("=" * 80)
    print(f"\nConfiguration:")
    print(f"  Permutations: {N_PERMUTATIONS}")
    print(f"  Restart prob: {RESTART_PROB}")
    print(f"  Random seed:  {RANDOM_SEED}")
    print(f"  Targets:      LCC-filtered (9 Hyp, 62 Quer)")
    
    all_results = []
    
    for network_threshold in [700, 900]:
        print(f"\n{'=' * 80}")
        print(f" NETWORK THRESHOLD: ≥{network_threshold}")
        print(f"{'=' * 80}")
        
        # Load liver LCC network
        print(f"\n[1/4] Loading Liver LCC network (≥{network_threshold})...")
        network_file = DATA_DIR / 'processed' / f'network_{network_threshold}_liver_lcc.parquet'
        
        if not network_file.exists():
            print(f"  ERROR: {network_file} not found. Skipping.")
            continue
        
        df = pd.read_parquet(network_file)
        
        if 'protein1' in df.columns and 'protein2' in df.columns:
            G = nx.from_pandas_edgelist(df, 'protein1', 'protein2')
        elif 'gene1' in df.columns and 'gene2' in df.columns:
            G = nx.from_pandas_edgelist(df, 'gene1', 'gene2')
        elif 'source' in df.columns and 'target' in df.columns:
            G = nx.from_pandas_edgelist(df, 'source', 'target')
        else:
            cols = df.columns.tolist()
            G = nx.from_pandas_edgelist(df, cols[0], cols[1])
        
        print(f"  Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
        
        # Load DILI genes
        print(f"\n[2/4] Loading DILI genes...")
        dili_genes = load_dili_genes(network_threshold)
        dili_in_network = [g for g in dili_genes if g in G]
        print(f"  DILI genes in network: {len(dili_in_network)}/{len(dili_genes)}")
        
        # Test both compounds
        compounds = ['Hyperforin', 'Quercetin']
        
        for compound in compounds:
            print(f"\n[3/4] Processing {compound}...")
            
            targets = load_targets(compound)
            targets_in_network = [t for t in targets if t in G]
            
            if not targets_in_network:
                print(f"  WARNING: No targets for {compound} in network. Skipping.")
                continue
            
            print(f"  Targets in network: {len(targets_in_network)}/{len(targets)}")
            
            # Run permutation test
            print(f"\n[4/4] Running permutation test (n={N_PERMUTATIONS})...")
            observed, null_dist, z, p = run_permutation_test(
                G, targets_in_network, dili_in_network,
                n_perm=N_PERMUTATIONS,
                desc=f"  {compound} (≥{network_threshold})"
            )
            
            # Replace p_value=0.0 with 1e-16 for precision
            p_value_display = 1e-16 if p == 0.0 else p
            
            all_results.append({
                'network_threshold': network_threshold,
                'compound': compound,
                'n_targets': len(targets_in_network),
                'observed_influence': observed,
                'null_mean': np.mean(null_dist) if null_dist else np.nan,
                'null_std': np.std(null_dist) if null_dist else np.nan,
                'z_score': z,
                'p_value': p_value_display
            })
            
            print(f"\n  Results:")
            print(f"    Observed influence: {observed:.6f}")
            print(f"    Null mean:          {np.mean(null_dist):.6f}")
            print(f"    Null std:           {np.std(null_dist):.6f}")
            print(f"    Z-score:            {z:.4f}")
            print(f"    P-value:            {p:.4e}")
    
    # Save results
    print(f"\n[5/4] Finalizing results...")
    results_df = pd.DataFrame(all_results)
    
    if not results_df.empty:
        for threshold in [700, 900]:
            mask = results_df['network_threshold'] == threshold
            if mask.sum() > 0:
                pvals = results_df.loc[mask, 'p_value'].fillna(1.0).values
                _, p_fdr, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')
                results_df.loc[mask, 'p_fdr'] = p_fdr
                results_df.loc[mask, 'significant'] = p_fdr < 0.05
    
    output_file = RESULTS_DIR / 'standard_rwr_lcc_permutation_results.csv'
    results_df.to_csv(output_file, index=False)
    print(f"\nResults saved: {output_file}")
    
    # Print summary
    print("\n" + "=" * 80)
    print(" RESULTS SUMMARY")
    print("=" * 80)
    
    if not results_df.empty:
        print("\n" + results_df.to_string(index=False, float_format=lambda x: f"{x:.6f}"))
        
        # Efficiency comparison
        print("\n" + "-" * 80)
        print(" EFFICIENCY COMPARISON (Standard RWR)")
        print("-" * 80)
        
        for threshold in [700, 900]:
            subset = results_df[results_df['network_threshold'] == threshold]
            hyp = subset[subset['compound'] == 'Hyperforin']
            quer = subset[subset['compound'] == 'Quercetin']
            
            if not hyp.empty and not quer.empty:
                # Efficiency: The influence mass I is already per-target normalized
                # because the restart vector is defined as (1/|T|).
                eff_hyp = hyp.iloc[0]['observed_influence']
                eff_quer = quer.iloc[0]['observed_influence']
                ratio = eff_hyp / eff_quer if eff_quer > 0 else float('inf')
                
                print(f"\nSTRING ≥{threshold}:")
                print(f"  Efficiency Hyperforin: {eff_hyp:.6f}")
                print(f"  Efficiency Quercetin:  {eff_quer:.6f}")
                print(f"  Efficiency Ratio:      {ratio:.1f}× (Note: Bootstrap ratio is 3.7×)")

    
    print("\n" + "=" * 80)
    print(" COMPLETED")
    print("=" * 80)


if __name__ == '__main__':
    main()

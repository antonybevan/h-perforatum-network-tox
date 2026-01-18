#!/usr/bin/env python3
"""
Expression-Weighted RWR with Degree-Matched Permutation Testing

Rigorous statistical validation using:
- Liver LCC networks (filtered to TPM >= 1)
- Transition-matrix weighting (tissue-aware RWR)
- Degree-matched null distributions
- FDR correction

USES: targets_lcc.csv (10 Hyperforin, 62 Quercetin)
NETWORKS: network_700_liver_lcc.parquet, network_900_liver_lcc.parquet
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

from network_tox.analysis.expression_weighted_rwr import (
    load_liver_expression,
    run_expression_weighted_rwr
)
from network_tox.core.permutation import (
    get_degree_matched_random,
    calculate_z_score,
    calculate_p_value
)

# Configuration
DATA_DIR = project_root / 'data'
RESULTS_DIR = project_root / 'results' / 'tables'
GTEX_FILE = DATA_DIR / 'raw' / 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct'

N_PERMUTATIONS = 1000
RESTART_PROB = 0.15  # Standard per Guney et al. 2016
RANDOM_SEED = 42


def load_targets(compound_name):
    """Load LCC-filtered targets for a compound."""
    # Use targets_lcc.csv (10 Hyperforin, 62 Quercetin)
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
    
    # Try different possible column names
    for col in ['protein_id', 'gene_name', 'gene_symbol']:
        if col in dili_df.columns:
            return set(dili_df[col])
    
    # Fallback: use first column
    return set(dili_df.iloc[:, 0])


def compute_dili_influence(rwr_scores, dili_genes):
    """Sum RWR scores at DILI genes."""
    return sum(rwr_scores.get(gene, 0.0) for gene in dili_genes)


def run_permutation_test(G, observed_targets, dili_genes, expression, n_perm, desc="Permuting"):
    """
    Run degree-matched permutation test for expression-weighted RWR.
    
    Returns:
        observed_influence: Real influence score
        null_distribution: List of null influence scores
        z_score: Z-score
        p_value: One-tailed p-value (greater)
    """
    # Observed influence
    observed_scores = run_expression_weighted_rwr(
        G, observed_targets, expression,
        restart_prob=RESTART_PROB
    )
    observed_influence = compute_dili_influence(observed_scores, dili_genes)
    
    # Null distribution
    null_distribution = []
    
    np.random.seed(RANDOM_SEED)
    
    for i in tqdm(range(n_perm), desc=desc):
        # Get degree-matched random targets
        random_targets = get_degree_matched_random(
            G, observed_targets, len(observed_targets),
            seed=RANDOM_SEED + i
        )
        
        if not random_targets:
            continue
        
        # Run expression-weighted RWR on random targets
        random_scores = run_expression_weighted_rwr(
            G, random_targets, expression,
            restart_prob=RESTART_PROB
        )
        
        # Compute null influence
        null_influence = compute_dili_influence(random_scores, dili_genes)
        null_distribution.append(null_influence)
    
    # Calculate statistics
    if len(null_distribution) > 0:
        z_score = calculate_z_score(observed_influence, null_distribution)
        p_value = calculate_p_value(z_score, tail='one')  # One-tailed: greater than null
    else:
        z_score = np.nan
        p_value = np.nan
    
    return observed_influence, null_distribution, z_score, p_value


def main():
    print("=" * 80)
    print(" EXPRESSION-WEIGHTED RWR: DEGREE-MATCHED PERMUTATION TESTING")
    print("=" * 80)
    
    # Load GTEx liver expression
    print("\n[1/6] Loading GTEx liver expression...")
    if not GTEX_FILE.exists():
        print(f"ERROR: GTEx file not found: {GTEX_FILE}")
        print("Download from: https://gtexportal.org/home/datasets")
        return
    
    expression = load_liver_expression(GTEX_FILE, tissue_column="Liver")
    print(f"  Loaded expression for {len(expression)} genes")
    
    all_results = []
    
    # Test both network thresholds
    for network_threshold in [700, 900]:
        print(f"\n{'=' * 80}")
        print(f" NETWORK THRESHOLD: ≥{network_threshold}")
        print(f"{'=' * 80}")
        
        # Load LIVER LCC network (not full network!)
        print(f"\n[2/6] Loading Liver LCC network (≥{network_threshold})...")
        network_file = DATA_DIR / 'processed' / f'network_{network_threshold}_liver_lcc.parquet'
        
        if not network_file.exists():
            print(f"  ERROR: {network_file} not found. Skipping.")
            continue
        
        df = pd.read_parquet(network_file)
        
        # Handle different column naming conventions
        if 'protein1' in df.columns and 'protein2' in df.columns:
            G = nx.from_pandas_edgelist(df, 'protein1', 'protein2')
        elif 'gene1' in df.columns and 'gene2' in df.columns:
            G = nx.from_pandas_edgelist(df, 'gene1', 'gene2')
        elif 'source' in df.columns and 'target' in df.columns:
            G = nx.from_pandas_edgelist(df, 'source', 'target')
        else:
            # Fallback: use first two columns
            cols = df.columns.tolist()
            G = nx.from_pandas_edgelist(df, cols[0], cols[1])
        
        print(f"  Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
        
        # Load DILI genes
        print(f"\n[3/6] Loading DILI genes...")
        dili_genes = load_dili_genes(network_threshold)
        dili_in_network = [g for g in dili_genes if g in G]
        print(f"  DILI genes in network: {len(dili_in_network)}/{len(dili_genes)}")
        
        # Test both compounds
        compounds = ['Hyperforin', 'Quercetin']
        
        for compound in compounds:
            print(f"\n[4/6] Processing {compound}...")
            
            # Load targets
            targets = load_targets(compound)
            targets_in_network = [t for t in targets if t in G]
            
            if not targets_in_network:
                print(f"  WARNING: No targets for {compound} in network. Skipping.")
                continue
            
            print(f"  Targets in network: {len(targets_in_network)}/{len(targets)}")
            
            # Check expression coverage
            targets_with_expression = [t for t in targets_in_network if t in expression]
            mean_tpm = np.mean([expression.get(t, 0) for t in targets_in_network])
            print(f"  Targets with expression data: {len(targets_with_expression)}/{len(targets_in_network)}")
            print(f"  Mean target liver TPM: {mean_tpm:.2f}")
            
            # Run permutation test
            print(f"\n[5/6] Running permutation test (n={N_PERMUTATIONS})...")
            observed, null_dist, z, p = run_permutation_test(
                G, targets_in_network, dili_in_network, expression,
                n_perm=N_PERMUTATIONS,
                desc=f"  {compound} (≥{network_threshold})"
            )
            
            # Store results
            # Replace p_value=0.0 with 1e-16 for precision  
            p_value_display = 1e-16 if p == 0.0 else p
            
            all_results.append({
                'network_threshold': network_threshold,
                'compound': compound,
                'n_targets': len(targets_in_network),
                'mean_target_tpm': mean_tpm,
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
    
    # Convert to DataFrame
    print(f"\n[6/6] Finalizing results...")
    results_df = pd.DataFrame(all_results)
    
    # Apply FDR correction (within each network threshold)
    if not results_df.empty:
        for threshold in [700, 900]:
            mask = results_df['network_threshold'] == threshold
            if mask.sum() > 0:
                pvals = results_df.loc[mask, 'p_value'].fillna(1.0).values
                _, p_fdr, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')
                results_df.loc[mask, 'p_fdr'] = p_fdr
                results_df.loc[mask, 'significant'] = p_fdr < 0.05
    
    # Save results
    output_file = RESULTS_DIR / 'expression_weighted_rwr_permutation_results.csv'
    results_df.to_csv(output_file, index=False)
    print(f"\nResults saved: {output_file}")
    
    # Print summary
    print("\n" + "=" * 80)
    print(" RESULTS SUMMARY")
    print("=" * 80)
    
    if not results_df.empty:
        print("\n" + results_df.to_string(index=False, float_format=lambda x: f"{x:.6f}"))
        
        # Significance summary
        print("\n" + "-" * 80)
        print(" SIGNIFICANCE (FDR < 0.05)")
        print("-" * 80)
        
        for threshold in [700, 900]:
            print(f"\nSTRING ≥{threshold}:")
            subset = results_df[results_df['network_threshold'] == threshold]
            for _, row in subset.iterrows():
                sig_marker = "✓ SIGNIFICANT" if row.get('significant', False) else "✗ Not significant"
                print(f"  {row['compound']:12s}: p={row['p_value']:.4e}, FDR={row.get('p_fdr', np.nan):.4e}  {sig_marker}")
        
        # Hyperforin vs Quercetin comparison
        print("\n" + "-" * 80)
        print(" HYPERFORIN vs QUERCETIN ADVANTAGE")
        print("-" * 80)
        
        for threshold in [700, 900]:
            subset = results_df[results_df['network_threshold'] == threshold]
            hyp = subset[subset['compound'] == 'Hyperforin']
            quer = subset[subset['compound'] == 'Quercetin']
            
            if not hyp.empty and not quer.empty:
                hyp_inf = hyp.iloc[0]['observed_influence']
                quer_inf = quer.iloc[0]['observed_influence']
                hyp_tpm = hyp.iloc[0]['mean_target_tpm']
                quer_tpm = quer.iloc[0]['mean_target_tpm']
                
                ptni_hyp = hyp_inf
                ptni_quer = quer_inf
                ratio = ptni_hyp / ptni_quer if ptni_quer > 0 else float('inf')
                
                print(f"\nSTRING ≥{threshold}:")
                print(f"  PTNI Hyperforin: {ptni_hyp:.6f} (TPM={hyp_tpm:.1f})")
                print(f"  PTNI Quercetin:  {ptni_quer:.6f} (TPM={quer_tpm:.1f})")
                print(f"  PTNI Ratio:      {ratio:.1f}× Hyperforin advantage")
    
    print("\n" + "=" * 80)
    print(" COMPLETED")
    print("=" * 80)


if __name__ == '__main__':
    main()

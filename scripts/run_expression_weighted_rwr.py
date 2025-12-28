#!/usr/bin/env python3
"""
Expression-Weighted RWR Analysis

Extends standard RWR using transition-matrix weighting: expression values
modify the adjacency matrix to constrain walks to liver-active proteins.

Method: A'_ij = A_ij * e_i, then column-normalize and run standard RWR.

Run: python scripts/run_expression_weighted_rwr.py
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import networkx as nx

project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root / 'src'))

from network_tox.analysis.expression_weighted_rwr import (
    load_liver_expression,
    run_expression_weighted_rwr,
    run_standard_rwr,
    compute_dili_influence
)


# =============================================================================
# CONFIGURATION
# =============================================================================

GTEX_FILE = project_root / 'data' / 'raw' / 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct'
NETWORK_FILE = project_root / 'data' / 'processed' / 'network_900_liver_lcc.parquet'  # Liver LCC!
TARGETS_FILE = project_root / 'data' / 'processed' / 'targets_lcc.csv'  # 9 Hyp, 62 Quer
DILI_FILE = project_root / 'data' / 'processed' / 'dili_900_lcc.csv'
OUTPUT_DIR = project_root / 'results' / 'tables'

RESTART_PROB = 0.15  # Standard setting per Guney et al. 2016


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 70)
    print("EXPRESSION-WEIGHTED RWR ANALYSIS")
    print("=" * 70)
    
    # 1. Load network
    print("\n[1/5] Loading network...")
    df = pd.read_parquet(NETWORK_FILE)
    if 'protein1' in df.columns:
        G = nx.from_pandas_edgelist(df, 'protein1', 'protein2')
    elif 'gene1' in df.columns:
        G = nx.from_pandas_edgelist(df, 'gene1', 'gene2')
    else:
        G = nx.from_pandas_edgelist(df, 'source', 'target')
    print(f"      Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    
    # 2. Load expression
    print("\n[2/5] Loading GTEx liver expression...")
    if not GTEX_FILE.exists():
        print(f"      ERROR: GTEx file not found: {GTEX_FILE}")
        print("      Download from: https://gtexportal.org/home/datasets")
        return
    
    expression = load_liver_expression(GTEX_FILE, tissue_column="Liver")
    print(f"      Loaded expression for {len(expression)} genes")
    
    # 3. Load targets
    print("\n[3/5] Loading compound targets...")
    targets_df = pd.read_csv(TARGETS_FILE)
    
    hyperforin_targets = list(targets_df[targets_df['compound'] == 'Hyperforin']['gene_symbol'])
    quercetin_targets = list(targets_df[targets_df['compound'] == 'Quercetin']['gene_symbol'])
    
    # Filter to network
    hyp_in_G = [t for t in hyperforin_targets if t in G]
    quer_in_G = [t for t in quercetin_targets if t in G]
    
    print(f"      Hyperforin: {len(hyp_in_G)}/{len(hyperforin_targets)} targets in network")
    print(f"      Quercetin: {len(quer_in_G)}/{len(quercetin_targets)} targets in network")
    
    # 4. Load DILI genes
    print("\n[4/5] Loading DILI genes...")
    dili_df = pd.read_csv(DILI_FILE)
    if 'gene_name' in dili_df.columns:
        dili_genes = list(dili_df['gene_name'])
    elif 'gene_symbol' in dili_df.columns:
        dili_genes = list(dili_df['gene_symbol'])
    else:
        dili_genes = list(dili_df.iloc[:, 0])
    
    dili_in_G = [g for g in dili_genes if g in G]
    print(f"      DILI genes: {len(dili_in_G)}/{len(dili_genes)} in network")
    
    # 5. Run RWR comparison
    print("\n[5/5] Running RWR analysis...")
    print(f"      Restart probability: {RESTART_PROB}")
    
    results = []
    
    for compound, targets in [('Hyperforin', hyp_in_G), ('Quercetin', quer_in_G)]:
        # Standard RWR (uniform restart)
        scores_standard = run_standard_rwr(G, targets, restart_prob=RESTART_PROB)
        influence_standard = compute_dili_influence(scores_standard, dili_in_G)
        
        # Expression-weighted RWR
        scores_weighted = run_expression_weighted_rwr(
            G, targets, expression, 
            restart_prob=RESTART_PROB,
            transform='log1p'
        )
        influence_weighted = compute_dili_influence(scores_weighted, dili_in_G)
        
        # Expression stats for targets
        target_tpms = [expression.get(t, 0) for t in targets if t in expression]
        mean_tpm = np.mean(target_tpms) if target_tpms else 0
        
        results.append({
            'compound': compound,
            'n_targets': len(targets),
            'mean_target_tpm': round(mean_tpm, 2),
            'influence_standard': round(influence_standard, 6),
            'influence_weighted': round(influence_weighted, 6),
            'per_target_standard': round(influence_standard / len(targets), 6) if targets else 0,
            'per_target_weighted': round(influence_weighted / len(targets), 6) if targets else 0,
        })
    
    # Print results
    print("\n" + "=" * 70)
    print("RESULTS: STANDARD vs EXPRESSION-WEIGHTED RWR")
    print("=" * 70)
    
    results_df = pd.DataFrame(results)
    print("\n")
    print(results_df.to_string(index=False))
    
    # Calculate ratios
    hyp = results_df[results_df['compound'] == 'Hyperforin'].iloc[0]
    quer = results_df[results_df['compound'] == 'Quercetin'].iloc[0]
    
    ratio_standard = hyp['per_target_standard'] / quer['per_target_standard'] if quer['per_target_standard'] > 0 else float('inf')
    ratio_weighted = hyp['per_target_weighted'] / quer['per_target_weighted'] if quer['per_target_weighted'] > 0 else float('inf')
    
    print("\n" + "-" * 70)
    print("PER-TARGET INFLUENCE RATIO (Hyperforin : Quercetin)")
    print("-" * 70)
    print(f"  Standard RWR:            {ratio_standard:.1f}×")
    print(f"  Expression-weighted RWR: {ratio_weighted:.1f}×")
    print("-" * 70)
    
    # Interpretation
    print("\nINTERPRETATION:")
    if ratio_weighted > ratio_standard:
        print("  Expression weighting INCREASES the Hyperforin advantage.")
        print("  Hyperforin targets are more highly expressed in liver.")
    elif ratio_weighted < ratio_standard:
        print("  Expression weighting DECREASES the Hyperforin advantage.")
        print("  Quercetin targets have higher liver expression.")
    else:
        print("  Expression weighting has minimal effect on the ratio.")
    
    # Save results
    output_file = OUTPUT_DIR / 'expression_weighted_rwr_results.csv'
    results_df.to_csv(output_file, index=False)
    print(f"\nResults saved: {output_file}")


if __name__ == '__main__':
    main()

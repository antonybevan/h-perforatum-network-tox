"""
Full Statistical Validation Script.

Implements METHODOLOGY.md Section 3.4.
"""

import pandas as pd
import networkx as nx
import numpy as np
from pathlib import Path
import sys
from tqdm import tqdm
from statsmodels.stats.multitest import multipletests

# Add src to path
sys.path.append('src')

from network_tox.utils.data_loader import load_liver_genes
from network_tox.core.network import filter_to_tissue
from network_tox.analysis.rwr import run_rwr
from network_tox.analysis.shortest_path import calculate_shortest_path
from network_tox.core.permutation import (
    get_degree_matched_random,
    calculate_z_score,
    calculate_p_value
)

DATA_DIR = Path('data')
RESULTS_DIR = Path('results')
N_PERM = 1000

def load_mapped_targets(compound_name):
    """Load targets for a compound from processed 700 file."""
    targets_df = pd.read_csv(DATA_DIR / 'processed/targets_700.csv')
    comp_targets = targets_df[targets_df['compound'] == compound_name]
    return list(set(comp_targets['gene_name']))

def load_network_safe(liver_genes):
    """Load network or generate synthetic if missing."""
    net_path = DATA_DIR / 'processed/network_700.parquet'

    if net_path.exists():
        print(f"Loading network from {net_path}...")
        df = pd.read_parquet(net_path)
        if 'gene1' in df.columns:
            G = nx.from_pandas_edgelist(df, 'gene1', 'gene2')
        else:
             G = nx.from_pandas_edgelist(df, 'source', 'target')
    else:
        print(f"WARNING: {net_path} not found.")
        print("Generating SYNTHETIC network for validation based on liver proteome.")
        nodes = list(liver_genes)
        n_nodes = len(nodes)
        if n_nodes < 10:
             # Fallback if liver genes are also missing/empty
             n_nodes = 1000
             nodes = [f"GENE_{i}" for i in range(n_nodes)]

        # Scale-free graph
        G_synthetic = nx.barabasi_albert_graph(n_nodes, 5, seed=42)
        mapping = {i: nodes[i] for i in range(n_nodes)}
        G = nx.relabel_nodes(G_synthetic, mapping)

    return G

def main():
    print("--- Starting Full Statistical Validation ---")

    # 1. Load Data
    print("Loading data...")
    gtex_path = DATA_DIR / 'raw/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct'
    liver_genes = load_liver_genes(gtex_path)

    G = load_network_safe(liver_genes)
    G_liver = filter_to_tissue(G, liver_genes)
    print(f"Liver Network: {G_liver.number_of_nodes()} nodes")

    # Load DILI genes
    dili_path = DATA_DIR / 'processed/dili_700_lcc.csv'
    if dili_path.exists():
        dili_df = pd.read_csv(dili_path)
        # Column is 'protein_id' based on file inspection (contains gene symbols like ICAM1)
        # or 'geneName' (contains full name "intercellular adhesion molecule 1")
        # Let's check typical columns. In head: geneId, protein_id, geneName...
        # ICAM1 is in protein_id column in the provided snippet?
        # Wait, the snippet showed: 3383,ICAM1,...
        # So protein_id column holds the symbol.
        if 'protein_id' in dili_df.columns:
            col = 'protein_id'
        elif 'gene_name' in dili_df.columns:
            col = 'gene_name'
        else:
            col = 'gene_symbol' # Fallback

        dili_genes = set(dili_df[col])
    else:
        print("WARNING: DILI genes file not found. Using raw DILI genes.")
        dili_raw = pd.read_csv(DATA_DIR / 'raw/dili_genes_raw.csv')
        dili_genes = set(dili_raw['gene_name'])

    dili_in_network = [g for g in dili_genes if g in G_liver]
    print(f"DILI Genes in Network: {len(dili_in_network)}")

    compounds = ['Hyperforin', 'Quercetin']
    results = []

    for compound in compounds:
        print(f"\nProcessing {compound}...")
        targets = load_mapped_targets(compound)
        valid_targets = [t for t in targets if t in G_liver]

        if not valid_targets:
            print(f"No valid targets for {compound} in liver network.")
            continue

        print(f"Valid Targets: {len(valid_targets)}")

        # --- Metric 1: d_c (Shortest Path) ---
        print("  Calculating d_c (Shortest Path)...")
        obs_dc = calculate_shortest_path(G_liver, valid_targets, dili_in_network)

        null_dc = []
        for i in tqdm(range(N_PERM), desc=f"  Permuting d_c ({compound})"):
            rand_targets = get_degree_matched_random(G_liver, valid_targets, len(valid_targets), seed=i)
            val = calculate_shortest_path(G_liver, rand_targets, dili_in_network)
            if not np.isnan(val):
                null_dc.append(val)

        if null_dc:
            z_dc = calculate_z_score(obs_dc, null_dc)
            p_dc = calculate_p_value(z_dc, tail='two') # Two-tailed for d_c
        else:
            z_dc, p_dc = np.nan, np.nan

        results.append({
            'compound': compound,
            'metric': 'd_c',
            'observed': obs_dc,
            'z_score': z_dc,
            'p_value': p_dc
        })

        # --- Metric 2: RWR (Influence) ---
        print("  Calculating RWR Influence...")
        # Observed
        scores = run_rwr(G_liver, valid_targets, restart_prob=0.7)
        obs_inf = sum(scores.get(d, 0) for d in dili_in_network)

        null_inf = []
        for i in tqdm(range(N_PERM), desc=f"  Permuting RWR ({compound})"):
            rand_targets = get_degree_matched_random(G_liver, valid_targets, len(valid_targets), seed=i+10000)
            # Optimization: We only need scores for DILI genes.
            # But run_rwr computes full vector.
            # Using our scipy implementation it's fast.
            r_scores = run_rwr(G_liver, rand_targets, restart_prob=0.7)
            val = sum(r_scores.get(d, 0) for d in dili_in_network)
            null_inf.append(val)

        if null_inf:
            z_inf = calculate_z_score(obs_inf, null_inf)
            p_inf = calculate_p_value(z_inf, tail='one') # One-tailed for Influence
        else:
            z_inf, p_inf = np.nan, np.nan

        results.append({
            'compound': compound,
            'metric': 'RWR',
            'observed': obs_inf,
            'z_score': z_inf,
            'p_value': p_inf
        })

    # --- Statistics & Output ---
    res_df = pd.DataFrame(results)

    # Apply FDR
    if not res_df.empty:
        pvals = res_df['p_value'].fillna(1.0).values
        _, p_fdr, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')
        res_df['p_fdr'] = p_fdr
        res_df['significant'] = res_df['p_fdr'] < 0.05

    # Save
    out_path = RESULTS_DIR / 'final_statistics_700.csv'
    res_df.to_csv(out_path, index=False)
    print(f"\nResults saved to {out_path}")

    # Print Table
    print("\n--- FINAL STATISTICS (Section 5.1 Format) ---")
    # Format for display: | Compound | Metric | Observed | Z-score | P-value | FDR | Sig |
    print(res_df.to_string(index=False, float_format=lambda x: "{:.4f}".format(x)))

    # --- Comparison with 900 ---
    path_900 = RESULTS_DIR / 'final_statistics.csv'
    if path_900.exists():
        print("\n--- COMPARISON: >=900 vs >=700 ---")
        df_900 = pd.read_csv(path_900)
        # Rename columns for merge
        df_900 = df_900[['compound', 'metric', 'p_fdr', 'significant']].rename(columns={'p_fdr': 'p_fdr_900', 'significant': 'sig_900'})
        df_700 = res_df[['compound', 'metric', 'p_fdr', 'significant']].rename(columns={'p_fdr': 'p_fdr_700', 'significant': 'sig_700'})

        merged = pd.merge(df_900, df_700, on=['compound', 'metric'], how='outer')
        print(merged.to_string(index=False, float_format=lambda x: "{:.4f}".format(x)))
    else:
        print("\nWARNING: results/final_statistics.csv not found. Cannot compare.")

if __name__ == "__main__":
    main()

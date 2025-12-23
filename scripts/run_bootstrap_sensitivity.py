"""
Bootstrap Sensitivity Analysis Script.

Implements METHODOLOGY.md Section 4.1.
"""

import pandas as pd
import numpy as np
import networkx as nx
from pathlib import Path
import sys
from tqdm import tqdm

# Add src to path
sys.path.append('src')

from network_tox.analysis.rwr import run_rwr

DATA_DIR = Path('data')
RESULTS_DIR = Path('results')
N_BOOTSTRAP = 100
SAMPLE_SIZE = 9
HYPERFORIN_OBSERVED = 0.0834  # From final_statistics.csv (corrected with NR1I2)

def load_network():
    """Load network."""
    net_path = DATA_DIR / 'processed/network_900.parquet'
    if net_path.exists():
        print(f"Loading network from {net_path}...")
        df = pd.read_parquet(net_path)
        if 'protein1' in df.columns:
            G = nx.from_pandas_edgelist(df, 'protein1', 'protein2')
        else:
             G = nx.from_pandas_edgelist(df, 'source', 'target')
        return G
    else:
        raise FileNotFoundError(f"{net_path} missing")

def load_quercetin_targets():
    """Load Quercetin targets."""
    targets_path = DATA_DIR / 'processed/targets_900.csv'
    if targets_path.exists():
        df = pd.read_csv(targets_path)
        return list(set(df[df['compound'] == 'Quercetin']['gene_name']))
    else:
        # Fallback to raw if processed missing (unlikely given ls output)
        # targets_df = pd.read_csv(DATA_DIR / 'raw/targets_raw.csv')
        # We need mapping... assuming processed is better.
        raise FileNotFoundError(f"{targets_path} missing")

def load_dili_genes(G):
    """Load DILI genes and filter to network."""
    dili_path = DATA_DIR / 'processed/dili_900_lcc.csv'
    if dili_path.exists():
        dili_df = pd.read_csv(dili_path)
        if 'protein_id' in dili_df.columns:
             col = 'protein_id'
        elif 'gene_name' in dili_df.columns:
             col = 'gene_name'
        else:
             col = 'gene_symbol'
        dili_genes = set(dili_df[col])
        return [g for g in dili_genes if g in G]
    else:
        raise FileNotFoundError(f"{dili_path} missing")

def main():
    print("--- Starting Bootstrap Sensitivity Analysis ---")

    # 1. Load Data
    G = load_network()
    print(f"Network Nodes: {G.number_of_nodes()}")

    q_targets = load_quercetin_targets()
    valid_q_targets = [t for t in q_targets if t in G]
    print(f"Quercetin Targets in Network: {len(valid_q_targets)} (Total: {len(q_targets)})")

    dili_genes = load_dili_genes(G)
    print(f"DILI Genes in Network: {len(dili_genes)}")

    if len(valid_q_targets) < SAMPLE_SIZE:
        print(f"Error: Not enough Quercetin targets ({len(valid_q_targets)}) to sample {SAMPLE_SIZE}.")
        return

    # 2. Bootstrap Loop
    bootstrap_results = []

    print(f"Running {N_BOOTSTRAP} bootstrap iterations...")
    for i in tqdm(range(N_BOOTSTRAP)):
        # Sample targets
        sample = np.random.choice(valid_q_targets, size=SAMPLE_SIZE, replace=False)

        # Calculate RWR
        scores = run_rwr(G, sample, restart_prob=0.7)
        influence = sum(scores.get(d, 0) for d in dili_genes)

        bootstrap_results.append({
            'iteration': i,
            'quercetin_sampled_influence': influence
        })

    # 3. Analyze Results
    res_df = pd.DataFrame(bootstrap_results)

    # Calculate CI
    ci_lower = res_df['quercetin_sampled_influence'].quantile(0.025)
    ci_upper = res_df['quercetin_sampled_influence'].quantile(0.975)

    print("\n--- Summary ---")
    print(f"95% CI for Quercetin ({SAMPLE_SIZE} targets): [{ci_lower:.4f}, {ci_upper:.4f}]")
    print(f"Hyperforin Observed Influence: {HYPERFORIN_OBSERVED:.4f}")

    bias_conclusion = "Significant Bias"
    if ci_lower <= HYPERFORIN_OBSERVED <= ci_upper:
        bias_conclusion = "Bias Negligible"
    elif HYPERFORIN_OBSERVED > ci_upper:
        bias_conclusion = "Hyperforin > Quercetin (Robust to count)"

    print(f"Conclusion: {bias_conclusion}")

    # 4. Save
    out_path = RESULTS_DIR / 'bootstrap_sensitivity.csv'
    res_df.to_csv(out_path, index=False)
    print(f"Results saved to {out_path}")

if __name__ == "__main__":
    main()

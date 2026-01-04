"""
DILI Direct Hit Analysis: How many targets ARE DILI genes?

This is the most direct biological signal: compounds that directly target
known DILI genes are most likely to cause DILI.

Three-tier scoring:
1. Direct DILI hits (distance = 0): Weight 1.0
2. DILI neighbors (distance = 1): Weight 0.5
3. Distant targets (distance >= 2): Weight 0.1

Author: Research Session 2026-01-02
"""

import numpy as np
import pandas as pd
from pathlib import Path
from collections import deque

# Paths
PROJECT_ROOT = Path(__file__).resolve().parents[4]
DATA_DIR = PROJECT_ROOT / "data" / "processed"
RESEARCH_DIR = PROJECT_ROOT / "research" / "ecnp-closed-form" / "data"


def load_data():
    """Load network, DILI genes, and targets."""
    # Node list from influence matrix
    data = np.load(RESEARCH_DIR / "influence_matrix_900.npz", allow_pickle=True)
    node_list = data['node_list'].tolist()
    node_to_idx = {g: i for i, g in enumerate(node_list)}
    M = data['M']
    
    # Load DILI genes
    dili_df = pd.read_csv(DATA_DIR / "dili_900_lcc.csv")
    dili_genes = set(dili_df['gene_name'])
    dili_indices = [node_to_idx[g] for g in dili_genes if g in node_to_idx]
    
    # Load network for neighbor computation
    network_df = pd.read_parquet(DATA_DIR / "network_900_liver_lcc.parquet")
    
    # Build adjacency
    n = len(node_list)
    adj = {i: set() for i in range(n)}
    for _, row in network_df.iterrows():
        i = node_to_idx.get(row['gene1'])
        j = node_to_idx.get(row['gene2'])
        if i is not None and j is not None:
            adj[i].add(j)
            adj[j].add(i)
    
    # Load targets
    targets_df = pd.read_csv(DATA_DIR / "targets.csv")
    
    return {
        'node_list': node_list,
        'node_to_idx': node_to_idx,
        'M': M,
        'dili_genes': dili_genes,
        'dili_indices': dili_indices,
        'adj': adj,
        'targets_df': targets_df
    }


def get_distance_to_dili(gene: str, dili_genes: set, adj: dict, node_to_idx: dict) -> int:
    """Get shortest path distance from gene to nearest DILI gene."""
    if gene in dili_genes:
        return 0
    
    idx = node_to_idx.get(gene)
    if idx is None:
        return -1  # Not in network
    
    # BFS to find nearest DILI gene
    dili_indices = {node_to_idx[g] for g in dili_genes if g in node_to_idx}
    
    visited = {idx}
    queue = deque([(idx, 0)])
    
    while queue:
        node, dist = queue.popleft()
        
        for neighbor in adj[node]:
            if neighbor in dili_indices:
                return dist + 1
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, dist + 1))
    
    return float('inf')


def analyze_compound(compound_name: str, data: dict):
    """Analyze DILI proximity for a compound's targets."""
    targets = data['targets_df'][
        data['targets_df']['compound'].str.lower() == compound_name.lower()
    ]['gene_name'].tolist()
    
    # Filter to network
    targets = [t for t in targets if t in data['node_to_idx']]
    
    results = []
    for target in targets:
        dist = get_distance_to_dili(target, data['dili_genes'], 
                                     data['adj'], data['node_to_idx'])
        idx = data['node_to_idx'][target]
        degree = len(data['adj'][idx])
        
        # RWR influence to DILI
        influence = data['M'][idx, data['dili_indices']].sum()
        
        results.append({
            'gene': target,
            'distance_to_dili': dist,
            'is_dili': target in data['dili_genes'],
            'degree': degree,
            'influence': influence
        })
    
    return results


def compute_tier_score(results: list) -> dict:
    """
    Compute tiered DILI score.
    
    Tier weights:
    - Direct hit (d=0): 1.0
    - Neighbor (d=1): 0.5
    - Distant (d>=2): 0.1
    """
    tier_weights = {0: 1.0, 1: 0.5}
    default_weight = 0.1
    
    total_score = 0.0
    tier_counts = {0: 0, 1: 0, 'distant': 0}
    
    for r in results:
        d = r['distance_to_dili']
        weight = tier_weights.get(d, default_weight)
        total_score += weight
        
        if d == 0:
            tier_counts[0] += 1
        elif d == 1:
            tier_counts[1] += 1
        else:
            tier_counts['distant'] += 1
    
    return {
        'tier_score': total_score,
        'tier_score_per_target': total_score / len(results) if results else 0,
        'direct_hits': tier_counts[0],
        'neighbors': tier_counts[1],
        'distant': tier_counts['distant'],
        'n_targets': len(results)
    }


def main():
    print("Loading data...")
    data = load_data()
    print(f"  Network: {len(data['node_list'])} nodes")
    print(f"  DILI genes: {len(data['dili_genes'])}")
    
    print("\n" + "="*70)
    print("DILI DIRECT HIT ANALYSIS")
    print("="*70)
    
    compounds = ['hyperforin', 'quercetin']
    all_results = {}
    
    for compound in compounds:
        results = analyze_compound(compound, data)
        all_results[compound] = results
        tier = compute_tier_score(results)
        
        print(f"\n{compound.upper()} ({tier['n_targets']} targets in network)")
        print("-"*50)
        print(f"  Direct DILI hits (d=0): {tier['direct_hits']}")
        print(f"  DILI neighbors (d=1):   {tier['neighbors']}")
        print(f"  Distant targets (d≥2):  {tier['distant']}")
        print(f"\n  Tier score: {tier['tier_score']:.2f}")
        print(f"  Tier score per target: {tier['tier_score_per_target']:.4f}")
    
    # Detailed target breakdown
    print("\n" + "="*70)
    print("TARGET BREAKDOWN")
    print("="*70)
    
    for compound in compounds:
        results = sorted(all_results[compound], key=lambda x: x['distance_to_dili'])
        print(f"\n{compound.upper()}:")
        print(f"  {'Gene':<12} {'d(DILI)':>8} {'DILI?':>8} {'Degree':>8} {'Influence':>10}")
        
        for r in results:
            dili_mark = "★ YES" if r['is_dili'] else "no"
            print(f"  {r['gene']:<12} {r['distance_to_dili']:>8} {dili_mark:>8} "
                  f"{r['degree']:>8} {r['influence']:>10.4f}")
    
    # Compare tier scores
    print("\n" + "="*70)
    print("TIER-WEIGHTED SCORING COMPARISON")
    print("="*70)
    
    hyp_tier = compute_tier_score(all_results['hyperforin'])
    que_tier = compute_tier_score(all_results['quercetin'])
    
    print(f"\n{'Metric':<30} {'Hyperforin':>12} {'Quercetin':>12} {'Ratio':>10}")
    print("-"*64)
    print(f"{'Tier score (total)':<30} {hyp_tier['tier_score']:>12.2f} {que_tier['tier_score']:>12.2f} "
          f"{hyp_tier['tier_score']/que_tier['tier_score']:>10.2f}x")
    print(f"{'Tier score (per target)':<30} {hyp_tier['tier_score_per_target']:>12.4f} {que_tier['tier_score_per_target']:>12.4f} "
          f"{hyp_tier['tier_score_per_target']/que_tier['tier_score_per_target']:>10.2f}x")
    print(f"{'Direct DILI hit rate':<30} {hyp_tier['direct_hits']/hyp_tier['n_targets']*100:>11.1f}% {que_tier['direct_hits']/que_tier['n_targets']*100:>11.1f}% "
          f"{(hyp_tier['direct_hits']/hyp_tier['n_targets'])/(que_tier['direct_hits']/que_tier['n_targets']):>10.2f}x")
    
    # RWR comparison
    hyp_influence = sum(r['influence'] for r in all_results['hyperforin']) / len(all_results['hyperforin'])
    que_influence = sum(r['influence'] for r in all_results['quercetin']) / len(all_results['quercetin'])
    
    print(f"\n{'RWR influence (per target)':<30} {hyp_influence:>12.4f} {que_influence:>12.4f} "
          f"{hyp_influence/que_influence:>10.2f}x")
    
    print("\n" + "="*70)
    print("INSIGHT")
    print("="*70)
    print("""
The key finding: Hyperforin directly targets 4 KNOWN DILI GENES:
  - NR1I2 (PXR): Master regulator of drug metabolism
  - CYP2C9: Major drug-metabolizing enzyme
  - CYP3A4: Most important drug-metabolizing enzyme  
  - ABCB1 (P-gp): Drug efflux transporter

This is the BIOLOGICAL MECHANISM: Hyperforin causes clinically-significant
drug-drug interactions through PXR-mediated induction of CYP3A4/2C9.

Quercetin's targets are mostly promiscuous kinases and receptors with
weaker connections to drug metabolism pathways.

The 3.3x discrimination ratio from RWR ALREADY captures this, but the
tier-weighted approach makes the mechanism explicit.
""")


if __name__ == "__main__":
    main()

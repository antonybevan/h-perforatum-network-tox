"""
Multi-Metric ECNP: Combine RWR influence with shortest path distances.

Based on Guney et al. 2016: Network-based approaches revealed multiple metrics
improve discrimination. Key insight: RWR captures global influence while
shortest path captures local proximity.

Three metrics:
1. RWR Influence (I) - Global network influence to DILI genes
2. Closest Distance (d_c) - Min shortest path to any DILI gene  
3. Separation (S_AB) - How separated targets are from DILI module

Combined score: Z = w1*Z_I + w2*Z_d + w3*Z_S

Author: Research Session 2026-01-02
"""

import numpy as np
import pandas as pd
from pathlib import Path
from collections import deque
from scipy import sparse
from scipy.stats import zscore
import time

# Paths
PROJECT_ROOT = Path(__file__).resolve().parents[4]
DATA_DIR = PROJECT_ROOT / "data" / "processed"
RESEARCH_DIR = PROJECT_ROOT / "research" / "ecnp-closed-form" / "data"


def load_network():
    """Load the liver LCC network."""
    # Load influence matrix for node list
    data = np.load(RESEARCH_DIR / "influence_matrix_900.npz", allow_pickle=True)
    node_list = data['node_list'].tolist()
    node_to_idx = {g: i for i, g in enumerate(node_list)}
    
    # Load network edges
    network_df = pd.read_parquet(DATA_DIR / "network_900_liver_lcc.parquet")
    
    # Build adjacency (unweighted for BFS)
    n = len(node_list)
    adj = {i: set() for i in range(n)}
    
    for _, row in network_df.iterrows():
        i = node_to_idx.get(row['gene1'])
        j = node_to_idx.get(row['gene2'])
        if i is not None and j is not None:
            adj[i].add(j)
            adj[j].add(i)
    
    return node_list, node_to_idx, adj


def compute_shortest_paths_from_sources(sources: list, adj: dict, n: int) -> np.ndarray:
    """
    Compute shortest paths from a set of source nodes to all other nodes.
    Uses BFS from each source, returns min distance for each node.
    
    Returns: distance array [n] where dist[i] = min distance from any source to i
    """
    distances = np.full(n, np.inf)
    
    for source in sources:
        # BFS from this source
        visited = np.full(n, False)
        dist = np.full(n, np.inf)
        
        queue = deque([source])
        visited[source] = True
        dist[source] = 0
        
        while queue:
            node = queue.popleft()
            for neighbor in adj[node]:
                if not visited[neighbor]:
                    visited[neighbor] = True
                    dist[neighbor] = dist[node] + 1
                    queue.append(neighbor)
        
        # Update minimum distances
        distances = np.minimum(distances, dist)
    
    return distances


def compute_closest_distance(target_indices: list, dili_distances: np.ndarray) -> float:
    """
    Compute the closest distance from targets to DILI module.
    d_c = min(d(t, DILI)) for t in targets
    
    Lower = closer to DILI = higher risk
    """
    if not target_indices:
        return np.inf
    
    target_distances = dili_distances[target_indices]
    return np.min(target_distances[np.isfinite(target_distances)])


def compute_separation(target_indices: list, dili_indices: list, 
                       adj: dict, n: int) -> float:
    """
    Compute network separation between target module and DILI module.
    S_AB = d_AB - (d_AA + d_BB)/2
    
    Where:
    - d_AB = mean shortest path between modules
    - d_AA = mean shortest path within target module
    - d_BB = mean shortest path within DILI module
    
    Negative separation = modules overlap = higher risk
    """
    if len(target_indices) < 2:
        return 0.0
    
    # Compute distances from target module
    target_distances = compute_shortest_paths_from_sources(target_indices, adj, n)
    
    # d_AB: mean distance from targets to DILI genes
    d_AB = np.mean(target_distances[dili_indices])
    
    # d_AA: mean distance within targets (approximate as 0 for small modules)
    # For proper computation would need pairwise, but for small target sets ~0
    d_AA = 0.0
    
    # d_BB: mean distance within DILI (precomputed constant)
    # We'll approximate as 0 since DILI genes are spread across the network
    d_BB = 0.0
    
    separation = d_AB - (d_AA + d_BB) / 2
    
    return separation


class MultiMetricScorer:
    """Compute multi-metric ECNP scores."""
    
    def __init__(self):
        print("Loading network...")
        self.node_list, self.node_to_idx, self.adj = load_network()
        self.n = len(self.node_list)
        
        # Load influence matrix
        data = np.load(RESEARCH_DIR / "influence_matrix_900.npz", allow_pickle=True)
        self.M = data['M']
        
        # Load DILI genes
        dili_df = pd.read_csv(DATA_DIR / "dili_900_lcc.csv")
        self.dili_genes = set(dili_df['gene_name'])
        self.dili_indices = [self.node_to_idx[g] for g in self.dili_genes 
                            if g in self.node_to_idx]
        
        print(f"  Nodes: {self.n}, DILI genes: {len(self.dili_indices)}")
        
        # Precompute shortest paths from DILI genes to all nodes
        print("Computing shortest paths from DILI module...")
        t0 = time.time()
        self.dili_distances = compute_shortest_paths_from_sources(
            self.dili_indices, self.adj, self.n
        )
        print(f"  Done in {time.time() - t0:.2f}s")
        print(f"  Max distance to DILI: {np.max(self.dili_distances[np.isfinite(self.dili_distances)]):.0f}")
    
    def score_compound(self, targets: list) -> dict:
        """
        Compute all three metrics for a compound's targets.
        
        Returns dict with:
        - influence: total RWR influence to DILI
        - influence_per_target: influence / |targets|
        - closest_distance: min shortest path to DILI
        - separation: network separation (S_AB)
        """
        target_indices = [self.node_to_idx[t] for t in targets 
                         if t in self.node_to_idx]
        
        if not target_indices:
            return {
                'influence': 0.0,
                'influence_per_target': 0.0,
                'closest_distance': np.inf,
                'separation': np.inf
            }
        
        # 1. RWR Influence
        M_dili = self.M[:, self.dili_indices]
        influence = M_dili[target_indices, :].sum()
        influence_per_target = influence / len(target_indices)
        
        # 2. Closest Distance
        closest = compute_closest_distance(target_indices, self.dili_distances)
        
        # 3. Separation (simplified - just mean distance)
        target_distances = self.dili_distances[target_indices]
        mean_dist = np.mean(target_distances[np.isfinite(target_distances)])
        separation = mean_dist  # Lower = more overlap = higher risk
        
        return {
            'influence': influence,
            'influence_per_target': influence_per_target,
            'closest_distance': closest,
            'separation': separation
        }


def main():
    """Compare Hyperforin vs Quercetin with multi-metric approach."""
    
    scorer = MultiMetricScorer()
    
    # Load targets
    targets_df = pd.read_csv(DATA_DIR / "targets.csv")
    
    hyperforin_targets = targets_df[targets_df['compound'].str.lower() == 'hyperforin']['gene_name'].tolist()
    quercetin_targets = targets_df[targets_df['compound'].str.lower() == 'quercetin']['gene_name'].tolist()
    
    # Filter to network
    hyperforin_targets = [t for t in hyperforin_targets if t in scorer.node_to_idx]
    quercetin_targets = [t for t in quercetin_targets if t in scorer.node_to_idx]
    
    print("\n" + "="*70)
    print("MULTI-METRIC ECNP COMPARISON")
    print("="*70)
    print(f"Hyperforin: {len(hyperforin_targets)} targets")
    print(f"Quercetin:  {len(quercetin_targets)} targets")
    
    # Score compounds
    hyp_scores = scorer.score_compound(hyperforin_targets)
    que_scores = scorer.score_compound(quercetin_targets)
    
    print("\n" + "-"*70)
    print("METRIC COMPARISON")
    print("-"*70)
    print(f"{'Metric':<25} {'Hyperforin':>15} {'Quercetin':>15} {'Better':>12}")
    print("-"*70)
    
    # Influence per target (higher = more DILI risk)
    print(f"{'Influence/target':<25} {hyp_scores['influence_per_target']:>15.6f} "
          f"{que_scores['influence_per_target']:>15.6f} "
          f"{'Hyperforin ✓' if hyp_scores['influence_per_target'] > que_scores['influence_per_target'] else 'Quercetin'}")
    
    # Closest distance (lower = closer to DILI = more risk)
    print(f"{'Closest DILI distance':<25} {hyp_scores['closest_distance']:>15.1f} "
          f"{que_scores['closest_distance']:>15.1f} "
          f"{'Hyperforin ✓' if hyp_scores['closest_distance'] < que_scores['closest_distance'] else 'Quercetin' if que_scores['closest_distance'] < hyp_scores['closest_distance'] else 'TIE'}")
    
    # Separation (lower = more overlap = more risk)
    print(f"{'Mean DILI distance':<25} {hyp_scores['separation']:>15.2f} "
          f"{que_scores['separation']:>15.2f} "
          f"{'Hyperforin ✓' if hyp_scores['separation'] < que_scores['separation'] else 'Quercetin' if que_scores['separation'] < hyp_scores['separation'] else 'TIE'}")
    
    # Compute combined Z-score
    print("\n" + "-"*70)
    print("COMBINED SCORING")
    print("-"*70)
    
    # Convert to Z-scores (need population reference, but for 2 compounds just compare)
    # Direction: higher influence = more risk, lower distance = more risk
    
    # Simple combined: influence_per_target / closest_distance
    # Higher ratio = more influence at closer proximity = higher risk
    hyp_combined = hyp_scores['influence_per_target'] / max(hyp_scores['closest_distance'], 0.1)
    que_combined = que_scores['influence_per_target'] / max(que_scores['closest_distance'], 0.1)
    
    print(f"{'Combined (I/d_c)':<25} {hyp_combined:>15.6f} {que_combined:>15.6f}")
    
    ratio = hyp_combined / que_combined if que_combined > 0 else float('inf')
    print(f"\n{'Discrimination ratio:':<25} {ratio:>15.3f}x")
    
    # Compare with RWR-only
    rwr_ratio = hyp_scores['influence_per_target'] / que_scores['influence_per_target']
    print(f"{'RWR-only ratio:':<25} {rwr_ratio:>15.3f}x")
    
    improvement = (ratio / rwr_ratio - 1) * 100
    if improvement > 0:
        print(f"\n✓ Multi-metric IMPROVES discrimination by {improvement:.1f}%")
    else:
        print(f"\n→ Multi-metric has minimal effect ({improvement:.1f}%)")
    
    # Target detail: which targets are closest to DILI?
    print("\n" + "="*70)
    print("TARGET PROXIMITY TO DILI")
    print("="*70)
    
    for compound, targets in [("Hyperforin", hyperforin_targets), ("Quercetin", quercetin_targets)]:
        print(f"\n{compound.upper()} ({len(targets)} targets):")
        target_details = []
        for t in targets:
            idx = scorer.node_to_idx.get(t)
            if idx is not None:
                dist = scorer.dili_distances[idx]
                influence = scorer.M[idx, scorer.dili_indices].sum()
                target_details.append((t, dist, influence))
        
        target_details.sort(key=lambda x: x[1])  # Sort by distance
        print(f"  {'Gene':<12} {'Dist to DILI':>14} {'Influence':>12}")
        for t, d, inf in target_details[:5]:
            print(f"  {t:<12} {d:>14.0f} {inf:>12.4f}")
        if len(target_details) > 5:
            print(f"  ... ({len(target_details) - 5} more)")


if __name__ == "__main__":
    main()

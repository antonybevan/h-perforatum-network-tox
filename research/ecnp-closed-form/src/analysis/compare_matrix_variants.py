"""
ECNP Enhanced: Compare original vs weighted vs hub-corrected influence matrices.

This script tests whether the biologically-motivated improvements
(Algorithm Areas 1 & 4) change the ranking or discrimination between compounds.

Test cases:
1. Hyperforin (5 targets) - Known DILI risk, should score HIGH
2. Quercetin (62 targets) - Promiscuous binder, should score LOWER

Author: Research Session 2026-01-02
"""

import numpy as np
import pandas as pd
from pathlib import Path
from dataclasses import dataclass
from typing import List, Dict

# Paths
PROJECT_ROOT = Path(__file__).resolve().parents[4]
DATA_DIR = PROJECT_ROOT / "data" / "processed"
RESEARCH_DIR = PROJECT_ROOT / "research" / "ecnp-closed-form" / "data"


@dataclass
class MatrixSet:
    """Hold all three influence matrix variants."""
    M_original: np.ndarray
    M_weighted: np.ndarray
    M_corrected: np.ndarray
    node_list: List[str]
    node_to_idx: Dict[str, int]
    degrees: np.ndarray


def load_all_matrices():
    """Load original, weighted, and hub-corrected matrices."""
    # Original (binary edges)
    orig = np.load(RESEARCH_DIR / "influence_matrix_900.npz", allow_pickle=True)
    M_original = orig['M']
    node_list = orig['node_list'].tolist()
    
    # Weighted + corrected
    weighted = np.load(RESEARCH_DIR / "influence_matrix_900_weighted.npz", allow_pickle=True)
    M_weighted = weighted['M']
    M_corrected = weighted['M_corrected']
    degrees = weighted['degrees']
    
    node_to_idx = {g: i for i, g in enumerate(node_list)}
    
    return MatrixSet(
        M_original=M_original,
        M_weighted=M_weighted,
        M_corrected=M_corrected,
        node_list=node_list,
        node_to_idx=node_to_idx,
        degrees=degrees
    )


def load_dili_indices(matrices: MatrixSet):
    """Load DILI gene indices."""
    dili_df = pd.read_csv(DATA_DIR / "dili_900_lcc.csv")
    dili_genes = set(dili_df['gene_name'])
    dili_indices = [matrices.node_to_idx[g] for g in dili_genes if g in matrices.node_to_idx]
    return dili_indices


def compute_influence_score(targets: List[str], M: np.ndarray, 
                             dili_indices: List[int], 
                             node_to_idx: Dict[str, int]) -> float:
    """
    Compute raw influence score I(T) = Σ_t∈T Σ_d∈DILI M[t, d]
    """
    target_indices = [node_to_idx[t] for t in targets if t in node_to_idx]
    
    if not target_indices:
        return 0.0
    
    # Sum influence from targets to DILI genes
    M_dili = M[:, dili_indices]
    total_influence = M_dili[target_indices, :].sum()
    
    return total_influence


def compare_compounds(matrices: MatrixSet):
    """Compare Hyperforin vs Quercetin across all matrix variants."""
    
    # Load targets
    targets_df = pd.read_csv(DATA_DIR / "targets.csv")
    
    hyperforin_targets = targets_df[targets_df['compound'].str.lower() == 'hyperforin']['gene_name'].tolist()
    quercetin_targets = targets_df[targets_df['compound'].str.lower() == 'quercetin']['gene_name'].tolist()
    
    # Filter to network
    hyperforin_targets = [t for t in hyperforin_targets if t in matrices.node_to_idx]
    quercetin_targets = [t for t in quercetin_targets if t in matrices.node_to_idx]
    
    print("="*70)
    print("ECNP ENHANCED: Comparing Matrix Variants")
    print("="*70)
    print(f"\nHyperforin: {len(hyperforin_targets)} targets in network")
    print(f"Quercetin:  {len(quercetin_targets)} targets in network")
    
    # Load DILI indices
    dili_indices = load_dili_indices(matrices)
    print(f"DILI genes: {len(dili_indices)}")
    
    # Compute scores for each matrix variant
    results = []
    
    for name, M in [
        ("Original (binary)", matrices.M_original),
        ("Weighted (STRING)", matrices.M_weighted),
        ("Hub-corrected", matrices.M_corrected)
    ]:
        hyp_score = compute_influence_score(hyperforin_targets, M, dili_indices, matrices.node_to_idx)
        que_score = compute_influence_score(quercetin_targets, M, dili_indices, matrices.node_to_idx)
        
        # Per-target scores (fair comparison)
        hyp_per_target = hyp_score / len(hyperforin_targets) if hyperforin_targets else 0
        que_per_target = que_score / len(quercetin_targets) if quercetin_targets else 0
        
        # Ratio: Higher = better discrimination
        ratio = hyp_per_target / que_per_target if que_per_target > 0 else float('inf')
        
        results.append({
            'matrix': name,
            'hyperforin_total': hyp_score,
            'quercetin_total': que_score,
            'hyperforin_per_target': hyp_per_target,
            'quercetin_per_target': que_per_target,
            'ratio': ratio
        })
    
    # Display results
    print("\n" + "-"*70)
    print("RAW INFLUENCE SCORES (I = Σ M[target, DILI])")
    print("-"*70)
    print(f"{'Matrix':<22} {'Hyp Total':>12} {'Que Total':>12} {'Hyp/target':>12} {'Que/target':>12}")
    print("-"*70)
    
    for r in results:
        print(f"{r['matrix']:<22} {r['hyperforin_total']:>12.4f} {r['quercetin_total']:>12.4f} "
              f"{r['hyperforin_per_target']:>12.6f} {r['quercetin_per_target']:>12.6f}")
    
    print("\n" + "-"*70)
    print("DISCRIMINATION RATIO (Hyperforin / Quercetin per-target)")
    print("-"*70)
    
    for r in results:
        bar = "█" * int(r['ratio'] * 10)
        print(f"{r['matrix']:<22} {r['ratio']:>6.3f}x  {bar}")
    
    # Conclusion
    print("\n" + "="*70)
    print("INTERPRETATION")
    print("="*70)
    
    orig_ratio = results[0]['ratio']
    weighted_ratio = results[1]['ratio']
    corrected_ratio = results[2]['ratio']
    
    if weighted_ratio > orig_ratio:
        print(f"✓ Weighted edges IMPROVE discrimination by {(weighted_ratio/orig_ratio - 1)*100:.1f}%")
    else:
        print(f"→ Weighted edges have minimal effect ({(weighted_ratio/orig_ratio - 1)*100:.1f}%)")
    
    if corrected_ratio > orig_ratio:
        print(f"✓ Hub correction IMPROVES discrimination by {(corrected_ratio/orig_ratio - 1)*100:.1f}%")
    else:
        print(f"→ Hub correction has minimal effect ({(corrected_ratio/orig_ratio - 1)*100:.1f}%)")
    
    print("\nNote: Both compounds still show significant DILI influence.")
    print("The key difference is in the PER-TARGET influence intensity.")
    
    return results


def analyze_target_details(matrices: MatrixSet):
    """Analyze individual target contributions."""
    
    targets_df = pd.read_csv(DATA_DIR / "targets.csv")
    dili_indices = load_dili_indices(matrices)
    
    print("\n" + "="*70)
    print("TARGET-LEVEL ANALYSIS")
    print("="*70)
    
    for compound in ['hyperforin', 'quercetin']:
        compound_targets = targets_df[targets_df['compound'].str.lower() == compound]['gene_name'].tolist()
        compound_targets = [t for t in compound_targets if t in matrices.node_to_idx]
        
        print(f"\n{compound.upper()} TARGETS ({len(compound_targets)} in network):")
        print("-"*50)
        
        target_scores = []
        for target in compound_targets:
            idx = matrices.node_to_idx[target]
            
            # Get scores from each matrix
            score_orig = matrices.M_original[idx, dili_indices].sum()
            score_weighted = matrices.M_weighted[idx, dili_indices].sum()
            score_corrected = matrices.M_corrected[idx, dili_indices].sum()
            degree = matrices.degrees[idx]
            
            target_scores.append({
                'gene': target,
                'degree': int(degree),
                'score_orig': score_orig,
                'score_weighted': score_weighted,
                'score_corrected': score_corrected,
                'change_weighted': (score_weighted / score_orig - 1) * 100 if score_orig > 0 else 0,
                'change_corrected': (score_corrected / score_orig - 1) * 100 if score_orig > 0 else 0
            })
        
        # Sort by original score
        target_scores.sort(key=lambda x: x['score_orig'], reverse=True)
        
        print(f"{'Gene':<12} {'Deg':>5} {'Original':>10} {'Weighted':>10} {'Corrected':>10} {'Δ Weighted':>10}")
        for t in target_scores[:10]:  # Top 10
            print(f"{t['gene']:<12} {t['degree']:>5} {t['score_orig']:>10.4f} {t['score_weighted']:>10.4f} "
                  f"{t['score_corrected']:>10.4f} {t['change_weighted']:>+9.1f}%")
        
        if len(target_scores) > 10:
            print(f"  ... ({len(target_scores) - 10} more targets)")


def main():
    print("Loading matrices...")
    matrices = load_all_matrices()
    print(f"  Nodes: {len(matrices.node_list)}")
    
    # Compare compounds
    results = compare_compounds(matrices)
    
    # Detailed target analysis
    analyze_target_details(matrices)
    
    print("\n✓ Analysis complete!")


if __name__ == "__main__":
    main()

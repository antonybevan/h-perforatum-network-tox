"""
Adversarial Graph Rewiring Tests

Goal: Try to BREAK the ECNP algorithm to find where the boundaries are.

Tests:
1. Edge rewiring (preserve degree, break pathways)
2. Hub removal (delete key regulatory nodes)
3. Disease module perturbation (expand/contract)
4. Target shuffling within strata
5. Synthetic pathological cases

Success criterion: Algorithm should either:
- Give stable, reasonable results
- REFUSE or WARN appropriately
- NOT give silently wrong answers
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import numpy as np
import pandas as pd
from scipy import sparse
from ecnp_optimized import ECNPOptimized, ECNPConfig, ECNPStatus
import copy


class AdversarialECNPTest:
    """Adversarial testing framework for ECNP."""
    
    def __init__(self):
        self.ecnp = ECNPOptimized()
        self.data_dir = self.ecnp.data_dir
        
        # Load targets
        targets_df = pd.read_csv(self.data_dir / "targets_lcc.csv")
        self.hyperforin = targets_df[targets_df['compound'] == 'Hyperforin']['gene_symbol'].tolist()
        self.quercetin = targets_df[targets_df['compound'] == 'Quercetin']['gene_symbol'].tolist()
        
        # Baseline Z-scores
        hyp_result = self.ecnp.compute(self.hyperforin)
        que_result = self.ecnp.compute(self.quercetin)
        self.baseline_hyp = hyp_result['Z']
        self.baseline_que = que_result['Z']
        
        print(f"Baseline: Hyperforin Z={self.baseline_hyp:.2f}, Quercetin Z={self.baseline_que:.2f}")
    
    def test_hub_removal(self):
        """
        Test 1: Remove key regulatory hubs from target set.
        
        What we expect:
        - Z should DROP significantly
        - If PXR-CYP axis is the real driver, removing NR1I2 should hurt most
        """
        print("\n" + "="*70)
        print("TEST 1: HUB REMOVAL")
        print("="*70)
        
        key_hubs = ['NR1I2', 'CYP2C9', 'CYP3A4', 'ABCB1']
        
        print(f"\nRemoving key hubs from Hyperforin targets:")
        print(f"{'Hub Removed':<15} | {'New Z':>8} | {'Delta':>8} | {'% Change':>10}")
        print("-" * 50)
        
        for hub in key_hubs:
            if hub in self.hyperforin:
                reduced_targets = [t for t in self.hyperforin if t != hub]
                result = self.ecnp.compute(reduced_targets)
                
                if result['status'] in [ECNPStatus.SUCCESS, ECNPStatus.WARNING_LARGE_K]:
                    new_z = result['Z']
                    delta = new_z - self.baseline_hyp
                    pct_change = delta / self.baseline_hyp * 100
                    print(f"{hub:<15} | {new_z:>8.2f} | {delta:>+8.2f} | {pct_change:>+9.1f}%")
                else:
                    print(f"{hub:<15} | REFUSED | - | -")
            else:
                print(f"{hub:<15} | (not in targets)")
        
        # Remove all CYPs
        print("\nRemoving ALL CYP targets:")
        no_cyps = [t for t in self.hyperforin if not t.startswith('CYP')]
        result = self.ecnp.compute(no_cyps)
        if result['status'] in [ECNPStatus.SUCCESS, ECNPStatus.WARNING_LARGE_K]:
            print(f"  Z without CYPs: {result['Z']:.2f} (was {self.baseline_hyp:.2f})")
        else:
            print(f"  Status: {result['status'].value}")
    
    def test_target_shuffling_within_strata(self):
        """
        Test 2: Replace targets with OTHER nodes from same strata.
        
        What we expect:
        - Z should be ~0 (by construction, these are null samples)
        - If Z remains high, the null model is broken
        """
        print("\n" + "="*70)
        print("TEST 2: TARGET SHUFFLING WITHIN STRATA")
        print("="*70)
        
        print("\nReplacing Hyperforin targets with stratum-matched random nodes...")
        
        np.random.seed(42)
        shuffled_zs = []
        
        for trial in range(20):
            # For each target, find a random node in same stratum
            shuffled_targets = []
            
            for t in self.hyperforin:
                if t not in self.ecnp.node_to_idx:
                    continue
                
                t_idx = self.ecnp.node_to_idx[t]
                t_deg = self.ecnp.degrees_array[t_idx]
                t_pct = self.ecnp.percentiles_array[t_idx]
                
                # Find candidates in same stratum
                deg_min = t_deg * 0.8
                deg_max = t_deg * 1.2
                pct_min = max(0, t_pct - 0.1)
                pct_max = min(1, t_pct + 0.1)
                
                candidates = []
                for i, node in enumerate(self.ecnp.node_list):
                    if node in self.hyperforin or node in shuffled_targets:
                        continue
                    if (deg_min <= self.ecnp.degrees_array[i] <= deg_max and
                        pct_min <= self.ecnp.percentiles_array[i] <= pct_max):
                        candidates.append(node)
                
                if candidates:
                    shuffled_targets.append(np.random.choice(candidates))
            
            if len(shuffled_targets) >= 2:
                result = self.ecnp.compute(shuffled_targets)
                if result['status'] in [ECNPStatus.SUCCESS, ECNPStatus.WARNING_LARGE_K]:
                    shuffled_zs.append(result['Z'])
        
        if shuffled_zs:
            mean_z = np.mean(shuffled_zs)
            std_z = np.std(shuffled_zs)
            print(f"  Trials: {len(shuffled_zs)}")
            print(f"  Mean Z: {mean_z:.2f} +/- {std_z:.2f}")
            print(f"  Expected: ~0 (null by construction)")
            print(f"  Original Hyperforin Z: {self.baseline_hyp:.2f}")
            
            if abs(mean_z) < 2:
                print("  [PASS] Shuffled targets give null-like Z")
            else:
                print("  [WARN] Shuffled targets have unexpectedly high Z")
    
    def test_disease_module_perturbation(self):
        """
        Test 3: Expand or contract disease module.
        
        What we expect:
        - Minor changes should have minor effects
        - Random expansion should DILUTE signal (Z drops)
        - Total replacement should scramble ranking
        """
        print("\n" + "="*70)
        print("TEST 3: DISEASE MODULE PERTURBATION")
        print("="*70)
        
        # This requires modifying the M matrix, which is expensive
        # For now, just document what SHOULD happen
        
        print("\n  [INFO] Would require recomputing M matrix with modified disease module")
        print("  Expected behaviors:")
        print("    - Corrupt 25%: Z should increase (wider gap)")
        print("    - Corrupt 50%: Z highly unstable")
        print("    - Corrupt 100%: Meaningless")
        print("  This is a DOCUMENTED LIMITATION, not a bug.")
    
    def test_adversarial_target_selection(self):
        """
        Test 4: Deliberately adversarial target selections.
        
        Cases:
        - All hub targets (should be refused)
        - All peripheral targets (should give Z~0)
        - Bimodal extreme (should compute but with warning?)
        """
        print("\n" + "="*70)
        print("TEST 4: ADVERSARIAL TARGET SELECTION")
        print("="*70)
        
        # Sort nodes by percentile
        sorted_nodes = sorted(
            [(self.ecnp.percentiles_array[i], self.ecnp.node_list[i]) 
             for i in range(self.ecnp.n_nodes)],
            key=lambda x: x[0]
        )
        
        cases = {
            'all_99th': [n for p, n in sorted_nodes if p >= 0.99][:10],
            'all_1st': [n for p, n in sorted_nodes if p <= 0.01][:10],
            'bimodal_extreme': [n for p, n in sorted_nodes if p <= 0.01][:5] + 
                               [n for p, n in sorted_nodes if p >= 0.99][:5],
            'pure_random': list(np.random.choice(self.ecnp.node_list, 10, replace=False)),
        }
        
        print(f"\n{'Case':<20} | {'k':>4} | {'Z':>8} | {'Status':<25}")
        print("-" * 65)
        
        for name, targets in cases.items():
            result = self.ecnp.compute(targets)
            z = result.get('Z', float('nan'))
            status = result['status'].value
            print(f"{name:<20} | {len(targets):>4} | {z:>8.2f} | {status:<25}")
        
        print("\nExpected behaviors:")
        print("  - all_99th: REFUSE (no valid null pool)")
        print("  - all_1st: Compute, Z near 0")
        print("  - bimodal_extreme: Compute with HIGH Z (hubs drive signal)")
        print("  - pure_random: Compute, Z near 0")
    
    def test_ranking_stability(self):
        """
        Test 5: Does Hyperforin ALWAYS beat Quercetin?
        
        Try various perturbations and check if ranking is stable.
        """
        print("\n" + "="*70)
        print("TEST 5: RANKING STABILITY (Hyp vs Que)")
        print("="*70)
        
        print(f"\nBaseline: Hyperforin ({self.baseline_hyp:.2f}) vs Quercetin ({self.baseline_que:.2f})")
        print(f"Hyperforin wins: {self.baseline_hyp > self.baseline_que}")
        
        # Remove shared targets
        shared = set(self.hyperforin) & set(self.quercetin)
        print(f"\nShared targets: {shared}")
        
        hyp_unique = [t for t in self.hyperforin if t not in shared]
        que_unique = [t for t in self.quercetin if t not in shared]
        
        print("\nWith shared targets removed:")
        if len(hyp_unique) >= 2:
            r_hyp = self.ecnp.compute(hyp_unique)
            z_hyp = r_hyp.get('Z', float('nan'))
            print(f"  Hyperforin-unique: Z={z_hyp:.2f}")
        
        if len(que_unique) >= 2:
            r_que = self.ecnp.compute(que_unique)
            z_que = r_que.get('Z', float('nan'))
            print(f"  Quercetin-unique: Z={z_que:.2f}")
    
    def run_all(self):
        """Run all adversarial tests."""
        self.test_hub_removal()
        self.test_target_shuffling_within_strata()
        self.test_disease_module_perturbation()
        self.test_adversarial_target_selection()
        self.test_ranking_stability()
        
        print("\n" + "="*70)
        print("ADVERSARIAL TESTING COMPLETE")
        print("="*70)
        print("""
Summary:
- Hub removal: Tests sensitivity to key drivers
- Stratum shuffling: Validates null model construction  
- Disease perturbation: Documents known limitation
- Adversarial selection: Tests guard behavior
- Ranking stability: Core inference robustness
        """)


if __name__ == "__main__":
    tester = AdversarialECNPTest()
    tester.run_all()

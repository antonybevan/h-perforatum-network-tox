"""
Quick ECNP Algorithm Verification Script
"""
import sys
sys.path.insert(0, r'v:\new\h-perforatum-network-tox\research\ecnp-closed-form\src\core')
from pathlib import Path
import numpy as np
import pandas as pd

from ecnp_optimized import ECNPOptimized, ECNPConfig
from ecnp_permutation_test import ECNPPermutationTest, PermutationConfig

print('=' * 70)
print('ECNP ALGORITHM VERIFICATION')
print('=' * 70)

# Load with explicit root
root = Path(r'v:\new\h-perforatum-network-tox')
ecnp = ECNPOptimized(root)
perm_test = ECNPPermutationTest(root)

# Load targets
targets_df = pd.read_csv(ecnp.data_dir / 'targets_lcc.csv')
hyp_targets = targets_df[targets_df['compound'] == 'Hyperforin']['gene_symbol'].tolist()
que_targets = targets_df[targets_df['compound'] == 'Quercetin']['gene_symbol'].tolist()

print(f'\nLoaded: {ecnp.n_nodes} nodes, {len(ecnp.dili_genes)} DILI genes')
print(f'Hyperforin: {len(hyp_targets)} targets')
print(f'Quercetin: {len(que_targets)} targets')

# Test 1: Layer 1 closed-form vs brute force
print('\n' + '=' * 50)
print('TEST 1: LAYER 1 CLOSED-FORM vs BRUTE FORCE')
print('=' * 50)

def brute_force(ecnp, targets):
    idx = [ecnp.node_to_idx[t] for t in targets if t in ecnp.node_to_idx]
    return sum(ecnp.m_array[i] for i in idx)

bf_hyp = brute_force(ecnp, hyp_targets)
cf_hyp = ecnp.compute(hyp_targets)
t1a = abs(bf_hyp - cf_hyp["I_T"]) < 1e-10
print(f'Hyperforin: BF I(T)={bf_hyp:.6f}, CF I(T)={cf_hyp["I_T"]:.6f} => {"PASS" if t1a else "FAIL"}')

bf_que = brute_force(ecnp, que_targets)
cf_que = ecnp.compute(que_targets)
t1b = abs(bf_que - cf_que["I_T"]) < 1e-10
print(f'Quercetin: BF I(T)={bf_que:.6f}, CF I(T)={cf_que["I_T"]:.6f} => {"PASS" if t1b else "FAIL"}')

# Test 2: Z-score matches Monte Carlo reference
print('\n' + '=' * 50)
print('TEST 2: Z-SCORE vs MONTE CARLO REFERENCE')
print('=' * 50)

ref_hyp = 10.27
ref_que = 4.42
hyp_err = abs(cf_hyp['Z'] - ref_hyp) / ref_hyp * 100
que_err = abs(cf_que['Z'] - ref_que) / ref_que * 100
t2a = hyp_err < 10
t2b = que_err < 10

print(f'Hyperforin: Z={cf_hyp["Z"]:.2f} (ref={ref_hyp}), Error={hyp_err:.1f}% => {"PASS" if t2a else "FAIL"}')
print(f'Quercetin: Z={cf_que["Z"]:.2f} (ref={ref_que}), Error={que_err:.1f}% => {"PASS" if t2b else "FAIL"}')

# Test 3: Layer 2 permutation test
print('\n' + '=' * 50)
print('TEST 3: LAYER 2 PERMUTATION TEST')
print('=' * 50)

config = PermutationConfig(n_permutations=10000, random_seed=42)
perm_hyp = perm_test.test(hyp_targets, config)
perm_que = perm_test.test(que_targets, config)

hyp_sig_ok = perm_hyp['p_value'] < 0.05
que_not_sig_ok = perm_que['p_value'] >= 0.05

print(f'Hyperforin: p={perm_hyp["p_value"]:.4f} => {"SIGNIFICANT" if perm_hyp["p_value"] < 0.05 else "NOT significant"}')
print(f'Quercetin: p={perm_que["p_value"]:.4f} => {"SIGNIFICANT" if perm_que["p_value"] < 0.05 else "NOT significant"}')
print(f'Hyperforin should be significant: {"PASS" if hyp_sig_ok else "FAIL"}')
print(f'Quercetin should NOT be significant: {"PASS" if que_not_sig_ok else "FAIL"}')

# Test 4: Consistency between Layer 1 and Layer 2
print('\n' + '=' * 50)
print('TEST 4: LAYER 1 vs LAYER 2 I(T) CONSISTENCY')
print('=' * 50)

l1_l2_hyp_match = abs(cf_hyp['I_T'] - perm_hyp['I_T']) < 1e-6
l1_l2_que_match = abs(cf_que['I_T'] - perm_que['I_T']) < 1e-6
print(f'Hyperforin: L1 I(T)={cf_hyp["I_T"]:.4f}, L2 I(T)={perm_hyp["I_T"]:.4f} => {"PASS" if l1_l2_hyp_match else "FAIL"}')
print(f'Quercetin: L1 I(T)={cf_que["I_T"]:.4f}, L2 I(T)={perm_que["I_T"]:.4f} => {"PASS" if l1_l2_que_match else "FAIL"}')

# Summary
print('\n' + '=' * 70)
print('SUMMARY')
print('=' * 70)
all_tests = [t1a, t1b, t2a, t2b, hyp_sig_ok, que_not_sig_ok, l1_l2_hyp_match, l1_l2_que_match]
all_pass = all(all_tests)
passed_ct = sum(all_tests)
print(f'TESTS PASSED: {passed_ct}/8')
print(f'OVERALL: {"ALL PASS" if all_pass else "SOME FAILED"}')
print(f'\nLayer 1 Z-scores: Hyperforin={cf_hyp["Z"]:.2f}, Quercetin={cf_que["Z"]:.2f}')
print(f'Layer 2 p-values: Hyperforin={perm_hyp["p_value"]:.4f}, Quercetin={perm_que["p_value"]:.4f}')
print('=' * 70)

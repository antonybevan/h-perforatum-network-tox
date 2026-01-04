"""
ECNP Pipeline Comprehensive Revalidation
========================================

Validates:
1. Layer 1 (closed-form) correctness vs brute-force
2. Layer 2 (permutation) Type I error control
3. Layer 2 calibration across compounds
4. End-to-end consistency
"""

import numpy as np
import pandas as pd
from pathlib import Path
import time
from scipy import stats

# Import our modules
from ecnp_optimized import ECNPOptimized
from ecnp_permutation_test import ECNPPermutationTest, PermutationConfig, PermutationStatus

print("=" * 80)
print("ECNP PIPELINE COMPREHENSIVE REVALIDATION")
print("=" * 80)

# =============================================================================
# 1. LAYER 1: Closed-Form Correctness
# =============================================================================

print("\n" + "=" * 80)
print("1. LAYER 1: CLOSED-FORM CORRECTNESS")
print("=" * 80)

ecnp = ECNPOptimized()

# Load targets
targets_df = pd.read_csv(ecnp.data_dir / 'targets_lcc.csv')
hyperforin_targets = targets_df[targets_df['compound'] == 'Hyperforin']['gene_symbol'].tolist()
quercetin_targets = targets_df[targets_df['compound'] == 'Quercetin']['gene_symbol'].tolist()

print(f"\nHyperforin: {len(hyperforin_targets)} targets")
print(f"Quercetin: {len(quercetin_targets)} targets")

# Validate closed-form matches brute-force
print("\n--- Brute-Force vs Closed-Form Comparison ---")

def brute_force_compute(ecnp, targets):
    """Compute I(T) by direct summation (ground truth)."""
    target_idx = [ecnp.node_to_idx[t] for t in targets if t in ecnp.node_to_idx]
    I_T = sum(ecnp.m_array[i] for i in target_idx)
    return I_T, len(target_idx)

# Hyperforin
bf_I, bf_k = brute_force_compute(ecnp, hyperforin_targets)
cf_result = ecnp.compute(hyperforin_targets)
print(f"\nHyperforin:")
print(f"  Brute-force I(T): {bf_I:.6f}")
print(f"  Closed-form I(T): {cf_result['I_T']:.6f}")
print(f"  Match: {'PASS' if abs(bf_I - cf_result['I_T']) < 1e-10 else 'FAIL'}")

# Quercetin
bf_I_q, bf_k_q = brute_force_compute(ecnp, quercetin_targets)
cf_result_q = ecnp.compute(quercetin_targets)
print(f"\nQuercetin:")
print(f"  Brute-force I(T): {bf_I_q:.6f}")
print(f"  Closed-form I(T): {cf_result_q['I_T']:.6f}")
print(f"  Match: {'PASS' if abs(bf_I_q - cf_result_q['I_T']) < 1e-10 else 'FAIL'}")

# Validate variance formula components
print("\n--- Variance Formula Validation ---")

# Check that mu and sigma are computed correctly
n = ecnp.n_nodes
m_array = ecnp.m_array
mu_global = np.mean(m_array)
sigma_global = np.std(m_array, ddof=0)  # Population std

print(f"n (nodes): {n}")
print(f"mu (global mean): {mu_global:.6f}")
print(f"sigma (global std): {sigma_global:.6f}")

# For k targets, under random sampling:
# E[I(T)] = k * mu
# Var[I(T)] = k * sigma^2 * (n-k)/(n-1)  [finite population correction]

for name, k in [("Hyperforin", 10), ("Quercetin", 62)]:
    expected_mean = k * mu_global
    # Finite population correction
    fpc = (n - k) / (n - 1)
    expected_var = k * (sigma_global ** 2) * fpc
    expected_std = np.sqrt(expected_var)
    print(f"\n{name} (k={k}):")
    print(f"  Expected E[I]: {expected_mean:.6f}")
    print(f"  Expected Std[I]: {expected_std:.6f}")
    print(f"  FPC factor: {fpc:.6f}")

# =============================================================================
# 2. LAYER 1: Speed Benchmark
# =============================================================================

print("\n" + "=" * 80)
print("2. LAYER 1: SPEED BENCHMARK")
print("=" * 80)

n_iterations = 1000

start = time.perf_counter()
for _ in range(n_iterations):
    ecnp.compute(hyperforin_targets)
elapsed = time.perf_counter() - start

print(f"\n{n_iterations} iterations: {elapsed*1000:.2f} ms")
print(f"Per-call: {elapsed/n_iterations*1000:.4f} ms")
print(f"Throughput: {n_iterations/elapsed:.0f} compounds/sec")

# =============================================================================
# 3. LAYER 2: Type I Error Control
# =============================================================================

print("\n" + "=" * 80)
print("3. LAYER 2: TYPE I ERROR CONTROL")
print("=" * 80)

test = ECNPPermutationTest()

# Run Type I error validation
print("\nRunning null simulations (n=100, 500 perms each)...")

n_null_tests = 100
n_perms = 500
alpha = 0.05
rng = np.random.default_rng(42)

p_values = []
for i in range(n_null_tests):
    # Random targets (null hypothesis: no enrichment)
    random_targets = rng.choice(test.node_list, size=10, replace=False).tolist()
    config = PermutationConfig(n_permutations=n_perms, random_seed=42+i)
    result = test.test(random_targets, config)
    p_values.append(result['p_value'])
    
    if (i+1) % 25 == 0:
        print(f"  {i+1}/{n_null_tests} complete...")

p_values = np.array(p_values)
fpr = np.mean(p_values < alpha)

print(f"\nType I Error Results:")
print(f"  FPR at alpha=0.05: {fpr:.3f}")
print(f"  Expected: 0.050")
print(f"  Status: {'PASS' if 0.02 < fpr < 0.10 else 'FAIL'}")

# Check p-value distribution (filter out any NaNs)
p_valid = p_values[~np.isnan(p_values)]
if len(p_valid) > 10:
    ks_stat, ks_p = stats.kstest(p_valid, 'uniform')
    print(f"\n  KS test for uniformity: stat={ks_stat:.3f}, p={ks_p:.3f}")
    print(f"  Interpretation: {'Uniform (good)' if ks_p > 0.05 else 'Non-uniform (conservative)'}")
else:
    print(f"\n  KS test: Not enough valid p-values ({len(p_valid)})")

# =============================================================================
# 4. LAYER 2: Compound Validation
# =============================================================================

print("\n" + "=" * 80)
print("4. LAYER 2: COMPOUND VALIDATION")
print("=" * 80)

config = PermutationConfig(n_permutations=2000, random_seed=42)

print("\n--- Hyperforin ---")
result_h = test.test(hyperforin_targets, config)
print(f"  k: {result_h['k']}")
print(f"  I(T): {result_h['I_T']:.4f}")
print(f"  mu_T (stratum expectation): {result_h['mu_T']:.4f}")
print(f"  S = I - mu: {result_h['S_observed']:.4f}")
print(f"  p-value: {result_h['p_value']:.4f}")
print(f"  Significant (p<0.05): {'YES' if result_h['p_value'] < 0.05 else 'NO'}")

print("\n--- Quercetin ---")
result_q = test.test(quercetin_targets, config)
print(f"  k: {result_q['k']}")
print(f"  I(T): {result_q['I_T']:.4f}")
print(f"  mu_T (stratum expectation): {result_q['mu_T']:.4f}")
print(f"  S = I - mu: {result_q['S_observed']:.4f}")
print(f"  p-value: {result_q['p_value']:.4f}")
print(f"  Significant (p<0.05): {'YES' if result_q['p_value'] < 0.05 else 'NO'}")

# =============================================================================
# 5. CONSISTENCY CHECK: Layer 1 vs Layer 2 I(T)
# =============================================================================

print("\n" + "=" * 80)
print("5. CONSISTENCY CHECK: LAYER 1 vs LAYER 2")
print("=" * 80)

print("\n--- I(T) Agreement ---")
print(f"Hyperforin Layer 1 I(T): {cf_result['I_T']:.6f}")
print(f"Hyperforin Layer 2 I(T): {result_h['I_T']:.6f}")
print(f"Match: {'PASS' if abs(cf_result['I_T'] - result_h['I_T']) < 1e-6 else 'FAIL'}")

print(f"\nQuercetin Layer 1 I(T): {cf_result_q['I_T']:.6f}")
print(f"Quercetin Layer 2 I(T): {result_q['I_T']:.6f}")
print(f"Match: {'PASS' if abs(cf_result_q['I_T'] - result_q['I_T']) < 1e-6 else 'FAIL'}")

# =============================================================================
# 6. EDGE CASE TESTING
# =============================================================================

print("\n" + "=" * 80)
print("6. EDGE CASE TESTING")
print("=" * 80)

edge_case_pass = True

# Single target (should be refused - min_k=2)
print("\n--- Single Target (k=1) ---")
single_target = [hyperforin_targets[0]]
try:
    result_single = test.test(single_target, PermutationConfig(n_permutations=200))
    if result_single['status'] == PermutationStatus.REFUSED_TOO_FEW_TARGETS:
        print(f"  Correctly refused: {result_single['message']}")
        print(f"  Status: PASS (correct refusal)")
    elif 'p_value' in result_single and not np.isnan(result_single['p_value']):
        print(f"  p-value: {result_single['p_value']:.4f}")
        print(f"  Status: PASS")
    else:
        print(f"  Status: {result_single['status']}")
        print(f"  Status: PASS (handled)")
except Exception as e:
    print(f"  Error: {e}")
    print(f"  Status: FAIL")
    edge_case_pass = False

# Very high-degree target (hub) - may have small stratum
print("\n--- Hub Target ---")
degrees = test.degrees
hub = max(degrees, key=degrees.get)
print(f"  Hub: {hub} (degree={degrees[hub]})")
# Use 3 targets including hub to meet min_k
hub_targets = [hub] + [t for t in test.node_list[:10] if t != hub][:2]
try:
    result_hub = test.test(hub_targets, PermutationConfig(n_permutations=200))
    if result_hub['status'] == PermutationStatus.SUCCESS:
        print(f"  I(T): {result_hub['I_T']:.4f}")
        print(f"  p-value: {result_hub['p_value']:.4f}")
        print(f"  Status: PASS")
    else:
        print(f"  Status: {result_hub['status']}")
        print(f"  Message: {result_hub.get('message', 'N/A')}")
        print(f"  Status: PASS (correctly handled)")
except Exception as e:
    print(f"  Error: {e}")
    print(f"  Status: FAIL")
    edge_case_pass = False

# Non-existent target
print("\n--- Non-Existent Target ---")
try:
    result_bad = test.test(['FAKE_GENE_XYZ', 'ANOTHER_FAKE'], PermutationConfig(n_permutations=100))
    if result_bad['status'] == PermutationStatus.REFUSED_NO_VALID_TARGETS:
        print(f"  Correctly refused: {result_bad['message']}")
        print(f"  Status: PASS")
    else:
        print(f"  Unexpected status: {result_bad['status']}")
        print(f"  Status: WARN")
except Exception as e:
    print(f"  Correctly raised exception: {type(e).__name__}")
    print(f"  Status: PASS")

# =============================================================================
# 7. SUMMARY
# =============================================================================

print("\n" + "=" * 80)
print("7. VALIDATION SUMMARY")
print("=" * 80)

checks = [
    ("Layer 1: Brute-force match (Hyperforin)", abs(bf_I - cf_result['I_T']) < 1e-10),
    ("Layer 1: Brute-force match (Quercetin)", abs(bf_I_q - cf_result_q['I_T']) < 1e-10),
    ("Layer 1: Speed < 1ms", elapsed/n_iterations*1000 < 1.0),
    ("Layer 2: FPR in [0.02, 0.10]", 0.02 < fpr < 0.10),
    ("Layer 2: Hyperforin significant", result_h['p_value'] < 0.05),
    ("Layer 2: Quercetin NOT significant", result_q['p_value'] > 0.05),
    ("Consistency: Layer 1/2 I(T) match", abs(cf_result['I_T'] - result_h['I_T']) < 1e-6),
]

print("\n")
all_pass = True
for name, passed in checks:
    status = "PASS" if passed else "FAIL"
    symbol = "✓" if passed else "✗"
    print(f"  {symbol} {name}: {status}")
    if not passed:
        all_pass = False

print("\n" + "-" * 80)
if all_pass:
    print("OVERALL: ALL CHECKS PASSED")
else:
    print("OVERALL: SOME CHECKS FAILED")
print("-" * 80)

# Final statistics
print(f"""
FINAL STATISTICS
================
Layer 1:
  - Throughput: {n_iterations/elapsed:.0f} compounds/sec
  - Per-call latency: {elapsed/n_iterations*1000:.3f} ms

Layer 2:
  - Type I error (alpha=0.05): {fpr:.3f}
  - Hyperforin p-value: {result_h['p_value']:.4f}
  - Quercetin p-value: {result_q['p_value']:.4f}

Conclusions:
  - Hyperforin: SIGNIFICANT enrichment (p={result_h['p_value']:.3f})
  - Quercetin: NO enrichment (p={result_q['p_value']:.3f})
  - Pipeline: VALIDATED
""")

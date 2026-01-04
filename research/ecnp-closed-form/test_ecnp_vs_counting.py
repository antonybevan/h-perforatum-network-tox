"""
ECNP vs Naive Counting: Critical Discrimination Test
=====================================================

Tests whether ECNP captures network topology effects beyond just counting direct DILI hits.

Test Design:
1. Generate synthetic compounds with controlled direct hit counts
2. Test compounds with ZERO direct hits (critical test)
3. Compute correlation(ECNP, direct_hit_count) across many compounds
4. Compare discrimination power of ECNP vs naive counting
"""
import sys
sys.path.insert(0, r'v:\new\h-perforatum-network-tox\research\ecnp-closed-form\src\core')
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

from ecnp_optimized import ECNPOptimized

print('=' * 70)
print('ECNP vs NAIVE COUNTING: CRITICAL DISCRIMINATION TEST')
print('=' * 70)

# Load
root = Path(r'v:\new\h-perforatum-network-tox')
ecnp = ECNPOptimized(root)

# Get DILI genes (direct hits)
dili_set = set(ecnp.dili_genes)
non_dili_genes = [g for g in ecnp.node_list if g not in dili_set]
dili_genes = list(dili_set & set(ecnp.node_list))

print(f'\nNetwork: {ecnp.n_nodes} nodes')
print(f'DILI genes: {len(dili_genes)}')
print(f'Non-DILI genes: {len(non_dili_genes)}')

# Real compounds for reference
targets_df = pd.read_csv(ecnp.data_dir / 'targets_lcc.csv')
hyp_targets = set(targets_df[targets_df['compound'] == 'Hyperforin']['gene_symbol'].tolist())
que_targets = set(targets_df[targets_df['compound'] == 'Quercetin']['gene_symbol'].tolist())

hyp_direct = hyp_targets & dili_set
que_direct = que_targets & dili_set
print(f'\nHyperforin: {len(hyp_targets)} targets, {len(hyp_direct)} direct DILI hits')
print(f'Quercetin: {len(que_targets)} targets, {len(que_direct)} direct DILI hits')

# =============================================================================
# TEST 1: Synthetic compounds with varying direct hit counts
# =============================================================================
print('\n' + '=' * 70)
print('TEST 1: SYNTHETIC COMPOUNDS WITH VARYING DIRECT HITS')
print('=' * 70)

np.random.seed(42)
results = []

# Generate compounds with 0, 1, 2, 3, 4, 5 direct hits + random non-DILI targets
for n_direct in [0, 1, 2, 3, 4, 5]:
    for trial in range(10):  # 10 trials each
        # Pick n_direct DILI genes + (10 - n_direct) non-DILI genes = 10 total
        n_non_direct = 10 - n_direct
        
        if n_direct > 0:
            direct_targets = list(np.random.choice(dili_genes, n_direct, replace=False))
        else:
            direct_targets = []
        
        non_direct_targets = list(np.random.choice(non_dili_genes, n_non_direct, replace=False))
        targets = direct_targets + non_direct_targets
        
        result = ecnp.compute(targets)
        if result['status'].value == 'success':
            results.append({
                'n_direct': n_direct,
                'trial': trial,
                'Z': result['Z'],
                'I_T': result['I_T']
            })

df = pd.DataFrame(results)

# Compute correlation
corr, pval = stats.pearsonr(df['n_direct'], df['Z'])
print(f'\nCorrelation(direct_hits, Z): r = {corr:.3f}, p = {pval:.2e}')

# Summary by direct hit count
print('\n| Direct Hits | Mean Z | Std Z | n |')
print('|-------------|--------|-------|---|')
for n_direct in sorted(df['n_direct'].unique()):
    subset = df[df['n_direct'] == n_direct]
    print(f'| {n_direct} | {subset["Z"].mean():+.2f} | {subset["Z"].std():.2f} | {len(subset)} |')

# =============================================================================
# TEST 2: ZERO DIRECT HITS - Critical Test
# =============================================================================
print('\n' + '=' * 70)
print('TEST 2: ZERO DIRECT HITS (CRITICAL RWR TEST)')
print('=' * 70)

zero_hit_results = df[df['n_direct'] == 0]
mean_z_zero = zero_hit_results['Z'].mean()
std_z_zero = zero_hit_results['Z'].std()

print(f'\nCompounds with ZERO direct DILI hits:')
print(f'  Mean Z: {mean_z_zero:+.2f}')
print(f'  Std Z: {std_z_zero:.2f}')
print(f'  Range: [{zero_hit_results["Z"].min():.2f}, {zero_hit_results["Z"].max():.2f}]')

# What does this mean?
if mean_z_zero > 0:
    print(f'\n  FINDING: Compounds with NO direct hits still get positive Z-scores!')
    print(f'  This suggests RWR captures pathway proximity effects.')
else:
    print(f'\n  FINDING: Compounds with NO direct hits have Z ~ 0')
    print(f'  This means direct hits dominate; RWR adds minimal value.')

# =============================================================================
# TEST 3: High RWR influence from non-DILI targets
# =============================================================================
print('\n' + '=' * 70)
print('TEST 3: IDENTIFY HIGH-INFLUENCE NON-DILI NODES')
print('=' * 70)

# Find non-DILI genes with highest DILI influence (pathway hubs)
non_dili_influences = [(g, ecnp.m_array[ecnp.node_to_idx[g]]) for g in non_dili_genes]
non_dili_influences.sort(key=lambda x: -x[1])

print('\nTop 10 non-DILI genes by DILI influence:')
print('| Gene | DILI Influence | Percentile |')
print('|------|----------------|------------|')
for g, inf in non_dili_influences[:10]:
    pct = ecnp.percentile_ranks[g] * 100
    print(f'| {g} | {inf:.4f} | {pct:.1f}% |')

# Test a "pathway hub" compound - 10 highest-influence non-DILI genes
hub_targets = [g for g, _ in non_dili_influences[:10]]
hub_result = ecnp.compute(hub_targets)
print(f'\n"Pathway Hub" compound (top 10 non-DILI by influence):')
print(f'  Z-score: {hub_result["Z"]:.2f}')
print(f'  Direct DILI hits: 0')

# Compare to random non-DILI
random_non_dili = list(np.random.choice(non_dili_genes, 10, replace=False))
random_result = ecnp.compute(random_non_dili)
print(f'\nRandom non-DILI compound:')
print(f'  Z-score: {random_result["Z"]:.2f}')
print(f'  Direct DILI hits: 0')

# =============================================================================
# TEST 4: ECNP vs Naive Counting - Full Comparison
# =============================================================================
print('\n' + '=' * 70)
print('TEST 4: ECNP vs NAIVE COUNTING DISCRIMINATION')
print('=' * 70)

# Generate 100 synthetic compounds
np.random.seed(123)
full_results = []

for i in range(100):
    k = np.random.randint(5, 30)
    n_direct = min(np.random.randint(0, 6), k)  # 0-5 direct hits
    n_non = k - n_direct
    
    if n_direct > 0 and n_direct <= len(dili_genes):
        direct_t = list(np.random.choice(dili_genes, n_direct, replace=False))
    else:
        direct_t = []
        n_direct = 0
        n_non = k
    
    non_t = list(np.random.choice(non_dili_genes, n_non, replace=False))
    targets = direct_t + non_t
    
    result = ecnp.compute(targets)
    if result['status'].value == 'success':
        full_results.append({
            'compound': i,
            'k': len(targets),
            'n_direct': n_direct,
            'Z': result['Z'],
            'I_T': result['I_T'],
            'naive_score': n_direct  # Simple counting
        })

full_df = pd.DataFrame(full_results)

# Correlation analysis
corr_ecnp_direct, _ = stats.pearsonr(full_df['n_direct'], full_df['Z'])
corr_ecnp_naive, _ = stats.pearsonr(full_df['naive_score'], full_df['Z'])

print(f'\nCorrelation with direct hit count:')
print(f'  ECNP Z vs direct_hits: r = {corr_ecnp_direct:.3f}')
print(f'  Naive score = direct_hits (by definition)')

# Is there variance in Z unexplained by direct hits?
# Regress Z on n_direct and look at residual variance
from scipy.stats import linregress
slope, intercept, r_value, p_value, std_err = linregress(full_df['n_direct'], full_df['Z'])
residuals = full_df['Z'] - (slope * full_df['n_direct'] + intercept)
residual_var = residuals.var()
total_var = full_df['Z'].var()
unexplained_ratio = residual_var / total_var

print(f'\nVariance analysis:')
print(f'  R² of Z ~ direct_hits: {r_value**2:.3f}')
print(f'  Unexplained variance: {unexplained_ratio:.1%}')

if unexplained_ratio > 0.1:
    print(f'  FINDING: {unexplained_ratio:.1%} of Z variance NOT explained by direct hits')
    print(f'  This is the RWR "pathway effect" contribution.')
else:
    print(f'  FINDING: Direct hits explain {r_value**2:.1%} of Z variance')
    print(f'  RWR adds minimal value beyond counting.')

# =============================================================================
# SUMMARY
# =============================================================================
print('\n' + '=' * 70)
print('SUMMARY: DOES ECNP ADD VALUE BEYOND COUNTING?')
print('=' * 70)

print(f'''
Key Metrics:
  - Correlation(Z, direct_hits) = {corr_ecnp_direct:.2f}
  - R² of Z ~ direct_hits = {r_value**2:.2f}
  - Mean Z for 0 direct hits = {mean_z_zero:+.2f}
  - "Pathway Hub" Z (0 hits) = {hub_result["Z"]:.2f}
  
Interpretation:
''')

if r_value**2 < 0.70:
    print('  ✅ ECNP captures significant pathway effects beyond direct counting')
    print(f'     Only {r_value**2:.0%} of Z variance comes from direct hits')
elif r_value**2 < 0.90:
    print('  ⚠️ ECNP partially captures pathway effects')
    print(f'     {r_value**2:.0%} of Z variance comes from direct hits')
else:
    print('  ❌ ECNP is mostly just counting direct hits')
    print(f'     {r_value**2:.0%} of Z variance explained by counting')

if hub_result["Z"] > 2:
    print(f'  ✅ Pathway hubs generate high Z without direct hits (Z={hub_result["Z"]:.1f})')
else:
    print(f'  ⚠️ Pathway hubs generate low Z without direct hits (Z={hub_result["Z"]:.1f})')

print('=' * 70)

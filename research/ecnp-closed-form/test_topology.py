"""Test: Same direct hits, different network topology"""
import sys
sys.path.insert(0, r'v:\new\h-perforatum-network-tox\research\ecnp-closed-form\src\core')
from pathlib import Path
import numpy as np
from ecnp_optimized import ECNPOptimized

root = Path(r'v:\new\h-perforatum-network-tox')
ecnp = ECNPOptimized(root)

dili_set = set(ecnp.dili_genes)
non_dili = [g for g in ecnp.node_list if g not in dili_set]

# Sort non-DILI by influence
influences = [(g, ecnp.m_array[ecnp.node_to_idx[g]]) for g in non_dili]
influences.sort(key=lambda x: -x[1])
high_inf = [g for g,_ in influences[:100]]  # Top 100 by influence
low_inf = [g for g,_ in influences[-100:]]   # Bottom 100

np.random.seed(42)
dili_list = list(dili_set & set(ecnp.node_list))

print('SAME DIRECT HITS, DIFFERENT TOPOLOGY')
print('='*50)
print('| Direct | Non-DILI Type | Z-score |')
print('|--------|---------------|---------|')

for n_direct in [0, 1, 2, 3]:
    if n_direct > 0:
        direct = list(np.random.choice(dili_list, n_direct, replace=False))
    else:
        direct = []
    
    n_non = 10 - n_direct
    
    # High-influence non-DILI
    hi = direct + list(np.random.choice(high_inf, n_non, replace=False))
    r_hi = ecnp.compute(hi)
    
    # Low-influence non-DILI  
    lo = direct + list(np.random.choice(low_inf, n_non, replace=False))
    r_lo = ecnp.compute(lo)
    
    z_hi = r_hi['Z'] if r_hi['status'].value == 'success' else float('nan')
    z_lo = r_lo['Z'] if r_lo['status'].value == 'success' else float('nan')
    
    print(f'| {n_direct} | High-influence | {z_hi:+.2f} |')
    print(f'| {n_direct} | Low-influence  | {z_lo:+.2f} |')
    diff = z_hi - z_lo if not (np.isnan(z_hi) or np.isnan(z_lo)) else float('nan')
    print(f'| {n_direct} | DIFFERENCE     | {diff:+.2f} |')
    print('|--------|---------------|---------|')

print()
print('If ECNP captures topology, high-influence should score HIGHER')
print('than low-influence with SAME direct hit count.')

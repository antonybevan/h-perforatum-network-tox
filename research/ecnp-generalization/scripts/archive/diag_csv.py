import pandas as pd
from sklearn.metrics import roc_auc_score
from pathlib import Path

ROOT = Path(r'v:\new\h-perforatum-network-tox')
df = pd.read_csv(ROOT / 'research/ecnp-generalization/results/dilirank_706_with_ecnp.csv')
v = df[df.ecnp_z.notna()]

results = []

# Overall
pos = v[v.is_dili==1].ecnp_z
neg = v[v.is_dili==0].ecnp_z
np_count = int(v.is_dili.sum())
nn_count = len(v) - np_count
auc = roc_auc_score(v.is_dili, v.ecnp_z)
results.append({'subset': 'Overall', 'n': len(v), 'dili_pos': np_count, 'dili_neg': nn_count, 
                'mean_z_pos': pos.mean(), 'mean_z_neg': neg.mean(), 'auc': auc})

# By pool size
for pmin, pmax, lbl in [(1,100,'pool_1-100'), (101,300,'pool_101-300'), (301,600,'pool_301-600'), (601,5000,'pool_601+')]:
    s = v[(v.pool_size>=pmin) & (v.pool_size<=pmax)]
    n = len(s)
    np_count = int(s.is_dili.sum())
    nn_count = n - np_count
    if np_count >= 3 and nn_count >= 3:
        auc = roc_auc_score(s.is_dili, s.ecnp_z)
        results.append({'subset': lbl, 'n': n, 'dili_pos': np_count, 'dili_neg': nn_count, 
                        'mean_z_pos': s[s.is_dili==1].ecnp_z.mean(), 'mean_z_neg': s[s.is_dili==0].ecnp_z.mean(), 'auc': auc})

# By target count  
for kmin, kmax, lbl in [(1,2,'k_1-2'), (3,5,'k_3-5'), (6,10,'k_6-10'), (11,100,'k_11+')]:
    s = v[(v.n_targets_mapped>=kmin) & (v.n_targets_mapped<=kmax)]
    n = len(s)
    np_count = int(s.is_dili.sum())
    nn_count = n - np_count
    if np_count >= 3 and nn_count >= 3:
        auc = roc_auc_score(s.is_dili, s.ecnp_z)
        results.append({'subset': lbl, 'n': n, 'dili_pos': np_count, 'dili_neg': nn_count,
                        'mean_z_pos': s[s.is_dili==1].ecnp_z.mean(), 'mean_z_neg': s[s.is_dili==0].ecnp_z.mean(), 'auc': auc})

# Save results
results_df = pd.DataFrame(results)
output = ROOT / 'research/ecnp-generalization/results/performance_gap_diagnosis.csv'
results_df.to_csv(output, index=False)
print(f"Saved: {output}")

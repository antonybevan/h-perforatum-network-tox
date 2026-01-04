import pandas as pd
from sklearn.metrics import roc_auc_score
from pathlib import Path

ROOT = Path(r'v:\new\h-perforatum-network-tox')
df = pd.read_csv(ROOT / 'research/ecnp-generalization/results/dilirank_706_with_ecnp.csv')
v = df[df.ecnp_z.notna()]

print("="*60)
print("DIAGNOSTIC: 706-Drug Performance Gap")
print("="*60)

# Class balance
print(f"\n1. CLASS BALANCE")
print(f"   n = {len(v)}")
print(f"   DILI+: {v.is_dili.sum()} ({v.is_dili.mean()*100:.1f}%)")
print(f"   DILI-: {(v.is_dili==0).sum()} ({(1-v.is_dili.mean())*100:.1f}%)")

# ECNP separation
pos = v[v.is_dili==1].ecnp_z
neg = v[v.is_dili==0].ecnp_z
print(f"\n2. ECNP SEPARATION")
print(f"   DILI+ mean Z: {pos.mean():.3f}")
print(f"   DILI- mean Z: {neg.mean():.3f}")
print(f"   Separation:   {pos.mean()-neg.mean():.3f} (should be >0)")

# AUC by pool size
print(f"\n3. AUC BY POOL SIZE (network coverage)")
for pmin, pmax, lbl in [(1,100,'1-100'), (101,300,'101-300'), (301,600,'301-600'), (601,2000,'601+')]:
    s = v[(v.pool_size>=pmin) & (v.pool_size<=pmax)]
    n = len(s)
    np_count = int(s.is_dili.sum())
    nn_count = n - np_count
    if np_count >= 3 and nn_count >= 3:
        auc = roc_auc_score(s.is_dili, s.ecnp_z)
        print(f"   pool={lbl}: n={n}, {np_count}:{nn_count}, AUC={auc:.3f}")
    else:
        print(f"   pool={lbl}: n={n}, {np_count}:{nn_count}, (imbalanced)")

# AUC by target count
print(f"\n4. AUC BY TARGET COUNT")
for kmin, kmax, lbl in [(1,2,'k=1-2'), (3,5,'k=3-5'), (6,10,'k=6-10'), (11,50,'k>10')]:
    s = v[(v.n_targets_mapped>=kmin) & (v.n_targets_mapped<=kmax)]
    n = len(s)
    np_count = int(s.is_dili.sum())
    nn_count = n - np_count
    if np_count >= 3 and nn_count >= 3:
        auc = roc_auc_score(s.is_dili, s.ecnp_z)
        print(f"   {lbl}: n={n}, {np_count}:{nn_count}, AUC={auc:.3f}")
    else:
        print(f"   {lbl}: n={n}, {np_count}:{nn_count}, (imbalanced)")

# High coverage subset
print(f"\n5. HIGH COVERAGE SUBSET (pool > 200)")
high_cov = v[v.pool_size > 200]
if len(high_cov) >= 20:
    np_count = int(high_cov.is_dili.sum())
    nn_count = len(high_cov) - np_count
    if np_count >= 3 and nn_count >= 3:
        auc = roc_auc_score(high_cov.is_dili, high_cov.ecnp_z)
        print(f"   n={len(high_cov)}, {np_count}:{nn_count}, AUC={auc:.3f}")

print("\n" + "="*60)
print("DIAGNOSIS")
print("="*60)
print("""
The ECNP *does* separate DILI+ from DILI- (mean Z diff = 0.5).
But the model AUC is low (~0.58) which suggests:

1. HIGH CLASS IMBALANCE (67% DILI+) makes AUC misleading
2. Need to check if better coverage = better AUC
3. The 202-drug model may have used higher-coverage drugs only
""")

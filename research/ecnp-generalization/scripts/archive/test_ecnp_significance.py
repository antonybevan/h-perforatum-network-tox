"""
Test: Does ECNP make a STATISTICALLY SIGNIFICANT contribution?
==============================================================

Proper ablation with p-value:
1. Model WITHOUT ECNP
2. Model WITH ECNP
3. Bootstrap/Permutation test for significance
"""
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

ROOT = Path(r'v:\new\h-perforatum-network-tox')

# Load data
df = pd.read_csv(ROOT / 'research/ecnp-generalization/results/dilirank_706_with_ecnp.csv')
df = df[df['ecnp_z'].notna()].copy()

# Features
X_base = df[['n_targets_mapped', 'pool_size']].values  # Base features
X_ecnp = df[['n_targets_mapped', 'pool_size', 'ecnp_z']].values  # + ECNP
y = df['is_dili'].values

print("="*60)
print("ECNP SIGNIFICANCE TEST")
print("="*60)
print(f"Samples: {len(df)}, DILI+: {y.sum()}, DILI-: {(y==0).sum()}")

# Standardize
X_base_scaled = StandardScaler().fit_transform(X_base)
X_ecnp_scaled = StandardScaler().fit_transform(X_ecnp)

# 5-fold CV
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

def get_cv_auc(X, y, cv):
    aucs = []
    for train, test in cv.split(X, y):
        model = LogisticRegression(class_weight='balanced', max_iter=1000)
        model.fit(X[train], y[train])
        y_pred = model.predict_proba(X[test])[:, 1]
        aucs.append(roc_auc_score(y[test], y_pred))
    return np.mean(aucs)

auc_base = get_cv_auc(X_base_scaled, y, cv)
auc_ecnp = get_cv_auc(X_ecnp_scaled, y, cv)
delta = auc_ecnp - auc_base

print(f"\n--- OBSERVED ---")
print(f"AUC without ECNP: {auc_base:.4f}")
print(f"AUC with ECNP:    {auc_ecnp:.4f}")
print(f"Delta:            {delta:+.4f}")

# Permutation test: Is ECNP contribution significant?
print(f"\n--- PERMUTATION TEST (n=1000) ---")
n_perms = 1000
perm_deltas = []

for i in range(n_perms):
    # Shuffle ECNP column
    ecnp_perm = np.random.permutation(df['ecnp_z'].values)
    X_perm = np.column_stack([X_base, ecnp_perm])
    X_perm_scaled = StandardScaler().fit_transform(X_perm)
    
    auc_perm = get_cv_auc(X_perm_scaled, y, cv)
    perm_deltas.append(auc_perm - auc_base)
    
    if (i+1) % 200 == 0:
        print(f"  {i+1}/{n_perms} permutations done...")

perm_deltas = np.array(perm_deltas)

# P-value: how often does random ECNP beat real ECNP?
p_value = (perm_deltas >= delta).mean()

print(f"\n--- RESULTS ---")
print(f"ECNP contribution: {delta:+.4f} AUC")
print(f"Permutation null mean: {perm_deltas.mean():+.4f}")
print(f"Permutation null std:  {perm_deltas.std():.4f}")
print(f"P-value: {p_value:.4f}")

if p_value < 0.05:
    print(f"\n*** ECNP IS SIGNIFICANT (p < 0.05) ***")
else:
    print(f"\n*** ECNP IS NOT SIGNIFICANT (p = {p_value:.3f}) ***")

# Save
results = pd.DataFrame({
    'metric': ['auc_without_ecnp', 'auc_with_ecnp', 'delta', 'p_value', 'significant'],
    'value': [auc_base, auc_ecnp, delta, p_value, p_value < 0.05]
})
output = ROOT / 'research/ecnp-generalization/results/ecnp_significance_test.csv'
results.to_csv(output, index=False)
print(f"\nSaved: {output}")

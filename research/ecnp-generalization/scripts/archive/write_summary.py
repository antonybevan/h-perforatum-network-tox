import pandas as pd
from pathlib import Path

ROOT = Path(r'v:/new/h-perforatum-network-tox')

# Load data
df = pd.read_csv(ROOT / 'research/ecnp-generalization/results/dilirank_706_with_ecnp.csv')
df_valid = df[df['ecnp_z'].notna()]
results = pd.read_csv(ROOT / 'research/ecnp-generalization/results/model_comparison.csv')

# Print summary
print('='*60)
print('706-DRUG ECNP MODEL RESULTS')
print('='*60)
print(f'Total drugs: {len(df)}')
print(f'Valid ECNP: {len(df_valid)}')
print(f'DILI+: {df_valid.is_dili.sum()} ({df_valid.is_dili.mean()*100:.1f}%)')
print(f'DILI-: {(df_valid.is_dili==0).sum()}')
print()
print('Model Performance (5-fold CV AUC):')
print('-'*40)
for _, row in results.iterrows():
    print(f"{row['model']:25s}: {row['auc_mean_706']:.3f} +/- {row['auc_std_706']:.3f}")
best = results.loc[results['auc_mean_706'].idxmax()]
print(f"\nBest: {best['model']} with AUC = {best['auc_mean_706']:.3f}")

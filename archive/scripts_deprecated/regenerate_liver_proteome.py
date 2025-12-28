"""
Regenerate liver_proteome.csv from GTEx source data.
"""

import pandas as pd
import sys
sys.path.append('src')
from pathlib import Path
from network_tox.utils.data_loader import load_liver_genes

# Load directly from GTEx raw file
print('Loading liver genes from GTEx...')
gtex_path = Path('data/raw/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct')
liver_genes = load_liver_genes(gtex_path)

# Save to CSV
output_path = Path('data/processed/liver_proteome.csv')
df = pd.DataFrame({'gene_symbol': sorted(list(liver_genes))})
df.to_csv(output_path, index=False)

print(f'Saved {len(liver_genes)} liver-expressed genes to {output_path}')
print(f'First 10: {list(df["gene_symbol"].head(10))}')
print('Done!')

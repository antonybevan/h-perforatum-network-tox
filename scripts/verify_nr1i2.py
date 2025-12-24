import pandas as pd
import networkx as nx
import sys
sys.path.append('src')
from network_tox.utils.data_loader import load_liver_genes
from network_tox.core.network import filter_to_tissue
from pathlib import Path

# Load processed targets
df_targets = pd.read_csv('data/processed/targets.csv')  # Unfiltered targets
hf_targets = df_targets[df_targets['compound']=='Hyperforin']['gene_name'].tolist()

# Load liver network
gtex_path = Path('data/raw/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct')
liver_genes = load_liver_genes(gtex_path)
df_net = pd.read_parquet('data/processed/network_900.parquet')
G = nx.from_pandas_edgelist(df_net, 'protein1', 'protein2')
G_liver = filter_to_tissue(G, liver_genes)

# Filter
in_liver = [t for t in hf_targets if t in G_liver]

print(f'Total Hyperforin targets in processed file: {len(hf_targets)}')
print(f'Targets in liver network: {len(in_liver)}')
print(f'Targets: {sorted(in_liver)}')
print(f'NR1I2 included: {"NR1I2" in in_liver}')

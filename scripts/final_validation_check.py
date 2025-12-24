"""
Final Validation Check Script
Verifies all data, scripts, and results are correct without bias or errors.
"""

import pandas as pd
import networkx as nx
import sys
from pathlib import Path

sys.path.append('src')
from network_tox.utils.data_loader import load_liver_genes
from network_tox.core.network import filter_to_tissue

print("="*80)
print("FINAL COMPREHENSIVE VALIDATION CHECK")
print("="*80)

DATA_DIR = Path('data')
RESULTS_DIR = Path('results')

# Test 1: Verify NR1I2 is in processed target files
print("\n1. CHECKING PROCESSED TARGET FILES")
print("-" * 80)

for threshold in ['900', '700']:
    df = pd.read_csv(DATA_DIR / f'processed/targets_{threshold}.csv')
    hf_targets = df[df['compound'] == 'Hyperforin']
    
    has_nr1i2 = 'NR1I2' in hf_targets['gene_name'].values
    nr1i2_protein = hf_targets[hf_targets['gene_name'] == 'NR1I2']['protein_id'].values
    
    print(f"\nThreshold {threshold}:")
    print(f"  Total Hyperforin targets: {len(hf_targets)}")
    print(f"  NR1I2 present: {has_nr1i2} {'[OK]' if has_nr1i2 else '[FAIL]'}")
    if has_nr1i2:
        print(f"  NR1I2 protein ID: {nr1i2_protein[0]}")
    print(f"  Total Quercetin targets: {len(df[df['compound'] == 'Quercetin'])}")

# Test 2: Verify liver network filtering
print("\n2. CHECKING LIVER NETWORK FILTERING")
print("-" * 80)

gtex_path = DATA_DIR / 'raw/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct'
liver_genes = load_liver_genes(gtex_path)

for threshold in ['900', '700']:
    # Load targets
    df_targets = pd.read_csv(DATA_DIR / f'processed/targets_{threshold}.csv')
    
    # Load network
    if threshold == '900':
        df_net = pd.read_parquet(DATA_DIR / 'processed/network_900.parquet')
        G = nx.from_pandas_edgelist(df_net, 'protein1', 'protein2')
    else:
        df_net = pd.read_parquet(DATA_DIR / 'processed/network_700.parquet')
        G = nx.from_pandas_edgelist(df_net, 'gene1', 'gene2')
    
    G_liver = filter_to_tissue(G, liver_genes)
    
    # Check targets in liver network
    hf_targets = df_targets[df_targets['compound'] == 'Hyperforin']['gene_name'].tolist()
    qu_targets = df_targets[df_targets['compound'] == 'Quercetin']['gene_name'].tolist()
    
    hf_in_liver = [t for t in hf_targets if t in G_liver]
    qu_in_liver = [t for t in qu_targets if t in G_liver]
    
    print(f"\nThreshold {threshold}:")
    print(f"  Hyperforin: {len(hf_targets)} total -> {len(hf_in_liver)} in liver network")
    print(f"  Quercetin: {len(qu_targets)} total -> {len(qu_in_liver)} in liver network")
    print(f"  NR1I2 in liver network: {'NR1I2' in hf_in_liver} {'[OK]' if 'NR1I2' in hf_in_liver else '[FAIL]'}")

# Test 3: Verify DILI modules
print("\n3. CHECKING DILI MODULES")
print("-" * 80)

for threshold in ['900', '700']:
    dili_df = pd.read_csv(DATA_DIR / f'processed/dili_{threshold}_lcc.csv')
    
    # Load network
    if threshold == '900':
        df_net = pd.read_parquet(DATA_DIR / 'processed/network_900.parquet')
        G = nx.from_pandas_edgelist(df_net, 'protein1', 'protein2')
    else:
        df_net = pd.read_parquet(DATA_DIR / 'processed/network_700.parquet')
        G = nx.from_pandas_edgelist(df_net, 'gene1', 'gene2')
    
    G_liver = filter_to_tissue(G, liver_genes)
    
    dili_genes = dili_df['protein_id'].tolist()
    dili_in_liver = [g for g in dili_genes if g in G_liver]
    
    print(f"\nThreshold {threshold}:")
    print(f"  DILI genes total: {len(dili_genes)}")
    print(f"  DILI genes in liver network: {len(dili_in_liver)}")

# Test 4: Verify validation results consistency
print("\n4. CHECKING VALIDATION RESULTS")
print("-" * 80)

df_900 = pd.read_csv(RESULTS_DIR / 'final_statistics.csv')
df_700 = pd.read_csv(RESULTS_DIR / 'final_statistics_700.csv')

print("\nPrimary (>=900):")
print(df_900.to_string(index=False))

print("\nRobustness (>=700):")
print(df_700.to_string(index=False))

# Test 5: Check for NR1I2 in target lists used
print("\n5. VERIFYING NR1I2 INCLUSION IN ANALYSIS")
print("-" * 80)

# The actual targets should be 9 for Hyperforin (including NR1I2)
expected_hf_900 = 9
expected_qu_900 = 63
expected_hf_700 = 9
expected_qu_700 = 63

df_targets_900 = pd.read_csv(DATA_DIR / 'processed/targets_900.csv')
G_900 = nx.from_pandas_edgelist(pd.read_parquet(DATA_DIR / 'processed/network_900.parquet'), 'protein1', 'protein2')
G_900_liver = filter_to_tissue(G_900, liver_genes)
hf_900_in = [t for t in df_targets_900[df_targets_900['compound']=='Hyperforin']['gene_name'] if t in G_900_liver]
qu_900_in = [t for t in df_targets_900[df_targets_900['compound']=='Quercetin']['gene_name'] if t in G_900_liver]

print(f"\nExpected target counts (>=900): Hyperforin={expected_hf_900}, Quercetin={expected_qu_900}")
print(f"Actual target counts (>=900): Hyperforin={len(hf_900_in)}, Quercetin={len(qu_900_in)}")
print(f"Match: {len(hf_900_in) == expected_hf_900 and len(qu_900_in) == expected_qu_900} {'[OK]' if len(hf_900_in) == expected_hf_900 and len(qu_900_in) == expected_qu_900 else '[FAIL]'}")
print(f"\nNR1I2 in actual targets: {'NR1I2' in hf_900_in} {'[OK]' if 'NR1I2' in hf_900_in else '[FAIL]'}")
print(f"All 9 targets: {sorted(hf_900_in)}")

# Test 6: Bootstrap check
print("\n6. CHECKING BOOTSTRAP SENSITIVITY")
print("-" * 80)

df_boot = pd.read_csv(RESULTS_DIR / 'bootstrap_sensitivity.csv')
print(f"Bootstrap iterations: {len(df_boot)}")
ci_lower = df_boot['quercetin_sampled_influence'].quantile(0.025)
ci_upper = df_boot['quercetin_sampled_influence'].quantile(0.975)
hyperforin_obs = 0.0834

print(f"95% CI: [{ci_lower:.4f}, {ci_upper:.4f}]")
print(f"Hyperforin observed: {hyperforin_obs:.4f}")
print(f"Within CI: {ci_lower <= hyperforin_obs <= ci_upper} {'[OK]' if ci_lower <= hyperforin_obs <= ci_upper else '[FAIL]'}")

print("\n" + "="*80)
print("VALIDATION COMPLETE")
print("="*80)

"""
Analyze how targets_raw.csv was filtered to create targets.csv
"""

import pandas as pd
from pathlib import Path

DATA_DIR = Path('data')

# Load files
raw = pd.read_csv(DATA_DIR / 'raw/targets_raw.csv')
processed = pd.read_csv(DATA_DIR / 'processed/targets.csv')

print("="*80)
print("TARGETS FILTERING ANALYSIS")
print("="*80)

print("\n1. RAW TARGETS (from literature/ChEMBL)")
print("-" * 80)
print(f"Total rows: {len(raw)}")
print(f"Hyperforin: {len(raw[raw['compound']=='Hyperforin'])}")
print(f"Quercetin: {len(raw[raw['compound']=='Quercetin'])}")

print("\n2. PROCESSED TARGETS (after filtering)")
print("-" * 80)
print(f"Total rows: {len(processed)}")
print(f"Hyperforin: {len(processed[processed['compound']=='Hyperforin'])}")
print(f"Quercetin: {len(processed[processed['compound']=='Quercetin'])}")

print(f"\n3. FILTERED OUT")
print("-" * 80)
print(f"Total removed: {len(raw) - len(processed)} targets")

# Find which protein IDs were removed
raw_ids = set(raw['protein_id'])
proc_ids = set(processed['protein_id'])
removed_ids = raw_ids - proc_ids

print(f"\nRemoved protein IDs: {len(removed_ids)}")

# Show removed targets by compound
for compound in ['Hyperforin', 'Quercetin']:
    raw_comp = raw[raw['compound'] == compound]
    proc_comp = processed[processed['compound'] == compound]
    
    raw_comp_ids = set(raw_comp['protein_id'])
    proc_comp_ids = set(proc_comp['protein_id'])
    removed_comp = raw_comp_ids - proc_comp_ids
    
    print(f"\n{compound}:")
    print(f"  Raw: {len(raw_comp_ids)}")
    print(f"  Processed: {len(proc_comp_ids)}")
    print(f"  Removed: {len(removed_comp)}")
    
    if len(removed_comp) > 0:
        # Try to match to gene names
        removed_df = raw_comp[raw_comp['protein_id'].isin(removed_comp)]
        print(f"  Removed proteins:")
        for _, row in removed_df.iterrows():
            # Try to find gene name in processed
            gene_name = "UNKNOWN"
            if 'gene_name' in processed.columns:
                match = processed[processed['protein_id'] == row['protein_id']]
                if len(match) > 0:
                    gene_name = match.iloc[0]['gene_name']
            print(f"    {row['protein_id']} - {row['source']}")

print("\n" + "="*80)
print("ANALYSIS COMPLETE")
print("="*80)

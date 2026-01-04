"""
Map UniProt IDs to Gene Symbols
================================

Fetch gene symbols from UniProt API for all unique DrugBank targets.
Then recompute ECNP for 706 drugs using gene symbol nodes.
"""
import pandas as pd
import numpy as np
from pathlib import Path
import requests
import time
import ast
import warnings
warnings.filterwarnings('ignore')

ROOT = Path(r'v:\new\h-perforatum-network-tox')

print("="*70)
print("MAPPING UNIPROT TO GENE SYMBOLS")
print("="*70)

# =============================================================================
# LOAD UNIQUE TARGETS
# =============================================================================

print("\n--- Loading Targets ---")

df = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'results' / 'dilirank_with_targets.csv')

# Get unique UniProt IDs
unique_targets = set()
for targets in df['targets']:
    try:
        t_list = ast.literal_eval(targets)
        unique_targets.update(t_list)
    except:
        pass

print(f"Unique UniProt IDs: {len(unique_targets)}")

# Load existing mapping
existing_map = {}
existing_path = ROOT / 'data' / 'external' / 'uniprot_mapping.csv'
if existing_path.exists():
    with open(existing_path, 'r') as f:
        for line in f:
            if line.startswith('#') or ',' not in line:
                continue
            parts = line.strip().split(',')
            if len(parts) == 2:
                existing_map[parts[0]] = parts[1]
    print(f"Existing mappings: {len(existing_map)}")

# =============================================================================
# FETCH FROM UNIPROT API
# =============================================================================

print("\n--- Fetching from UniProt API ---")

def fetch_gene_symbol(uniprot_id):
    """Fetch gene symbol from UniProt."""
    try:
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
        r = requests.get(url, timeout=10)
        if r.ok:
            data = r.json()
            # Get gene name
            genes = data.get('genes', [])
            if genes:
                # Primary gene name
                gene_name = genes[0].get('geneName', {}).get('value')
                if gene_name:
                    return gene_name
            return None
        return None
    except:
        return None

# Only fetch for IDs not in existing map
to_fetch = unique_targets - set(existing_map.keys())
print(f"Need to fetch: {len(to_fetch)}")

new_map = {}
errors = []

for i, uid in enumerate(to_fetch):
    if i > 0 and i % 100 == 0:
        print(f"  Progress: {i}/{len(to_fetch)} | Found: {len(new_map)}")
    
    gene = fetch_gene_symbol(uid)
    if gene:
        new_map[uid] = gene
    else:
        errors.append(uid)
    
    time.sleep(0.15)  # Rate limit

print(f"\nFetched: {len(new_map)} new mappings")
print(f"Errors: {len(errors)}")

# =============================================================================
# COMBINE AND SAVE
# =============================================================================

print("\n--- Saving Combined Mapping ---")

full_map = {**existing_map, **new_map}
print(f"Total mappings: {len(full_map)}")

# Save
output = ROOT / 'research' / 'ecnp-generalization' / 'results' / 'uniprot_to_gene.csv'
with open(output, 'w') as f:
    f.write("uniprot_id,gene_symbol\n")
    for uid, gene in sorted(full_map.items()):
        f.write(f"{uid},{gene}\n")
print(f"Saved: {output}")

# =============================================================================
# CHECK OVERLAP WITH NETWORK
# =============================================================================

print("\n--- Checking Network Overlap ---")

import sys
sys.path.insert(0, str(ROOT / 'research' / 'ecnp-generalization' / 'scripts'))
try:
    from ecnp_robust import RobustECNP
    ecnp = RobustECNP()
    lcc_nodes = set(ecnp.node_to_idx.keys())
    
    mapped_genes = set(full_map.values())
    overlap = mapped_genes & lcc_nodes
    
    print(f"Mapped genes: {len(mapped_genes)}")
    print(f"In network: {len(overlap)}")
    print(f"Overlap: {len(overlap)/len(mapped_genes)*100:.1f}%")
except Exception as e:
    print(f"Could not load RobustECNP: {e}")

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "="*70)
print("SUMMARY")
print("="*70)

print(f"""
Unique UniProt IDs: {len(unique_targets)}
Existing mappings: {len(existing_map)}
New mappings fetched: {len(new_map)}
Total mappings: {len(full_map)}

Coverage: {len(full_map)/len(unique_targets)*100:.1f}%
""")

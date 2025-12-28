#!/usr/bin/env python3
"""Check if targets.csv was correctly filtered from raw targets."""

import pandas as pd
import networkx as nx
from pathlib import Path

print("=" * 80)
print(" TARGET FILTERING VERIFICATION")
print("=" * 80)

# Load raw targets
print("\n[1/4] Loading raw targets...")
raw = pd.read_csv('data/raw/targets_raw.csv')
hyp_raw = sorted(set(raw[raw['compound'] == 'Hyperforin']['gene_name']))
quer_raw = sorted(set(raw[raw['compound'] == 'Quercetin']['gene_name']))
print(f"  Raw Hyperforin: {len(hyp_raw)} targets")
print(f"  Raw Quercetin: {len(quer_raw)} targets")

# Load networks
print("\n[2/4] Loading networks...")
df700 = pd.read_parquet('data/processed/network_700.parquet')
G700 = nx.from_pandas_edgelist(df700, 'gene1', 'gene2')
print(f"  Network 700: {G700.number_of_nodes()} nodes, {G700.number_of_edges()} edges")

df900 = pd.read_parquet('data/processed/network_900.parquet')
G900 = nx.from_pandas_edgelist(df900, 'protein1', 'protein2')
print(f"  Network 900: {G900.number_of_nodes()} nodes, {G900.number_of_edges()} edges")

# Check filtering
print("\n[3/4] Checking target presence in networks...")

print("\nHyperforin:")
hyp_in_700 = [t for t in hyp_raw if t in G700]
hyp_in_900 = [t for t in hyp_raw if t in G900]
hyp_in_both = [t for t in hyp_raw if t in G700 and t in G900]
print(f"  In network 700: {len(hyp_in_700)}/{len(hyp_raw)}")
print(f"  In network 900: {len(hyp_in_900)}/{len(hyp_raw)}")
print(f"  In BOTH: {len(hyp_in_both)}/{len(hyp_raw)}")
print(f"  Targets in both: {hyp_in_both}")

print("\nQuercetin:")
quer_in_700 = [t for t in quer_raw if t in G700]
quer_in_900 = [t for t in quer_raw if t in G900]
quer_in_both = [t for t in quer_raw if t in G700 and t in G900]
print(f"  In network 700: {len(quer_in_700)}/{len(quer_raw)}")
print(f"  In network 900: {len(quer_in_900)}/{len(quer_raw)}")
print(f"  In BOTH: {len(quer_in_both)}/{len(quer_raw)}")

# Load processed targets
print("\n[4/4] Comparing with targets.csv...")
proc = pd.read_csv('data/processed/targets.csv')
hyp_proc = sorted(set(proc[proc['compound'] == 'Hyperforin']['gene_name']))
quer_proc = sorted(set(proc[proc['compound'] == 'Quercetin']['gene_name']))

print(f"  Processed Hyperforin: {len(hyp_proc)} targets")
print(f"  Processed Quercetin: {len(quer_proc)} targets")

# Verification
print("\n" + "=" * 80)
print(" VERIFICATION RESULTS")
print("=" * 80)

hyp_match = set(hyp_in_both) == set(hyp_proc)
quer_match = set(quer_in_both) == set(quer_proc)

print(f"\nHyperforin:")
print(f"  Expected (in both networks): {len(hyp_in_both)}")
print(f"  Got (targets.csv): {len(hyp_proc)}")
print(f"  Status: {'✓ CORRECT' if hyp_match else '✗ INCORRECT'}")

if not hyp_match:
    missing = set(hyp_in_both) - set(hyp_proc)
    extra = set(hyp_proc) - set(hyp_in_both)
    if missing:
        print(f"  Missing from targets.csv: {missing}")
    if extra:
        print(f"  Extra in targets.csv: {extra}")

print(f"\nQuercetin:")
print(f"  Expected (in both networks): {len(quer_in_both)}")
print(f"  Got (targets.csv): {len(quer_proc)}")
print(f"  Status: {'✓ CORRECT' if quer_match else '✗ INCORRECT'}")

if not quer_match:
    missing = set(quer_in_both) - set(quer_proc)
    extra = set(quer_proc) - set(quer_in_both)
    if missing:
        print(f"  Missing from targets.csv: {sorted(missing)}")
    if extra:
        print(f"  Extra in targets.csv: {sorted(extra)}")

print("\n" + "=" * 80)
if hyp_match and quer_match:
    print(" ✓✓✓ FILTERING WAS CORRECT")
    print(" targets.csv contains exactly the targets present in BOTH networks")
else:
    print(" ✗✗✗ FILTERING HAS ISSUES")
    print(" targets.csv does not match intersection of both networks")
print("=" * 80)

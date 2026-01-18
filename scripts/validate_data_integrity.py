#!/usr/bin/env python3
"""
Data Integrity Validation

Verifies traceability of data transformations.
"""

import pandas as pd
import networkx as nx
from pathlib import Path
import sys

DATA_DIR = Path('data')
ERRORS = []
WARNINGS = []

def error(msg):
    ERRORS.append(msg)
    print(f"ERROR: {msg}")

def warning(msg):
    WARNINGS.append(msg)
    print(f"WARNING: {msg}")

def ok(msg):
    print(f"OK: {msg}")

def section(title):
    print(f"\n{'='*60}")
    print(f" {title}")
    print(f"{'='*60}")

# Load all data files
section("LOADING DATA FILES")

raw_targets = pd.read_csv(DATA_DIR / 'raw' / 'targets_raw.csv')
proc_targets = pd.read_csv(DATA_DIR / 'processed' / 'targets.csv')
lcc_targets = pd.read_csv(DATA_DIR / 'processed' / 'targets_lcc.csv')

dili_raw = pd.read_csv(DATA_DIR / 'raw' / 'dili_genes_raw.csv')
dili_700 = pd.read_csv(DATA_DIR / 'processed' / 'dili_700_lcc.csv')
dili_900 = pd.read_csv(DATA_DIR / 'processed' / 'dili_900_lcc.csv')

liver = pd.read_csv(DATA_DIR / 'processed' / 'liver_proteome.csv')

n700 = pd.read_parquet(DATA_DIR / 'processed' / 'network_700.parquet')
n900 = pd.read_parquet(DATA_DIR / 'processed' / 'network_900.parquet')
n700_lcc = pd.read_parquet(DATA_DIR / 'processed' / 'network_700_liver_lcc.parquet')
n900_lcc = pd.read_parquet(DATA_DIR / 'processed' / 'network_900_liver_lcc.parquet')

# Parse mapping
mapping_lines = open(DATA_DIR / 'external' / 'uniprot_mapping.csv').readlines()
mapping = {}
for line in mapping_lines:
    if ',' in line and not line.startswith('#'):
        parts = line.strip().split(',')
        if len(parts) == 2:
            mapping[parts[0]] = parts[1]

ok(f"Loaded all files successfully")

# ============================================================
section("TARGETS VALIDATION")
# ============================================================

# Check raw counts
raw_hyp = len(raw_targets[raw_targets['compound'] == 'Hyperforin'])
raw_quer = len(raw_targets[raw_targets['compound'] == 'Quercetin'])
if raw_hyp != 14:
    error(f"Raw Hyperforin count should be 14, got {raw_hyp}")
else:
    ok(f"Raw Hyperforin: {raw_hyp}")
if raw_quer != 122:
    error(f"Raw Quercetin count should be 122, got {raw_quer}")
else:
    ok(f"Raw Quercetin: {raw_quer}")

# Check processed counts
proc_hyp = len(proc_targets[proc_targets['compound'] == 'Hyperforin'])
proc_quer = len(proc_targets[proc_targets['compound'] == 'Quercetin'])
if proc_hyp != 14:
    error(f"Processed Hyperforin count should be 14, got {proc_hyp}")
else:
    ok(f"Processed Hyperforin: {proc_hyp}")
if proc_quer != 87:
    error(f"Processed Quercetin count should be 87, got {proc_quer}")
else:
    ok(f"Processed Quercetin: {proc_quer}")

# Check LCC counts
lcc_hyp = len(lcc_targets[lcc_targets['compound'] == 'Hyperforin'])
lcc_quer = len(lcc_targets[lcc_targets['compound'] == 'Quercetin'])
if lcc_hyp != 10:
    error(f"LCC Hyperforin count should be 10, got {lcc_hyp}")
else:
    ok(f"LCC Hyperforin: {lcc_hyp}")
if lcc_quer != 62:
    error(f"LCC Quercetin count should be 62, got {lcc_quer}")
else:
    ok(f"LCC Quercetin: {lcc_quer}")

# Verify all processed proteins are in raw
proc_pids = set(proc_targets['protein_id'])
raw_pids = set(raw_targets['protein_id'])
if not proc_pids.issubset(raw_pids):
    missing = proc_pids - raw_pids
    error(f"Processed has proteins not in raw: {missing}")
else:
    ok("All processed proteins are in raw")

# Verify all processed proteins have mapping
for pid in proc_pids:
    if pid not in mapping:
        error(f"Processed protein {pid} has no mapping")
ok("All processed proteins have mapping")

# Verify LCC genes are in processed
lcc_genes = set(lcc_targets['gene_symbol'])
proc_genes = set(proc_targets['gene_name'])
if not lcc_genes.issubset(proc_genes):
    missing = lcc_genes - proc_genes
    error(f"LCC has genes not in processed: {missing}")
else:
    ok("All LCC genes are in processed")

# ============================================================
section("FILTERING LOGIC VALIDATION")
# ============================================================

# Trace exactly what was lost
raw_to_proc_lost = raw_pids - proc_pids
expected_lost = 35  # 25 no mapping + 10 non-human
if len(raw_to_proc_lost) != expected_lost:
    error(f"Expected {expected_lost} lost from raw->proc, got {len(raw_to_proc_lost)}")
else:
    ok(f"Raw→Processed: lost {len(raw_to_proc_lost)} proteins (correct)")

# Verify non-human filtering
NON_HUMAN = {'Q91WR5', 'Q63344', 'Q9D6N1', 'Q965D5', 'Q965D6', 'Q965D7', 'P0AES6', 'P03468', 'P0DTD1'}
for pid in NON_HUMAN:
    if pid in proc_pids:
        error(f"Non-human protein {pid} should not be in processed")
ok("All non-human proteins correctly excluded")

# Verify the 10 human proteins are now INCLUDED
SHOULD_BE_INCLUDED = {'P08183', 'P15692'}  # ABCB1 and VEGFA (in processed, may be excluded from LCC)
for pid in SHOULD_BE_INCLUDED:
    if pid not in proc_pids:
        error(f"Human protein {pid} should be in processed")
ok("ABCB1 and VEGFA correctly included in processed")

# ============================================================
section("LCC FILTERING VALIDATION")
# ============================================================

# Build LCC network
G700_lcc = nx.from_pandas_edgelist(n700_lcc, 'gene1', 'gene2')
G900_lcc = nx.from_pandas_edgelist(n900_lcc, 'gene1', 'gene2')

lcc_700_nodes = set(G700_lcc.nodes())
lcc_900_nodes = set(G900_lcc.nodes())
lcc_both = lcc_700_nodes & lcc_900_nodes

# Verify Hyperforin LCC targets are in network
hyp_lcc_genes = set(lcc_targets[lcc_targets['compound'] == 'Hyperforin']['gene_symbol'])
for gene in hyp_lcc_genes:
    if gene not in lcc_both:
        error(f"Hyperforin LCC target {gene} not in liver LCC")
ok("All Hyperforin LCC targets are in liver network")

# Verify excluded Hyperforin targets are NOT in LCC
hyp_proc_genes = set(proc_targets[proc_targets['compound'] == 'Hyperforin']['gene_name'])
hyp_excluded = hyp_proc_genes - hyp_lcc_genes
expected_excluded = {'GRIN1', 'PMAIP1', 'TRPC6', 'VEGFA'}
if hyp_excluded != expected_excluded:
    error(f"Expected excluded Hyperforin: {expected_excluded}, got {hyp_excluded}")
else:
    ok(f"Hyperforin excluded targets correct: {hyp_excluded}")

# Verify excluded are truly not in LCC
for gene in hyp_excluded:
    if gene in lcc_both:
        error(f"Excluded gene {gene} is actually in LCC!")
ok("All excluded Hyperforin targets correctly not in liver LCC")

# ============================================================
section("DILI VALIDATION")
# ============================================================

if len(dili_raw) != 127:
    error(f"DILI raw should be 127, got {len(dili_raw)}")
else:
    ok(f"DILI raw: {len(dili_raw)}")

if len(dili_700) != 84:
    error(f"DILI 700 LCC should be 84, got {len(dili_700)}")
else:
    ok(f"DILI 700 LCC: {len(dili_700)}")

if len(dili_900) != 82:
    error(f"DILI 900 LCC should be 82, got {len(dili_900)}")
else:
    ok(f"DILI 900 LCC: {len(dili_900)}")

# Verify DILI LCC genes are in network
dili_700_genes = set(dili_700['gene_name'])
dili_900_genes = set(dili_900['gene_name'])

for gene in dili_700_genes:
    if gene not in lcc_700_nodes:
        error(f"DILI 700 gene {gene} not in network 700 LCC")

for gene in dili_900_genes:
    if gene not in lcc_900_nodes:
        error(f"DILI 900 gene {gene} not in network 900 LCC")

ok("All DILI LCC genes are in respective networks")

# ============================================================
section("NETWORK VALIDATION")
# ============================================================

if len(n700) != 236712:
    warning(f"Network 700 edges: {len(n700)} (expected 236712)")
else:
    ok(f"Network 700: {len(n700)} edges")

if len(n900) != 100383:
    warning(f"Network 900 edges: {len(n900)} (expected 100383)")
else:
    ok(f"Network 900: {len(n900)} edges")

if G700_lcc.number_of_nodes() != 9773:
    error(f"Network 700 LCC nodes: {G700_lcc.number_of_nodes()} (expected 9773)")
else:
    ok(f"Network 700 LCC: {G700_lcc.number_of_nodes()} nodes")

if G900_lcc.number_of_nodes() != 7677:
    error(f"Network 900 LCC nodes: {G900_lcc.number_of_nodes()} (expected 7677)")
else:
    ok(f"Network 900 LCC: {G900_lcc.number_of_nodes()} nodes")

# Verify 900 LCC is subset of 700 LCC
if not lcc_900_nodes.issubset(lcc_700_nodes):
    extra = lcc_900_nodes - lcc_700_nodes
    error(f"900 LCC has nodes not in 700 LCC: {extra}")
else:
    ok("900 LCC is strict subset of 700 LCC")

# ============================================================
section("LIVER PROTEOME VALIDATION")
# ============================================================

if len(liver) != 13496:
    warning(f"Liver proteome: {len(liver)} genes (expected 13496)")
else:
    ok(f"Liver proteome: {len(liver)} genes")

# Verify LCC nodes are in liver proteome
liver_genes = set(liver['gene_symbol'])
missing_from_liver = lcc_700_nodes - liver_genes
if missing_from_liver:
    error(f"LCC nodes not in liver proteome: {len(missing_from_liver)}")
else:
    ok("All LCC nodes are in liver proteome")

# ============================================================
section("FINAL SUMMARY")
# ============================================================

print()
if ERRORS:
    print(f"\nFAILED: {len(ERRORS)} errors found")
    for e in ERRORS:
        print(f"   - {e}")
    sys.exit(1)
else:
    print(f"\nALL CHECKS PASSED")
    if WARNINGS:
        print(f"   ({len(WARNINGS)} warnings)")
    sys.exit(0)

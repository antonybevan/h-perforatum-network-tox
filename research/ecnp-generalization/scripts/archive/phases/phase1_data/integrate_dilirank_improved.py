"""
Improved DILIrank Integration - Fuzzy Matching
===============================================

Uses fuzzy string matching to maximize DrugBank ↔ DILIrank overlap.
Target: n≥150 labeled compounds.
"""
import pandas as pd
import numpy as np
from pathlib import Path
import re
from difflib import SequenceMatcher

RESEARCH_ROOT = Path(r'v:\new\h-perforatum-network-tox')
DILIRANK_FILE = RESEARCH_ROOT / 'results' / 'tables' / 'dilirank_reference_set.csv'
DRUGBANK_FILE = RESEARCH_ROOT / 'research' / 'ecnp-generalization' / 'data' / 'curated' / 'drugbank_compounds.csv'
OUTPUT_DIR = RESEARCH_ROOT / 'research' / 'ecnp-generalization' / 'data' / 'curated'

# Load DILIrank
print("Loading DILIrank 2.0...")
dilirank = pd.read_csv(DILIRANK_FILE)
print(f"DILIrank entries: {len(dilirank)}")
print(f"Columns: {dilirank.columns.tolist()}")

# Check what label columns exist
label_cols = [c for c in dilirank.columns if 'dili' in c.lower() or 'concern' in c.lower() or 'class' in c.lower()]
print(f"Label columns: {label_cols}")

# Load DrugBank compounds
print("\nLoading DrugBank compounds...")
drugbank = pd.read_csv(DRUGBANK_FILE)
print(f"DrugBank compounds with ≥3 targets: {len(drugbank)}")

# Normalize names for matching
def normalize_name(name):
    """Normalize drug name for fuzzy matching."""
    if pd.isna(name):
        return ''
    s = str(name).lower()
    # Remove common salts/forms
    s = re.sub(r'\s*(hydrochloride|sulfate|sodium|potassium|acetate|maleate|fumarate|citrate|tartrate|phosphate|chloride)\s*$', '', s)
    # Remove special chars
    s = re.sub(r'[^a-z0-9]', '', s)
    return s

dilirank['name_norm'] = dilirank['name'].apply(normalize_name)
drugbank['name_norm'] = drugbank['drug_name'].apply(normalize_name)

# Create lookup from normalized DILIrank names
dilirank_lookup = {}
for idx, row in dilirank.iterrows():
    norm = row['name_norm']
    if norm:
        dilirank_lookup[norm] = row

# Match DrugBank to DILIrank
print("\nMatching DrugBank → DILIrank...")

def fuzzy_match(name_norm, lookup, threshold=0.85):
    """Find best fuzzy match above threshold."""
    if name_norm in lookup:
        return name_norm, 1.0
    
    best_match = None
    best_score = 0
    for candidate in lookup.keys():
        score = SequenceMatcher(None, name_norm, candidate).ratio()
        if score > best_score and score >= threshold:
            best_score = score
            best_match = candidate
    
    return best_match, best_score

matches = []
for idx, row in drugbank.iterrows():
    match_name, score = fuzzy_match(row['name_norm'], dilirank_lookup, threshold=0.80)
    
    if match_name:
        dili_row = dilirank_lookup[match_name]
        matches.append({
            'drugbank_id': row['drugbank_id'],
            'drug_name': row['drug_name'],
            'dilirank_name': dili_row['name'],
            'n_targets': row['n_targets'],
            'targets': row['targets'],
            'match_score': score,
            **{c: dili_row.get(c, None) for c in dilirank.columns if c not in ['name', 'name_norm']}
        })

matches_df = pd.DataFrame(matches)
print(f"Matched: {len(matches_df)} / {len(drugbank)} DrugBank compounds")

# Determine DILI label column
if 'ltkb_dili_class' in matches_df.columns:
    label_col = 'ltkb_dili_class'
elif 'dili_class' in matches_df.columns:
    label_col = 'dili_class'
else:
    # Find it
    for c in matches_df.columns:
        if 'class' in c.lower() or 'concern' in c.lower():
            label_col = c
            break
    else:
        label_col = None
        print("WARNING: Could not find DILI label column!")

if label_col:
    print(f"\nUsing label column: {label_col}")
    print(f"Label distribution:")
    print(matches_df[label_col].value_counts())
    
    # Create binary label
    # Positive = vMost-DILI-Concern, vLess-DILI-Concern
    # Negative = vNo-DILI-Concern
    # Exclude = Ambiguous
    
    matches_df['is_dili'] = matches_df[label_col].apply(lambda x: 
        1 if 'Most' in str(x) or 'Less' in str(x) else 
        (0 if 'No' in str(x) else np.nan)
    )
    
    # Filter to non-ambiguous
    labeled = matches_df[matches_df['is_dili'].notna()].copy()
    labeled['is_dili'] = labeled['is_dili'].astype(int)
    
    print(f"\nLabeled (non-ambiguous): {len(labeled)}")
    print(f"  DILI+ (Most + Less): {(labeled['is_dili'] == 1).sum()}")
    print(f"  DILI- (No-DILI): {(labeled['is_dili'] == 0).sum()}")
else:
    labeled = matches_df.copy()
    print("No label column found - saving all matches")

# Save
output_file = OUTPUT_DIR / 'labeled_compounds_improved.csv'
labeled.to_csv(output_file, index=False)
print(f"\nSaved: {output_file}")

# Summary
print("\n" + "="*60)
print("DATA PLUMBING SUMMARY")
print("="*60)
print(f"DrugBank compounds (≥3 targets): {len(drugbank)}")
print(f"Matched to DILIrank: {len(matches_df)}")
print(f"Non-ambiguous labels: {len(labeled)}")
if label_col:
    n_pos = (labeled['is_dili'] == 1).sum()
    n_neg = (labeled['is_dili'] == 0).sum()
    print(f"  Positives (DILI+): {n_pos}")
    print(f"  Negatives (DILI-): {n_neg}")
    target_reached = len(labeled) >= 150
    print(f"\nTarget n≥150: {'YES' if target_reached else 'NO'} ({len(labeled)})")

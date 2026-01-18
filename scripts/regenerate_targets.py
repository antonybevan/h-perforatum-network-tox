#!/usr/bin/env python3
"""
Regenerate targets.csv from raw data.
"""

import pandas as pd
from pathlib import Path

DATA_DIR = Path('data')

# Standard gene name mapping (use HGNC symbols)
GENE_STANDARDIZATION = {
    'MDR1': 'ABCB1',  # MDR1 is alias for ABCB1
}

# Known non-human proteins by UniProt ID or gene
NON_HUMAN_PATTERNS = {
    # Mouse UniProt prefixes
    'Q91',   # Mouse
    'Q9D',   # Mouse
    'Q63',   # Rat
    'Q965D', # Rat UGT family
}

NON_HUMAN_GENES = {
    'GyrA',     # E. coli DNA gyrase
    'HPV16E6', # Human papillomavirus (not human protein)
    'NSP5',     # SARS-CoV-2
    'NANA',     # Clostridium perfringens sialidase
    'Insr',     # Mouse insulin receptor
    'Rrm2',     # Mouse ribonucleoside reductase
    'Ugt1a6',   # Rat UGT
    'Ugt1a7',   # Rat UGT
    'Ugt1a8',   # Mouse UGT
    'Ugt1a9',   # Rat UGT
}


def load_mapping():
    """Load UniProt to gene symbol mapping."""
    mapping = {}
    mapping_file = DATA_DIR / 'external' / 'uniprot_mapping.csv'
    
    with open(mapping_file) as f:
        for line in f:
            if line.startswith('#') or ',' not in line:
                continue
            parts = line.strip().split(',')
            if len(parts) == 2:
                protein_id, gene = parts
                # Standardize gene name
                gene = GENE_STANDARDIZATION.get(gene, gene)
                mapping[protein_id] = gene
    
    return mapping


def is_non_human(protein_id: str, gene: str) -> bool:
    """Check if protein is non-human based on ID pattern or gene name."""
    # Check UniProt ID prefix
    for prefix in NON_HUMAN_PATTERNS:
        if protein_id.startswith(prefix):
            return True
    
    # Check gene name
    if gene in NON_HUMAN_GENES:
        return True
    
    # Lowercase gene names often indicate non-human
    if gene and gene[0].islower():
        return True
    
    return False


def main():
    print("Regenerating targets from raw data")
    
    # Load raw targets
    raw_file = DATA_DIR / 'raw' / 'targets_raw.csv'
    raw = pd.read_csv(raw_file)
    print(f"\n[1/4] Loaded raw targets: {len(raw)}")
    print(f"      Hyperforin: {len(raw[raw['compound']=='Hyperforin'])}")
    print(f"      Quercetin: {len(raw[raw['compound']=='Quercetin'])}")
    
    # Load mapping
    print("\n[2/4] Loading UniProt -> Gene mapping...")
    mapping = load_mapping()
    print(f"      Mapping entries: {len(mapping)}")
    
    # Apply filters
    print("\n[3/4] Applying filters...")
    results = []
    stats = {'no_mapping': 0, 'non_human': 0, 'kept': 0}
    
    for _, row in raw.iterrows():
        pid = row['protein_id']
        compound = row['compound']
        source = row['source']
        
        # Filter 1: Must have mapping
        if pid not in mapping:
            stats['no_mapping'] += 1
            continue
        
        gene = mapping[pid]
        
        # Filter 2: Must be human
        if is_non_human(pid, gene):
            stats['non_human'] += 1
            continue
        
        results.append({
            'compound': compound,
            'protein_id': pid,
            'source': source,
            'gene_name': gene
        })
        stats['kept'] += 1
    
    print(f"      No mapping: {stats['no_mapping']}")
    print(f"      Non-human: {stats['non_human']}")
    print(f"      Kept: {stats['kept']}")
    
    # Create DataFrame
    df = pd.DataFrame(results)
    
    # Remove duplicates (same compound + protein)
    before_dedup = len(df)
    df = df.drop_duplicates(subset=['compound', 'protein_id'], keep='first')
    print(f"      Duplicates removed: {before_dedup - len(df)}")
    
    # Save
    output_file = DATA_DIR / 'processed' / 'targets.csv'
    df.to_csv(output_file, index=False)
    print(f"\n[4/4] Saved: {output_file}")
    
    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"\nFinal counts:")
    print(f"  Hyperforin: {len(df[df['compound']=='Hyperforin'])}")
    print(f"  Quercetin: {len(df[df['compound']=='Quercetin'])}")
    print(f"  Total: {len(df)}")
    
    print("\nHyperforin targets:")
    hyp = df[df['compound']=='Hyperforin'][['protein_id', 'gene_name']]
    for _, row in hyp.iterrows():
        print(f"  {row['protein_id']} -> {row['gene_name']}")
    
    print("\n" + "=" * 80)


if __name__ == '__main__':
    main()

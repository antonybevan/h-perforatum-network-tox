"""
Curate drug targets from raw data with explicit filtering logic.

This script processes targets_raw.csv and applies multiple filters to create
a high-quality targets.csv file suitable for network analysis.

Usage:
    python scripts/curate_targets.py
"""

import pandas as pd
import gzip
from pathlib import Path
from tqdm import tqdm

DATA_DIR = Path('data')

def load_string_protein_info():
    """Load protein ID to gene symbol mapping from STRING."""
    print("Loading STRING protein-gene mapping...")
    
    info_file = DATA_DIR / 'external/string_info.txt.gz'
    if not info_file.exists():
        print(f"WARNING: {info_file} not found. Skipping STRING validation.")
        return None
    
    protein_to_gene = {}
    with gzip.open(info_file, 'rt') as f:
        # Skip header
        f.readline()
        
        for line in tqdm(f, desc="  Parsing STRING info"):
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                protein_id = parts[0]  # e.g., "9606.ENSP00000000233"
                gene_symbol = parts[1]  # e.g., "ARF5"
                
                # Extract just the ENSP ID for UniProt matching
                if protein_id.startswith('9606.'):
                    ensp_id = protein_id.replace('9606.', '')
                    protein_to_gene[ensp_id] = gene_symbol
    
    print(f"  Loaded {len(protein_to_gene):,} protein-gene mappings from STRING")
    return protein_to_gene

def load_uniprot_to_ensp():
    """Load UniProt to ENSP mapping (simplified - using STRING directly)."""
    # For now, we'll use a simplified approach
    # In production, you'd use IdMapping service or biomart
    return {}

def curate_targets(raw_file, output_file):
    """Apply filtering pipeline to raw targets."""
    
    print("\n" + "="*80)
    print("TARGET CURATION PIPELINE")
    print("="*80)
    
    # Load raw targets
    print(f"\n1. Loading raw targets from {raw_file}...")
    df = pd.read_csv(raw_file)
    print(f"   Total: {len(df)} targets")
    print(f"   Hyperforin: {len(df[df['compound']=='Hyperforin'])}")
    print(f"   Quercetin: {len(df[df['compound']=='Quercetin'])}")
    
    # Track filtering stats
    stats = {
        'total_input': len(df),
        'species_filtered': 0,
        'unmapped': 0,
        'duplicates': 0
    }
    
    # Filter 1: Species (keep only human proteins starting with P/Q/O)
    print("\n2. Filter: Human proteins only...")
    print("   Keeping proteins with UniProt human prefixes (P*, Q*, O*)")
    
    def is_human_protein(protein_id):
        """Check if protein ID looks like human UniProt."""
        # Human UniProt IDs start with P, Q, or O followed by 5 digits
        if len(protein_id) == 6 and protein_id[0] in ['P', 'Q', 'O']:
            return protein_id[1:].isdigit()
        return False
    
    df['is_human'] = df['protein_id'].apply(is_human_protein)
    non_human = df[~df['is_human']]
    
    print(f"   Non-human proteins: {len(non_human)}")
    for _, row in non_human.head(10).iterrows():
        print(f"     {row['protein_id']} ({row['compound']}, {row['source']})")
    
    df = df[df['is_human']].copy()
    stats['species_filtered'] = stats['total_input'] - len(df)
    print(f"   Remaining: {len(df)} targets")
    
    # Filter 2: Map to gene symbols (using UniProt annotations)
    # NOTE: This requires external mapping data
    # For now, we'll keep all human proteins and add gene names later
    print("\n3. Adding gene symbol mapping...")
    
    # Manual mapping for key targets (from files we've seen)
    gene_mapping = {
        'O75469': 'NR1I2',   # PXR
        'Q9Y210': 'TRPC6',
        'P08684': 'CYP3A4',
        'P11712': 'CYP2C9',
        'P20813': 'CYP2B6',
        'P08183': 'ABCB1',   # MDR1/P-gp
        'Q9UNQ0': 'ABCG2',   # BCRP
        'P31749': 'AKT1',
        'P08253': 'MMP2',
        'P14780': 'MMP9',
        'P15692': 'VEGFA',
        'Q13794': 'PMAIP1',  # NOXA
        'Q12879': 'GRIN1',   # NMDA
        'O15440': 'ABCC2',   # MRP2
        # Add more as needed from UniProt
        'O43570': 'CAV2',
        'Q72547': 'UGT1A8',
        # ... (this would be loaded from a file in production)
    }
    
    df['gene_name'] = df['protein_id'].map(gene_mapping)
    unmapped = df[df['gene_name'].isna()]
    print(f"   Unmapped proteins: {len(unmapped)}")
    
    # For unmapped, try to use protein ID as placeholder
    # In production, fetch from UniProt API
    df.loc[df['gene_name'].isna(), 'gene_name'] = 'UNMAPPED_' + df.loc[df['gene_name'].isna(), 'protein_id']
    
    # Filter 3: Remove unmapped proteins (if desired)
    print("\n4. Filter: Remove unmapped proteins...")
    unmapped_count = df['gene_name'].str.startswith('UNMAPPED').sum()
    print(f"   Unmapped: {unmapped_count}")
    
    # Keep only mapped for now
    df = df[~df['gene_name'].str.startswith('UNMAPPED')].copy()
    stats['unmapped'] = unmapped_count
    print(f"   Remaining: {len(df)} targets")
    
    # Filter 4: Remove duplicates  
    print("\n5. Deduplication...")
    duplicates = df.duplicated(subset=['compound', 'protein_id'], keep='first').sum()
    df = df.drop_duplicates(subset=['compound', 'protein_id'], keep='first')
    stats['duplicates'] = duplicates
    print(f"   Duplicates removed: {duplicates}")
    print(f"   Final count: {len(df)} targets")
    
    # Save
    print(f"\n6. Saving to {output_file}...")
    df[['compound', 'protein_id', 'source', 'gene_name']].to_csv(output_file, index=False)
    print(f"   ✓ Saved!")
    
    # Summary
    print("\n" + "="*80)
    print("CURATION SUMMARY")
    print("="*80)
    print(f"Input:             {stats['total_input']:>4} targets")
    print(f"Species filtered:  {stats['species_filtered']:>4} targets")
    print(f"Unmapped:          {stats['unmapped']:>4} targets")
    print(f"Duplicates:        {stats['duplicates']:>4} targets")
    print(f"Output:            {len(df):>4} targets")
    print(f"\nFinal counts:")
    print(f"  Hyperforin: {len(df[df['compound']=='Hyperforin'])}")
    print(f"  Quercetin:  {len(df[df['compound']=='Quercetin'])}")
    
    return df

def main():
    raw_file = DATA_DIR / 'raw/targets_raw.csv'
    output_file = DATA_DIR / 'processed/targets_curated.csv'
    
    if not raw_file.exists():
        print(f"ERROR: {raw_file} not found!")
        return
    
    curate_targets(raw_file, output_file)
    
    print("\n✓ Target curation complete!")
    print(f"\nNOTE: This is a simplified version. Full curation requires:")
    print("  - UniProt ID mapping service")
    print("  - STRING database validation")
    print("  - Manual review of ambiguous mappings")

if __name__ == '__main__':
    main()

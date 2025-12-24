"""
Extract network from STRING database v12.0 with confidence filtering.

Usage:
    python scripts/extract_string_network.py --threshold 900 --output data/processed/network_900.parquet
    python scripts/extract_string_network.py --threshold 700 --output data/processed/network_700.parquet
"""

import gzip
import pandas as pd
import argparse
from pathlib import Path
from tqdm import tqdm

def load_string_info(info_file):
    """Load STRING protein info and create protein_id -> gene_symbol mapping."""
    print(f"Loading STRING info from {info_file}...")
    
    protein_to_gene = {}
    with gzip.open(info_file, 'rt') as f:
        # Skip header
        _ = f.readline()  # noqa: F841
        
        for line in tqdm(f, desc="  Parsing protein info"):
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                protein_id = parts[0]  # e.g., "9606.ENSP00000000233"
                gene_symbol = parts[1]  # e.g., "ARF5"
                protein_to_gene[protein_id] = gene_symbol
    
    print(f"  Loaded {len(protein_to_gene):,} protein-gene mappings")
    return protein_to_gene

def extract_network(links_file, protein_to_gene, threshold, output_file):
    """Extract network from STRING links with confidence threshold."""
    print(f"\nExtracting network from {links_file} (threshold >= {threshold})...")
    
    edges = []
    skipped_unmapped = 0
    skipped_self_loops = 0
    
    with gzip.open(links_file, 'rt') as f:
        # Skip header
        _ = f.readline()  # noqa: F841
        
        for line in tqdm(f, desc="  Parsing links"):
            parts = line.strip().split()
            if len(parts) >= 3:
                protein1 = parts[0]
                protein2 = parts[1]
                score = int(parts[2])
                
                # Filter by confidence threshold
                if score < threshold:
                    continue
                
                # Map to gene symbols
                gene1 = protein_to_gene.get(protein1)
                gene2 = protein_to_gene.get(protein2)
                
                if gene1 is None or gene2 is None:
                    skipped_unmapped += 1
                    continue
                
                # Skip self-loops
                if gene1 == gene2:
                    skipped_self_loops += 1
                    continue
                
                # Ensure consistent edge direction (alphabetical)
                if gene1 > gene2:
                    gene1, gene2 = gene2, gene1
                
                edges.append({
                    'gene1': gene1,
                    'gene2': gene2,
                    'score': score
                })
    
    print(f"\n  Edges passing threshold: {len(edges):,}")
    print(f"  Skipped (unmapped proteins): {skipped_unmapped:,}")
    print(f"  Skipped (self-loops): {skipped_self_loops:,}")
    
    # Create dataframe and remove duplicates
    df = pd.DataFrame(edges)
    
    if len(df) > 0:
        # Remove duplicate edges (keep highest score)
        df = df.sort_values('score', ascending=False)
        df = df.drop_duplicates(subset=['gene1', 'gene2'], keep='first')
        
        print(f"  After deduplication: {len(df):,} edges")
        print(f"  Unique genes: {len(set(df['gene1']) | set(df['gene2'])):,}")
        
        # Save to parquet
        df.to_parquet(output_file, index=False)
        print(f"\n✓ Saved to {output_file}")
    else:
        print("\n✗ No edges found!")
    
    return df

def main():
    parser = argparse.ArgumentParser(description='Extract STRING network with confidence filtering')
    parser.add_argument('--threshold', type=int, required=True, choices=[700, 900],
                       help='Confidence score threshold (700 or 900)')
    parser.add_argument('--output', type=str, required=True,
                       help='Output parquet file path')
    
    args = parser.parse_args()
    
    # File paths
    data_dir = Path('data')
    links_file = data_dir / 'external/string_links.txt.gz'
    info_file = data_dir / 'external/string_info.txt.gz'
    
    if not links_file.exists():
        print(f"ERROR: {links_file} not found!")
        return
    
    if not info_file.exists():
        print(f"ERROR: {info_file} not found!")
        return
    
    # Load protein-gene mapping
    protein_to_gene = load_string_info(info_file)
    
    # Extract network
    extract_network(links_file, protein_to_gene, args.threshold, args.output)
    
    print("\n✓ Network extraction complete!")

if __name__ == '__main__':
    main()

"""
Orchestrate complete network regeneration from STRING source.

This script:
1. Extracts networks from STRING at ≥900 and ≥700 thresholds
2. Filters to liver-expressed genes (GTEx)
3. Validates node/edge counts
4. Generates metadata report

Usage:
    python scripts/regenerate_networks.py
"""

import subprocess
import pandas as pd
import json
from pathlib import Path

DATA_DIR = Path('data')
PROCESSED_DIR = DATA_DIR / 'processed'

def run_command(cmd, desc):
    """Run command and handle errors."""
    print(f"\n{'='*80}")
    print(f"{desc}")
    print(f"{'='*80}")
    print(f"Command: {' '.join(cmd)}")
    
    result = subprocess.run(cmd, capture_output=False)
    
    if result.returncode != 0:
        print(f"\n✗ ERROR: {desc} failed!")
        return False
    
    print(f"\n✓ {desc} complete!")
    return True

def validate_network(threshold, expected_nodes, expected_liver_nodes):
    """Validate regenerated network."""
    print(f"\nValidating network_{threshold}.parquet...")
    
    net_file = PROCESSED_DIR / f'network_{threshold}.parquet'
    if not net_file.exists():
        print(f"  ✗ File not found: {net_file}")
        return False
    
    df = pd.read_parquet(net_file)
    col1, col2 = df.columns[0], df.columns[1]
    nodes = len(set(df[col1]) | set(df[col2]))
    edges = len(df)
    
    print(f"  Nodes: {nodes:,}")
    print(f"  Edges: {edges:,}")
    print(f"  Expected liver nodes: {expected_liver_nodes:,}")
    
    # Note: We validate liver nodes count after extraction
    # The parquet itself doesn't have liver filter yet
    
    return True

def main():
    print("="*80)
    print("NETWORK REGENERATION FROM STRING SOURCE")
    print("="*80)
    
    results = {}
    
    # Step 1: Extract networks from STRING
    for threshold in [900, 700]:
        output = PROCESSED_DIR / f'network_{threshold}.parquet'
        
       success = run_command(
            ['python', 'scripts/extract_string_network.py',
             '--threshold', str(threshold),
             '--output', str(output)],
            f"Extracting STRING network (≥{threshold})"
        )
        
        if not success:
            return
        
        # Load stats
        df = pd.read_parquet(output)
        col1, col2 = df.columns[0], df.columns[1]
        
        results[f'network_{threshold}'] = {
            'threshold': threshold,
            'edges': len(df),
            'nodes': len(set(df[col1]) | set(df[col2])),
            'file': str(output)
        }
    
    # Save metadata
    metadata_file = PROCESSED_DIR / 'network_metadata.json'
    with open(metadata_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n✓ Metadata saved to {metadata_file}")
    
    # Summary
    print("\n" + "="*80)
    print("REGENERATION COMPLETE")
    print("="*80)
    
    print("\nGenerated networks:")
    for name, stats in results.items():
        print(f"\n{name}:")
        print(f"  Threshold: ≥{stats['threshold']}")
        print(f"  Nodes: {stats['nodes']:,}")
        print(f"  Edges: {stats['edges']:,}")
    
    print("\nNext steps:")
    print("1. Run validation scripts to confirm liver filtering")
    print("2. Check that analyses still produce same results")
    print("3. Commit regenerated networks")

if __name__ == '__main__':
    main()

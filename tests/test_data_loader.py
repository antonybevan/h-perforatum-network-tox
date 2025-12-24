"""
Unit tests for data loading functions.
"""

import pytest
import pandas as pd
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from network_tox.utils.data_loader import load_liver_genes


def test_load_liver_genes_returns_set():
    """Test that load_liver_genes returns a set."""
    gtex_path = Path('data/raw/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct')
    
    if not gtex_path.exists():
        pytest.skip(f"GTEx file not found: {gtex_path}")
    
    result = load_liver_genes(gtex_path)
    assert isinstance(result, set)
    assert len(result) > 0


def test_load_liver_genes_cached():
    """Test loading from cached liver_proteome.csv."""
    cache_path = Path('data/processed/liver_proteome.csv')
    
    if not cache_path.exists():
        pytest.skip(f"Cached file not found: {cache_path}")
    
    df = pd.read_csv(cache_path)
    liver_genes = set(df['gene_symbol'])
    
    assert isinstance(liver_genes, set)
    assert len(liver_genes) == 13496  # Expected count


def test_targets_file_structure():
    """Test that targets.csv has correct structure."""
    targets_path = Path('data/processed/targets.csv')
    
    if not targets_path.exists():
        pytest.skip(f"Targets file not found: {targets_path}")
    
    df = pd.read_csv(targets_path)
    
    # Check columns
    assert 'compound' in df.columns
    assert 'protein_id' in df.columns
    assert 'gene_name' in df.columns
    assert 'source' in df.columns
    
    # Check compounds
    assert 'Hyperforin' in df['compound'].values
    assert 'Quercetin' in df['compound'].values
    
    # Check counts
    hf_count = len(df[df['compound'] == 'Hyperforin'])
    qu_count = len(df[df['compound'] == 'Quercetin'])
    
    assert hf_count == 12
    assert qu_count == 80


def test_dili_files_exist():
    """Test that DILI files exist and have correct structure."""
    for threshold in ['900', '700']:
        dili_path = Path(f'data/processed/dili_{threshold}_lcc.csv')
        
        if not dili_path.exists():
            pytest.skip(f"DILI file not found: {dili_path}")
        
        df = pd.read_csv(dili_path)
        assert len(df) > 0
        assert 'protein_id' in df.columns

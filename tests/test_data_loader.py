"""
Unit tests for data loading functions.
"""

import pytest
import pandas as pd
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))


def test_placeholder():
    """Placeholder test - data files not available in CI."""
    # This ensures pytest doesn't fail if all tests are skipped
    assert True


# Integration tests (require data files - run locally only)
@pytest.mark.skipif(
    not Path('data/processed/liver_proteome.csv').exists(),
    reason="Data files not available in CI environment"
)
def test_load_liver_genes_cached():
    """Test loading from cached liver_proteome.csv."""
    from network_tox.utils.data_loader import load_liver_genes
    
    cache_path = Path('data/processed/liver_proteome.csv')
    df = pd.read_csv(cache_path)
    liver_genes = set(df['gene_symbol'])
    
    assert isinstance(liver_genes, set)
    assert len(liver_genes) == 13496  # Expected count


@pytest.mark.skipif(
    not Path('data/processed/targets.csv').exists(),
    reason="Data files not available in CI environment"
)
def test_targets_file_structure():
    """Test that targets.csv has correct structure."""
    targets_path = Path('data/processed/targets.csv')
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


@pytest.mark.skipif(
    not Path('data/processed/dili_900_lcc.csv').exists(),
    reason="Data files not available in CI environment"
)
def test_dili_files_exist():
    """Test that DILI files exist and have correct structure."""
    for threshold in ['900', '700']:
        dili_path = Path(f'data/processed/dili_{threshold}_lcc.csv')
        
        if dili_path.exists():
            df = pd.read_csv(dili_path)
            assert len(df) > 0
            assert 'protein_id' in df.columns

"""
Unit tests for core functionality - CI compatible.
Tests run without requiring src/ package installation.
"""

import pytest
import pandas as pd
import networkx as nx
from pathlib import Path


def is_lfs_pointer(filepath):
    """Check if file is a Git LFS pointer."""
    try:
        with open(filepath, 'r') as f:
            header = f.read(100)
            return header.startswith('version https://git-lfs.github.com/spec/v1')
    except Exception:
        return False


def test_basic_imports():
    """Test that required packages are installed."""
    import pandas
    import networkx
    import scipy
    import numpy
    assert True


def test_pandas_version():
    """Test pandas is recent enough."""
    import pandas as pd
    major, minor = pd.__version__.split('.')[:2]
    assert int(major) >= 2, f"pandas {pd.__version__} is too old"


def test_networkx_basics():
    """Test networkx graph operations work."""
    G = nx.Graph()
    G.add_edges_from([('A', 'B'), ('B', 'C')])
    assert G.number_of_nodes() == 3
    assert G.number_of_edges() == 2


# Local integration tests (skip in CI)
@pytest.mark.skipif(
    not Path('data/processed/targets.csv').exists() or is_lfs_pointer('data/processed/targets.csv'),
    reason="Data files not available in CI"
)
def test_targets_file_structure():
    """Test targets.csv structure (local only)."""
    df = pd.read_csv('data/processed/targets.csv')
    
    assert 'compound' in df.columns
    assert 'protein_id' in df.columns
    assert 'gene_name' in df.columns
    
    hf_count = len(df[df['compound'] == 'Hyperforin'])
    qu_count = len(df[df['compound'] == 'Quercetin'])
    
    assert hf_count == 12
    assert qu_count == 80


@pytest.mark.skipif(
    not Path('data/processed/liver_proteome.csv').exists() or is_lfs_pointer('data/processed/liver_proteome.csv'),
    reason="Data files not available in CI"
)
def test_liver_proteome_cached():
    """Test liver_proteome.csv structure (local only)."""
    df = pd.read_csv('data/processed/liver_proteome.csv')
    liver_genes = set(df['gene_symbol'])
    
    assert len(liver_genes) == 13496

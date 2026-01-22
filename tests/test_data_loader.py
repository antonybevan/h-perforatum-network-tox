"""Tests for core functionality."""

import pytest
import pandas as pd
import networkx as nx
from pathlib import Path


def is_lfs_pointer(filepath: Path) -> bool:
    """Check if file is a Git LFS pointer (not actual data).
    
    LFS pointers are small text files that look like:
        version https://git-lfs.github.com/spec/v1
        oid sha256:abc123...
        size 12345
    
    Returns True if file is a pointer, False if it's real data.
    """
    if not filepath.exists():
        return False
    
    try:
        with open(filepath, 'rb') as f:
            header = f.read(100)
            return b'version https://git-lfs' in header
    except Exception:
        return False


def data_file_ready(filepath: str) -> bool:
    """Check if a data file exists and contains real data (not LFS pointer)."""
    p = Path(filepath)
    if not p.exists():
        return False
    if is_lfs_pointer(p):
        return False
    return True


def skip_reason(filepath: str) -> str:
    """Generate informative skip reason for missing/pointer files."""
    p = Path(filepath)
    if not p.exists():
        return f"File not found: {filepath}"
    if is_lfs_pointer(p):
        return f"LFS pointer detected (run 'git lfs pull'): {filepath}"
    return ""


# Basic Package Tests

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


# Data Integration Tests

TARGETS_FILE = 'data/processed/targets_lcc.csv'  # LCC-filtered targets
LIVER_FILE = 'data/processed/liver_proteome.csv'
NETWORK_FILE = 'data/processed/network_900.parquet'


@pytest.mark.skipif(
    not data_file_ready(TARGETS_FILE),
    reason=skip_reason(TARGETS_FILE) or "Data not available"
)
def test_targets_file_structure():
    """Test targets_lcc.csv structure (LCC-filtered targets)."""
    df = pd.read_csv(TARGETS_FILE)
    
    assert 'compound' in df.columns
    assert 'gene_symbol' in df.columns
    assert 'in_lcc_900' in df.columns
    
    hf_count = len(df[df['compound'] == 'Hyperforin'])
    qu_count = len(df[df['compound'] == 'Quercetin'])
    
    assert hf_count == 10, f"Expected 10 Hyperforin targets (LCC-filtered), got {hf_count}"
    assert qu_count == 62, f"Expected 62 Quercetin targets (LCC-filtered), got {qu_count}"


@pytest.mark.skipif(
    not data_file_ready(LIVER_FILE),
    reason=skip_reason(LIVER_FILE) or "Data not available"
)
def test_liver_proteome_structure():
    """Test liver_proteome.csv structure."""
    df = pd.read_csv(LIVER_FILE)
    liver_genes = set(df['gene_symbol'])
    
    assert len(liver_genes) == 13496, f"Expected 13496 liver genes, got {len(liver_genes)}"


@pytest.mark.skipif(
    not data_file_ready(NETWORK_FILE),
    reason=skip_reason(NETWORK_FILE) or "Data not available"
)
def test_network_file_structure():
    """Test network_900.parquet structure."""
    df = pd.read_parquet(NETWORK_FILE)
    
    assert 'protein1' in df.columns
    assert 'protein2' in df.columns
    assert 'weight' in df.columns
    
    # Verify network size is reasonable
    assert len(df) > 50000, f"Network too small: {len(df)} edges"


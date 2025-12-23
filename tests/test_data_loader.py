"""Tests for data loader."""

import pandas as pd
from network_tox.utils import data_loader

def test_load_liver_genes(tmp_path):
    """Test loading and filtering liver genes."""
    # Create mock GTEx file
    d = {
        'Name': ['G1', 'G2', 'G3'],
        'Description': ['GENE1', 'GENE2', 'GENE3'],
        'Liver': [0.5, 1.5, 2.0],
        'Other': [10, 10, 10]
    }

    # Create file content with skip rows
    lines = ["#1.2", "100 100"]

    # Create dataframe and convert to tab-separated string
    df = pd.DataFrame(d)
    csv_content = df.to_csv(sep='\t', index=False)

    # Combine
    full_content = "\n".join(lines) + "\n" + csv_content

    p = tmp_path / "test_gtex.gct"
    p.write_text(full_content)

    genes = data_loader.load_liver_genes(str(p))

    assert 'GENE1' not in genes
    assert 'GENE2' in genes
    assert 'GENE3' in genes
    assert len(genes) == 2

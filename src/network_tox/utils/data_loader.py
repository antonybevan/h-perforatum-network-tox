"""Data loading utilities."""

import pandas as pd

def load_liver_genes(gtex_path):
    """
    Load liver genes from GTEx data with Median TPM > 1.

    Args:
        gtex_path: Path to GTEx GCT file

    Returns:
        Set of gene symbols (Description column)
    """
    # GCT format has 2 header lines to skip
    df = pd.read_csv(gtex_path, sep='\t', skiprows=2)

    # Filter for Liver > 1
    liver_genes = df[df['Liver'] > 1]['Description']

    return set(liver_genes)

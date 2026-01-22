import pandas as pd

def load_liver_genes(gtex_path):
    """
    Reads the GTEx GCT file, filters for Liver columns, and returns a set of gene symbols
    where the median TPM is greater than 1.

    Args:
        gtex_path (str): Path to the GTEx GCT file.

    Returns:
        set: A set of gene symbols (strings).
    """
    # Load the dataset, skipping the first 2 header lines
    df = pd.read_csv(gtex_path, sep='\t', skiprows=2)

    # Filter for columns containing 'Liver' and the gene symbol column 'Description'
    liver_cols = [col for col in df.columns if 'Liver' in col]

    if not liver_cols:
        raise ValueError("No columns containing 'Liver' found in the dataset.")

    # We also need the 'Description' column which contains gene symbols
    cols_to_keep = ['Description'] + liver_cols
    df_liver = df[cols_to_keep]

    # Calculate median TPM across liver columns (though usually there is only one "Liver" column in this file,
    # or if there are multiple, we take the median across them per row)
    # The instructions say "calculate the median TPM". Since we filtered for liver columns,
    # if there is more than one, we take the median. If there is one, it's just the value.
    # We convert to numeric just in case
    for col in liver_cols:
        df_liver.loc[:, col] = pd.to_numeric(df_liver[col], errors='coerce')

    # Calculate median across the liver columns
    # axis=1 operates on rows.
    df_liver.loc[:, 'median_tpm'] = df_liver[liver_cols].median(axis=1)

    # Filter where Median TPM > 1
    expressed_genes = df_liver[df_liver['median_tpm'] > 1]['Description']

    # Return as a set of unique gene symbols
    return set(expressed_genes.unique())

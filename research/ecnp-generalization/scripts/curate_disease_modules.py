"""
Disease Module Curation for ECNP Generalization Study

Downloads disease gene sets from DisGeNET or uses curated lists,
then filters to network LCC for consistency with ECNP analysis.

Disease Modules:
1. DILI - Drug-Induced Liver Injury (already have)
2. Cancer - Pan-cancer driver genes (COSMIC/IntOGen)
3. Type 2 Diabetes - T2D-associated genes
4. Alzheimer's Disease - AD-associated genes
5. Cardiovascular Disease - CVD-associated genes

Output: research/ecnp-generalization/data/{disease}_lcc.csv
"""

import pandas as pd
import numpy as np
from pathlib import Path

# Paths
PROJECT_ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = PROJECT_ROOT / "data" / "processed"
OUTPUT_DIR = PROJECT_ROOT / "research" / "ecnp-generalization" / "data"

# Load network node list for filtering
def load_network_nodes():
    """Get genes in LCC network."""
    edges = pd.read_parquet(DATA_DIR / "network_900_liver_lcc.parquet")
    nodes = set(edges['gene1'].tolist() + edges['gene2'].tolist())
    print(f"Network LCC: {len(nodes)} nodes")
    return nodes


def fetch_disgenet_genes(disease_name: str, disease_id: str, min_score: float = 0.3):
    """
    Fetch disease genes from DisGeNET.
    
    Note: DisGeNET requires API key for programmatic access.
    For this script, we'll use curated gene lists instead.
    """
    # DisGeNET API would require key - using curated lists below
    pass


def get_curated_disease_genes():
    """
    Curated disease gene sets from literature/databases.
    
    Sources:
    - DILI: Already in data/processed/dili_900_lcc.csv
    - Cancer: COSMIC Cancer Gene Census
    - T2D: GWAS Catalog + literature
    - Alzheimer's: AlzGene + GWAS
    - CVD: CARDIoGRAM + literature
    """
    
    # Cancer driver genes (COSMIC Cancer Gene Census subset)
    cancer_genes = [
        "TP53", "KRAS", "PIK3CA", "PTEN", "APC", "EGFR", "BRAF", "RB1",
        "CDKN2A", "MYC", "BRCA1", "BRCA2", "ATM", "VHL", "NF1", "NRAS",
        "CTNNB1", "SMAD4", "FBXW7", "NOTCH1", "IDH1", "JAK2", "KIT",
        "MET", "ALK", "RET", "FGFR1", "FGFR2", "FGFR3", "ERBB2", "CDK4",
        "CDK6", "CCND1", "CCNE1", "MDM2", "BCL2", "BCL6", "MLL", "EZH2",
        "DNMT3A", "TET2", "ASXL1", "SF3B1", "U2AF1", "SRSF2", "NPM1",
        "FLT3", "RUNX1", "CEBPA", "WT1", "GATA2", "KMT2A", "NSD1"
    ]
    
    # Type 2 Diabetes genes (GWAS + functional validation)
    t2d_genes = [
        "TCF7L2", "PPARG", "KCNJ11", "ABCC8", "SLC30A8", "HNF1A", "HNF4A",
        "GCK", "GLIS3", "MTNR1B", "IRS1", "IRS2", "CDKN2A", "CDKN2B",
        "IGF2BP2", "FTO", "HHEX", "CDKAL1", "WFS1", "JAZF1", "CDC123",
        "CAMK1D", "TSPAN8", "LGR5", "THADA", "ADAMTS9", "NOTCH2", "PROX1",
        "BCL11A", "ZBED3", "KLF14", "KCNQ1", "CENTD2", "AP3S2", "HNF1B",
        "ARAP1", "ZFAND6", "PRC1", "DUSP9", "HMGA2", "GRB14", "ANKRD55"
    ]
    
    # Alzheimer's Disease genes (AlzGene + GWAS)
    ad_genes = [
        "APP", "PSEN1", "PSEN2", "APOE", "TREM2", "CLU", "BIN1", "CR1",
        "PICALM", "MS4A6A", "CD33", "ABCA7", "EPHA1", "CD2AP", "SORL1",
        "HLA-DRB5", "PTK2B", "SLC24A4", "DSG2", "INPP5D", "MEF2C",
        "NME8", "ZCWPW1", "CELF1", "FERMT2", "CASS4", "ECHDC3", "ACE",
        "ADAM10", "APH1B", "NCSTN", "MAPT", "GRN", "TARDBP", "FUS",
        "CHMP2B", "VCP", "SQSTM1", "TBK1", "OPTN", "UBQLN2", "C9orf72"
    ]
    
    # Cardiovascular Disease genes (CARDIoGRAM + functional)
    cvd_genes = [
        "LDLR", "APOB", "PCSK9", "APOE", "APOA1", "CETP", "LIPC", "LPL",
        "ABCA1", "ABCG5", "ABCG8", "NPC1L1", "HMGCR", "ANGPTL3", "ANGPTL4",
        "LPA", "SORT1", "MRAS", "PHACTR1", "WDR12", "ADAMTS7", "ABO",
        "TCF21", "ZC3HC1", "EDNRA", "GUCY1A3", "MIA3", "SMAD3", "COL4A1",
        "COL4A2", "FN1", "VAMP5", "VAMP8", "GGCX", "PROCR", "F2", "F5",
        "F7", "F10", "SERPINC1", "PROS1", "PROC", "PLG", "FGA", "FGB"
    ]
    
    return {
        'cancer': cancer_genes,
        't2d': t2d_genes,
        'alzheimer': ad_genes,
        'cvd': cvd_genes
    }


def create_disease_module(disease_name: str, gene_list: list, network_nodes: set):
    """Create filtered disease module CSV."""
    # Filter to network
    in_network = [g for g in gene_list if g in network_nodes]
    
    df = pd.DataFrame({
        'gene_name': in_network,
        'disease': disease_name
    })
    
    print(f"{disease_name}: {len(gene_list)} genes, {len(in_network)} in network")
    
    return df


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    print("Loading network nodes...")
    network_nodes = load_network_nodes()
    
    print("\nCurating disease modules...")
    disease_genes = get_curated_disease_genes()
    
    # Copy DILI from existing
    dili = pd.read_csv(DATA_DIR / "dili_900_lcc.csv")
    dili.to_csv(OUTPUT_DIR / "dili_lcc.csv", index=False)
    print(f"DILI: {len(dili)} genes (copied from existing)")
    
    # Create other disease modules
    for disease, genes in disease_genes.items():
        df = create_disease_module(disease, genes, network_nodes)
        df.to_csv(OUTPUT_DIR / f"{disease}_lcc.csv", index=False)
    
    print("\n" + "="*50)
    print("Disease modules saved to:")
    for f in OUTPUT_DIR.glob("*.csv"):
        n = len(pd.read_csv(f))
        print(f"  {f.name}: {n} genes")


if __name__ == "__main__":
    main()

"""
DILI Gene Evidence Audit
=========================

For each of the 82 DILI genes:
1. Categorize evidence type
2. Flag and EXCLUDE genes with drug-specific evidence
3. Keep only disease-mechanism genes

Evidence Categories:
- GWAS: From genome-wide association studies (safe)
- CURATED: From databases like UniProt, ClinVar (safe if not drug-linked)
- MECHANISM: Liver biology genes (safe)
- DRUG-LINKED: Evidence derived from specific drug studies (EXCLUDE)
"""
import pandas as pd
from pathlib import Path

ROOT = Path(r'v:\new\h-perforatum-network-tox')

# Load DILI genes
genes_df = pd.read_csv(ROOT / 'data' / 'processed' / 'dili_900_lcc.csv')
print(f"Total DILI genes in LCC: {len(genes_df)}")

# =============================================================================
# GENE CLASSIFICATION BASED ON KNOWN BIOLOGY
# =============================================================================

# These genes are GWAS-validated or have known liver biology
# NOT derived from specific drug associations
GWAS_VALIDATED = {
    # Pharmacogenomics genes (GWAS validated for DILI susceptibility)
    'GSTM1', 'GSTT1', 'NAT2', 'CYP2E1', 'CYP2C9', 'CYP2C19',
    'MTHFR', 'POLG', 'GSTP1', 'GSTM2', 'UGT1A9',
    
    # HLA genes (immune susceptibility - validated by GWAS)
    'HLA-A', 'HLA-B', 'HLA-DQB1', 'HLA-DRB1',
}

# Liver biology/mechanism genes (not drug-specific)
LIVER_MECHANISM = {
    # Transporters
    'ABCB1',  # P-glycoprotein
    
    # Nuclear receptors
    'NR1I2', 'NR1I3', 'NR1H3', 'NR1H4', 'PPARA', 'AHR', 'ARNT',
    
    # Liver enzymes
    'ALB', 'GPT', 'FMO3', 'PON1', 'CYP2A6',
    
    # Oxidative stress
    'CAT', 'SOD1', 'SOD3', 'HMOX1', 'NFE2L2', 'GCLC',
    
    # Apoptosis/stress response
    'BAX', 'GADD45A', 'MAP3K5', 'ATG5',
    
    # Liver proteins/biomarkers
    'KRT18', 'ALB', 'TTR', 'TF', 'HPX', 'ALDOB', 'AMBP', 'ARG1',
    'GC', 'APOA1', 'APOE', 'APOH', 'FGA',
    
    # Inflammation (general, not drug-specific)
    'CCL2', 'CXCL1', 'CXCL10', 'HMGB1',
    
    # Other liver biology
    'CLU', 'LCN2', 'PLAT', 'PLG', 'RBP1', 'SPP1', 'C3',
    'CTNNB1', 'IGF1', 'MMP2', 'GSN', 'ENO1', 'TALDO1',
    
    # Growth factors
    'FLT1',
    
    # Cell adhesion
    'LGALS3',
    
    # Lipid metabolism
    'DGAT2',
    
    # Miscellaneous liver
    'BTD', 'MED1', 'PNP', 'PRKDC', 'PTGS2', 'SLPI', 'TBXA2R', 'TCTN1', 'HPD',
    'SNX18', 'COL3A1',
}

# POTENTIALLY DRUG-LINKED (need careful review)
# These are often identified from drug toxicity studies
DRUG_LINKED_RISK = {
    # Cytokines often measured in drug toxicity studies
    'IL1B', 'IL1A', 'IL4', 'IL6', 'IL11', 'IL17A', 'IL18', 'IL22',
    'TNF', 'IFNG', 'IFNA2',
    'IL1R2',
    
    # miRNAs often identified in drug toxicity biomarker studies
    'MIR10A', 'MIR10B', 'MIR122', 'MIR132', 'MIR141', 'MIR149',
    'MIR181C', 'MIR19A', 'MIR200C', 'MIR217', 'MIR218-1',
    'MIR29B2', 'MIR30A', 'MIR337', 'MIR34C', 'MIR367',
    'MIR410', 'MIR503', 'MIR592', 'MIR744', 'MIR764',
}

# =============================================================================
# CLASSIFY EACH GENE
# =============================================================================

def classify_gene(gene_name):
    """Classify gene by evidence type."""
    if gene_name in GWAS_VALIDATED:
        return 'GWAS', 'SAFE', 'GWAS-validated pharmacogenomics'
    elif gene_name in LIVER_MECHANISM:
        return 'MECHANISM', 'SAFE', 'Liver biology/mechanism'
    elif gene_name in DRUG_LINKED_RISK:
        return 'DRUG_LINKED', 'EXCLUDE', 'May be derived from drug studies'
    else:
        return 'UNKNOWN', 'REVIEW', 'Needs manual review'

# Classify all genes
results = []
for _, row in genes_df.iterrows():
    gene = row['gene_name']
    score = row['score']
    
    category, action, reason = classify_gene(gene)
    
    results.append({
        'gene_name': gene,
        'disgenet_score': score,
        'category': category,
        'action': action,
        'reason': reason
    })

audit_df = pd.DataFrame(results)

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "="*60)
print("DILI GENE EVIDENCE AUDIT")
print("="*60)

print("\nGene classification summary:")
print(audit_df['category'].value_counts())

print("\nAction summary:")
print(audit_df['action'].value_counts())

# SAFE genes
safe_genes = audit_df[audit_df['action'] == 'SAFE']
print(f"\nSAFE genes to keep: {len(safe_genes)}")

# EXCLUDE genes
exclude_genes = audit_df[audit_df['action'] == 'EXCLUDE']
print(f"EXCLUDE (drug-linked): {len(exclude_genes)}")
if len(exclude_genes) > 0:
    print("Genes to exclude:")
    for _, row in exclude_genes.iterrows():
        print(f"  - {row['gene_name']}: {row['reason']}")

# REVIEW genes
review_genes = audit_df[audit_df['action'] == 'REVIEW']
print(f"\nNEED REVIEW: {len(review_genes)}")
if len(review_genes) > 0:
    print("Genes needing review:")
    for _, row in review_genes.iterrows():
        print(f"  - {row['gene_name']}")

# =============================================================================
# SAVE AUDIT RESULTS
# =============================================================================

output = ROOT / 'research' / 'ecnp-generalization' / 'results' / 'dili_gene_audit.csv'
audit_df.to_csv(output, index=False)
print(f"\nSaved audit: {output}")

# =============================================================================
# CREATE CLEAN GENE LIST (EXCLUDING DRUG-LINKED)
# =============================================================================

clean_genes = audit_df[audit_df['action'].isin(['SAFE', 'REVIEW'])].copy()
print(f"\nClean gene list: {len(clean_genes)} genes (excluding {len(exclude_genes)} drug-linked)")

clean_output = ROOT / 'data' / 'processed' / 'dili_genes_clean.csv'

# Merge with original data
clean_genes_df = genes_df[genes_df['gene_name'].isin(clean_genes['gene_name'])].copy()
clean_genes_df.to_csv(clean_output, index=False)
print(f"Saved clean list: {clean_output}")

print("\n" + "="*60)
print("RECOMMENDATION")
print("="*60)
print(f"""
ORIGINAL: {len(genes_df)} genes
EXCLUDED: {len(exclude_genes)} genes (drug-linked evidence)
CLEAN:    {len(clean_genes_df)} genes

ACTION: Retrain ECNP with clean gene list to eliminate potential circularity.

Excluded genes are primarily:
  - Cytokines (IL1B, IL6, TNF) - often measured as drug toxicity markers
  - miRNAs (MIR122, etc.) - biomarkers identified in drug toxicity studies
  
These may be discovered BECAUSE of drug hepatotoxicity studies,
creating potential circular logic with DILIrank.
""")

"""
DrugBank XML Parser - Extract Drug-Target Pairs
================================================

Parses DrugBank XML to extract:
- Drug name and DrugBank ID
- Target proteins (gene symbols)
- Known actions (inhibitor, inducer, etc.)

Filters to targets present in Liver LCC network.
"""
import xml.etree.ElementTree as ET
from pathlib import Path
import pandas as pd
import sys

# Paths
RESEARCH_ROOT = Path(r'v:\new\h-perforatum-network-tox\research')
DRUGBANK_XML = RESEARCH_ROOT / 'ecnp-generalization' / 'data' / 'raw' / 'full database.xml'
OUTPUT_DIR = RESEARCH_ROOT / 'ecnp-generalization' / 'data' / 'curated'
NETWORK_DATA = RESEARCH_ROOT / 'ecnp-closed-form' / 'data' / 'node_list_900.csv'

# Check if alternative file name
if not DRUGBANK_XML.exists():
    alt_path = RESEARCH_ROOT / 'ecnp-generalization' / 'data' / 'raw' / 'drugbank_all_full_database.xml'
    if alt_path.exists():
        DRUGBANK_XML = alt_path
    else:
        # Try to find any XML file
        raw_dir = RESEARCH_ROOT / 'ecnp-generalization' / 'data' / 'raw'
        xml_files = list(raw_dir.glob('*.xml'))
        if xml_files:
            DRUGBANK_XML = xml_files[0]
        else:
            print(f"ERROR: Cannot find DrugBank XML file in {raw_dir}")
            print(f"Files in directory: {list(raw_dir.iterdir())}")
            sys.exit(1)

print(f"Using DrugBank XML: {DRUGBANK_XML}")
print(f"File size: {DRUGBANK_XML.stat().st_size / 1e9:.2f} GB")

# Load network nodes (Liver LCC)
print("\nLoading Liver LCC network nodes...")
nodes_df = pd.read_csv(NETWORK_DATA)
liver_genes = set(nodes_df['gene'].tolist())
print(f"Liver LCC genes: {len(liver_genes)}")

# Parse DrugBank XML
print("\nParsing DrugBank XML (this may take a few minutes)...")

# DrugBank namespace
ns = {'db': 'http://www.drugbank.ca'}

# Use iterparse for memory efficiency
drug_targets = []
drug_count = 0

context = ET.iterparse(str(DRUGBANK_XML), events=('end',))

for event, elem in context:
    if elem.tag == '{http://www.drugbank.ca}drug':
        drug_count += 1
        
        # Get drug info
        drugbank_id = None
        drug_name = None
        drug_type = elem.get('type', 'unknown')
        
        for child in elem:
            if child.tag == '{http://www.drugbank.ca}drugbank-id' and child.get('primary') == 'true':
                drugbank_id = child.text
            elif child.tag == '{http://www.drugbank.ca}name':
                drug_name = child.text
        
        # Get targets
        targets_elem = elem.find('db:targets', ns)
        if targets_elem is not None:
            for target in targets_elem.findall('db:target', ns):
                # Get gene name
                polypeptide = target.find('db:polypeptide', ns)
                if polypeptide is not None:
                    gene_name = polypeptide.find('db:gene-name', ns)
                    if gene_name is not None and gene_name.text:
                        gene = gene_name.text.strip().upper()
                        
                        # Get action
                        actions_elem = target.find('db:actions', ns)
                        actions = []
                        if actions_elem is not None:
                            for action in actions_elem.findall('db:action', ns):
                                if action.text:
                                    actions.append(action.text)
                        
                        drug_targets.append({
                            'drugbank_id': drugbank_id,
                            'drug_name': drug_name,
                            'drug_type': drug_type,
                            'gene_symbol': gene,
                            'actions': ';'.join(actions) if actions else '',
                            'in_liver_lcc': gene in liver_genes
                        })
        
        # Clear element to free memory
        elem.clear()
        
        if drug_count % 1000 == 0:
            print(f"  Processed {drug_count} drugs...")

print(f"\nTotal drugs processed: {drug_count}")

# Create dataframe
df = pd.DataFrame(drug_targets)
print(f"Total drug-target pairs: {len(df)}")
print(f"Pairs with targets in Liver LCC: {df['in_liver_lcc'].sum()}")

# Filter to Liver LCC targets
df_liver = df[df['in_liver_lcc']].copy()
print(f"\nFiltered to Liver LCC: {len(df_liver)} pairs")

# Summary per drug
drug_summary = df_liver.groupby(['drugbank_id', 'drug_name']).agg({
    'gene_symbol': list,
    'in_liver_lcc': 'count'
}).reset_index()
drug_summary.columns = ['drugbank_id', 'drug_name', 'targets', 'n_targets']
drug_summary = drug_summary[drug_summary['n_targets'] >= 3]  # Minimum 3 targets
print(f"Drugs with ≥3 targets in Liver LCC: {len(drug_summary)}")

# Save outputs
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Full pairs
df_liver.to_csv(OUTPUT_DIR / 'drugbank_target_pairs.csv', index=False)
print(f"\nSaved: {OUTPUT_DIR / 'drugbank_target_pairs.csv'}")

# Drug summary with target lists
drug_summary.to_csv(OUTPUT_DIR / 'drugbank_compounds.csv', index=False)
print(f"Saved: {OUTPUT_DIR / 'drugbank_compounds.csv'}")

# Quick stats
print(f"\n{'='*50}")
print("DRUGBANK CURATION SUMMARY")
print('='*50)
print(f"Total drugs in DrugBank: {drug_count}")
print(f"Total drug-target pairs: {len(df)}")
print(f"Pairs in Liver LCC: {len(df_liver)}")
print(f"Drugs with ≥3 Liver LCC targets: {len(drug_summary)}")
print(f"Target distribution:")
print(drug_summary['n_targets'].describe())

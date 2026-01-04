# Phase 1: Data Infrastructure Scripts

This folder contains scripts for curating compound-target and toxicity label databases.

## Scripts

- `download_drugbank.py` - Download and parse DrugBank XML
- `download_dilirank.py` - Download DILIrank toxicity labels
- `curate_targets.py` - Filter and map targets to Liver LCC
- `validate_coverage.py` - Check compound-target coverage

## Data Sources

| Source | URL |
|--------|-----|
| DrugBank | https://go.drugbank.com/releases/latest |
| DILIrank | https://www.fda.gov/science-research/liver-toxicity-knowledge-base-ltkb/drug-induced-liver-injury-rank-dilirank-dataset |
| ChEMBL | https://www.ebi.ac.uk/chembl/ |

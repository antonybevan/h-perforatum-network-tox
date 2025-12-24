# Network Generation Documentation

## Overview

This document describes the complete pipeline for generating protein-protein interaction networks from STRING v12.0 database with tissue-specific (liver) filtering.

---

## Data Sources

### Primary Source: STRING v12.0
- **Database:** STRING (Search Tool for the Retrieval of Interacting Genes/Proteins)
- **Version:** 12.0
- **Organism:** 9606 (Homo sapiens)
- **Files:**
  ```
  data/external/string_links.txt.gz    # Protein-protein interactions with confidence scores
  data/external/string_info.txt.gz     # Protein ID to gene symbol mapping
  ```

### Tissue Specificity: GTEx v8
- **Database:** Genotype-Tissue Expression (GTEx)
- **Version:** v8 (2017-06-05)
- **Tissue:** Liver
- **Filter:** Median TPM > 1 across liver samples
- **Files:**
  ```
  data/raw/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct
  data/processed/liver_proteome.csv    # 13,496 liver-expressed genes
  ```

---

## Pipeline Steps

### 1. STRING Extraction (`extract_string_network.py`)

**Purpose:** Extract protein-protein interactions from STRING at specific confidence thresholds.

**Process:**
1. **Load protein→gene mapping:**
   - Parse `string_info.txt.gz`
   - Create mapping: `9606.ENSP00000000233` → `ARF5`

2. **Filter interactions by confidence:**
   - Parse `string_links.txt.gz`
   - Keep only edges where `combined_score ≥ threshold`
   - Thresholds: 900 (high confidence) or 700 (medium-high confidence)

3. **Map to gene symbols:**
   - Convert protein IDs to human-readable gene symbols
   - Skip unmapped proteins (report count)

4. **Clean edges:**
   - Remove self-loops (gene1 == gene2)
   - Standardize edge direction (alphabetical order)
   - Deduplicate (keep highest score)

5. **Save as parquet:**
   - Columns: `gene1`, `gene2`, `score`
   - Format: Apache Parquet (efficient storage)

**Usage:**
```bash
python scripts/extract_string_network.py --threshold 900 --output data/processed/network_900.parquet
python scripts/extract_string_network.py --threshold 700 --output data/processed/network_700.parquet
```

**Output:**
- `network_900.parquet`: ~11,693 nodes, ~100,383 edges
- `network_700.parquet`: ~15,882 nodes, ~236,712 edges

---

### 2. Liver Filtering (`filter_liver_network.py`)

**Purpose:** Filter generic PPI network to liver-expressed genes and extract Largest Connected Component (LCC).

**Process:**
1. **Load network:**
   - Read parquet file from Step 1

2. **Load liver genes:**
   - Read `liver_proteome.csv` (13,496 genes)
   - Genes with median liver TPM > 1 in GTEx

3. **Filter nodes:**
   - Keep only genes in liver proteome
   - Remove genes not expressed in liver

4. **Extract LCC:**
   - Identify connected components
   - Keep only the largest connected component
   - Removes isolated nodes and small disconnected clusters

5. **Save filtered network:**
   - Update parquet with liver-filtered edges

**Usage:**
```bash
python scripts/filter_liver_network.py --input network_900.parquet --output network_900_liver.parquet
```

**Output:**
- `network_900.parquet` (after filter): ~7,677 nodes, ~XX edges
- `network_700.parquet` (after filter): ~9,773 nodes, ~XX edges

---

### 3. Pipeline Orchestration (`regenerate_networks.py`)

**Purpose:** Automate complete regeneration of both networks.

**Process:**
1. Extract STRING network at ≥900 threshold
2. Extract STRING network at ≥700 threshold
3. Generate metadata report
4. Validate node counts

**Usage:**
```bash
python scripts/regenerate_networks.py
```

**Output:**
- `data/processed/network_900.parquet`
- `data/processed/network_700.parquet`
- `data/processed/network_metadata.json`

---

## Network Statistics

### Network ≥900 (High Confidence)
| Metric | Value |
|--------|-------|
| STRING threshold | ≥900 |
| Total nodes | 11,693 |
| Total edges | 100,383 |
| Liver nodes (GTEx) | 7,677 |
| Liver edges | ~XX,XXX |

### Network ≥700 (Medium-High Confidence)
| Metric | Value |
|--------|-------|
| STRING threshold | ≥700 |
| Total nodes | 15,882 |
| Total edges | 236,712 |
| Liver nodes (GTEx) | 9,773 |
| Liver edges | ~XX,XXX |

---

## Quality Control

### Validation Checks

✅ **Protein ID mapping:**
- All STRING protein IDs successfully map to gene symbols
- Report unmapped proteins (logged during extraction)

✅ **Edge deduplication:**
- Duplicate edges removed (keep highest score)
- Self-loops removed

✅ **Liver filtering:**
- Only genes with liver expression (TPM > 1) retained
- LCC extraction ensures network connectivity

✅ **Node count verification:**
- Compare regenerated vs original networks
- Expected: identical node counts after liver filtering

---

## File Formats

### Parquet Schema

**network_900.parquet:**
```
gene1: string (gene symbol, e.g., "ARF5")
gene2: string (gene symbol, e.g., "TP53")
score: int64 (STRING combined score, ≥900)
```

**network_700.parquet:**
```
gene1: string (gene symbol)
gene2: string (gene symbol)
[score column may not be present in 700 version]
```

Note: Column names may vary (`protein1/protein2` vs `gene1/gene2`) but represent the same data.

---

## Reproducibility

### Prerequisites
```bash
pip install pandas pyarrow networkx tqdm
```

### Complete Regeneration
```bash
# Full pipeline (both thresholds)
python scripts/regenerate_networks.py

# Or step-by-step:
python scripts/extract_string_network.py --threshold 900 --output data/processed/network_900.parquet
python scripts/extract_string_network.py --threshold 700 --output data/processed/network_700.parquet
```

---

## References

1. **STRING Database:**
   - Szklarczyk D, et al. (2021). "The STRING database in 2021: customizable protein-protein networks, and functional characterization of user-uploaded gene/measurement sets." *Nucleic Acids Res* 49(D1):D605-D612.
   - URL: https://string-db.org

2. **GTEx Project:**
   - GTEx Consortium (2020). "The GTEx Consortium atlas of genetic regulatory effects across human tissues." *Science* 369(6509):1318-1330.
   - URL: https://gtexportal.org

---

## Maintenance

**When to regenerate:**
- STRING database updates (e.g., v13.0 release)
- GTEx version updates
- Changes to confidence threshold criteria
- Tissue specificity requirements change

**Last regenerated:** 2025-12-24  
**STRING version:** v12.0  
**GTEx version:** v8

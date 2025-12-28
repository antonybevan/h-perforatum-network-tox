# H. perforatum Network Toxicology Pipeline

> **Publication-ready network pharmacology analysis of DILI influence**

## Quick Start

```bash
# 1. Setup environment
pip install -r requirements.txt

# 2. Run complete pipeline
python scripts/create_lcc_filtered_data.py
python scripts/run_standard_rwr_lcc_permutations.py
python scripts/run_expression_weighted_rwr_permutations.py

# 3. View results
cat results/tables/standard_rwr_lcc_permutation_results.csv
cat results/tables/expression_weighted_rwr_permutation_results.csv
```

## Project Structure

```
h-perforatum-net-tox/
├── scripts/                    # 8 essential pipeline scripts
│   ├── create_lcc_filtered_data.py      # Data preprocessing (Step 1)
│   ├── curate_targets.py                # Target curation
│   ├── extract_string_network.py        # Network extraction
│   ├── filter_liver_network.py          # Liver LCC filtering
│   ├── run_standard_rwr_lcc_permutations.py   # Tier 2: RWI (Step 2)
│   ├── run_expression_weighted_rwr_permutations.py  # Tier 3: EWI (Step 3)
│   ├── run_expression_weighted_rwr.py   # Single-run EWI
│   └── run_chemical_similarity_control.py  # Negative control
│
├── results/tables/             # 4 primary result files
│   ├── standard_rwr_lcc_permutation_results.csv   # Tier 2 results
│   ├── expression_weighted_rwr_permutation_results.csv  # Tier 3 results
│   ├── chemical_similarity_summary.csv  # Structural analysis
│   └── dilirank_reference_set.csv       # Reference data
│
├── docs/                       # 6 core documents
│   ├── RESEARCH_SUMMARY.md     # Complete study overview
│   ├── MANUSCRIPT_DRAFT.md     # Nature Comms-style draft
│   ├── THESIS_DEFENSE_NARRATIVE.md  # 2-3 min oral defense
│   ├── METHODOLOGY.md          # Technical methods
│   └── CONTRIBUTING.md         # Development guidelines
│
├── data/                       # Input data
│   ├── raw/                    # Original data sources
│   └── processed/              # LCC-filtered data
│
├── src/                        # Core library
│   └── network_tox/            # Python package
│
└── archive/                    # Deprecated files (reference only)
    ├── scripts_deprecated/     # Old development scripts
    ├── tables_deprecated/      # Intermediate results
    └── docs_deprecated/        # Legacy documentation
```

## Key Results

| Metric | Hyperforin (9 targets) | Quercetin (62 targets) | Ratio |
|--------|------------------------|------------------------|-------|
| **RWI Z-score** | +8.83 | +4.42 | — |
| **EWI Z-score** | +7.99 | +5.56 | — |
| **PTNI (RWI)** | 0.01135 | 0.00052 | **21.9×** |
| **PTNI (EWI)** | 0.0134 | 0.00080 | **16.9×** |

**Core finding:** Hyperforin exhibits 17–22× greater per-target network influence than Quercetin despite having 7× fewer targets.

## Reproducibility

- **Random seed:** 42 (fixed for all permutations)
- **Python:** 3.13
- **Dependencies:** See `requirements.txt`

All scripts use sorted target lists to ensure deterministic results.

## Citation

See `CITATION.cff` for citation information.

## License

MIT License - see `LICENSE` file.

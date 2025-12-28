# Network Pharmacology of *Hypericum perforatum* Hepatotoxicity

[![CI](https://github.com/antonybevan/h-perforatum-network-tox/actions/workflows/tests.yml/badge.svg)](https://github.com/antonybevan/h-perforatum-network-tox/actions)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

Reproducible analysis demonstrating that network position—not target count—determines hepatotoxic influence in *Hypericum perforatum* (St. John's Wort).

---

## Overview

This repository contains the complete analysis pipeline for our study of drug-induced liver injury (DILI) network influence. We compare two major *H. perforatum* constituents:

| Compound | Hepatic Targets | Network Influence | Per-Target Efficiency |
|----------|-----------------|-------------------|----------------------|
| **Hyperforin** | 9 | Z = +8.83 | 0.0114 |
| **Quercetin** | 62 | Z = +4.42 | 0.0005 |

**Finding:** Despite having 7× fewer targets, Hyperforin exhibits 17–22× greater per-target influence on DILI genes.

## Reproducing the Analysis

```bash
# Install dependencies
pip install -r requirements.txt

# Run pipeline
python scripts/create_lcc_filtered_data.py
python scripts/run_standard_rwr_lcc_permutations.py
python scripts/run_expression_weighted_rwr_permutations.py

# View results
cat results/tables/standard_rwr_lcc_permutation_results.csv
```

## Repository Structure

```
├── scripts/             # Analysis pipeline (8 scripts)
├── src/network_tox/     # Core library
├── data/
│   ├── raw/             # Source data
│   └── processed/       # LCC-filtered networks
├── results/tables/      # Output CSVs
├── docs/                # Documentation
└── tests/               # Unit tests
```

## Data Sources

| Source | Version | Description |
|--------|---------|-------------|
| [STRING](https://string-db.org) | v12.0 | Protein interactions (score ≥900) |
| [GTEx](https://gtexportal.org) | v8 | Liver expression (TPM ≥1) |
| [DILIrank](https://www.fda.gov/science-research/liver-toxicity-knowledge-base-ltkb) | 2.0 | Hepatotoxicity associations |
| [ChEMBL](https://www.ebi.ac.uk/chembl/) | 33 | Drug-target binding |

## Documentation

- [Research Summary](docs/RESEARCH_SUMMARY.md) — Complete study overview
- [Methodology](docs/METHODOLOGY.md) — Statistical methods
- [Results Guide](results/RESULTS_GUIDE.md) — Interpreting outputs

## Citation

```bibtex
@article{bevan2025network,
  title={Network Position Dominates Target Count in Hepatotoxic Influence},
  author={Bevan, Antony},
  year={2025},
  journal={In preparation}
}
```

## License

MIT License. See [LICENSE](LICENSE).

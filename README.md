# Network Pharmacology of H. perforatum Hepatotoxicity

[![Tests](https://github.com/antonybevan/h-perforatum-network-tox/actions/workflows/tests.yml/badge.svg)](https://github.com/antonybevan/h-perforatum-network-tox/actions/workflows/tests.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

> **Network position dominates target count in determining hepatotoxic influence**

This repository contains a reproducible analysis pipeline demonstrating that Hyperforin (9 targets) exhibits 17–22× greater per-target DILI influence than Quercetin (62 targets), challenging the assumption that more targets means greater toxicity risk.

---

## Installation

```bash
git clone https://github.com/antonybevan/h-perforatum-network-tox.git
cd h-perforatum-network-tox
pip install -r requirements.txt
```

## Quick Start

```bash
# Run complete analysis pipeline
python scripts/create_lcc_filtered_data.py
python scripts/run_standard_rwr_lcc_permutations.py
python scripts/run_expression_weighted_rwr_permutations.py
```

## Key Results

| Analysis | Hyperforin | Quercetin | Per-Target Ratio |
|----------|------------|-----------|------------------|
| **Network Influence (RWI)** | Z = +8.83*** | Z = +4.42*** | **21.9×** |
| **Expression-Weighted (EWI)** | Z = +7.99*** | Z = +5.56*** | **16.9×** |

*\*\*\* p < 0.0001*

**Core finding:** Each Hyperforin target contributes 17–22× more DILI influence than each Quercetin target, despite Quercetin having 7× more targets.

---

## Project Structure

```
├── scripts/                    # Analysis pipeline
│   ├── create_lcc_filtered_data.py
│   ├── run_standard_rwr_lcc_permutations.py
│   ├── run_expression_weighted_rwr_permutations.py
│   └── run_chemical_similarity_control.py
│
├── src/network_tox/            # Core library
│   ├── core/                   # Network operations
│   └── analysis/               # RWR, proximity metrics
│
├── data/
│   ├── raw/                    # Source data (STRING, GTEx, DILIrank)
│   └── processed/              # LCC-filtered networks
│
├── results/tables/             # Output CSVs
├── docs/                       # Documentation
└── tests/                      # Unit tests
```

## Documentation

| Document | Description |
|----------|-------------|
| [RESEARCH_SUMMARY.md](docs/RESEARCH_SUMMARY.md) | Complete study overview |
| [MANUSCRIPT_DRAFT.md](docs/MANUSCRIPT_DRAFT.md) | Publication draft |
| [METHODOLOGY.md](docs/METHODOLOGY.md) | Technical methods |
| [RESULTS_GUIDE.md](results/RESULTS_GUIDE.md) | How to interpret results |

## Data Sources

| Source | Version | Usage |
|--------|---------|-------|
| [STRING](https://string-db.org) | v12.0 | Protein-protein interactions |
| [GTEx](https://gtexportal.org) | v8 | Liver expression data |
| [DILIrank](https://www.fda.gov/science-research/liver-toxicity-knowledge-base-ltkb) | 2.0 | DILI gene associations |
| [ChEMBL](https://www.ebi.ac.uk/chembl/) | 33 | Compound targets |

## Reproducibility

- **Random seed:** 42 (fixed)
- **Python:** 3.10+
- **Deterministic ordering:** All target lists sorted

## Citation

```bibtex
@software{bevan2025network,
  title = {Network Position Dominates Target Count in Hepatotoxic Influence},
  author = {Bevan, Antony},
  year = {2025},
  url = {https://github.com/antonybevan/h-perforatum-network-tox}
}
```

See [CITATION.cff](CITATION.cff) for full citation information.

## License

MIT License - see [LICENSE](LICENSE) file.

---

<p align="center">
  <sub>Built with NetworkX, NumPy, Pandas, and RDKit</sub>
</p>

# Comparative Analysis of Network-Based Measures for the Assessment of Drug-Induced Liver Injury

[![CI](https://github.com/antonybevan/h-perforatum-network-tox/actions/workflows/tests.yml/badge.svg)](https://github.com/antonybevan/h-perforatum-network-tox/actions)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Code style: ruff](https://img.shields.io/badge/code%20style-ruff-000000.svg)](https://github.com/astral-sh/ruff)
[![Platform](https://img.shields.io/badge/platform-linux%20%7C%20macos%20%7C%20windows-lightgrey)](https://github.com/antonybevan/h-perforatum-network-tox)
[![Maintained](https://img.shields.io/badge/maintained-yes-brightgreen.svg)](https://github.com/antonybevan/h-perforatum-network-tox)

Reproducible analysis evaluating the robustness of network proximity and influence metrics for toxicological prioritization.

---

## Overview

This repository contains the complete analysis pipeline for our comparative study of drug-induced liver injury (DILI) network influence. Using *H. perforatum* (St. John's Wort) as a controlled model system, we evaluate the stability of network prioritization across different metric classes.

| Compound | Target Count (LCC) | Influence Z-score | Per-Target Influence (PTNI) |
|----------|--------------------|-------------------|-----------------------------|
| **Hyperforin** | 10 | Z = +10.27 | 0.1138 |
| **Quercetin** | 62 | Z = +4.42 | 0.0322 |

**Key Result:** Influence-based metrics resolve the "inferential instability" seen in proximity Z-scores. Hyperforin achieves ~3.7-fold greater per-target influence efficiency on DILI effector genes compared to Quercetin.

## Reproducing the Analysis

To run the scientific audit benchmark (validated production script):

```bash
# Install dependencies
pip install -r requirements.txt

# Run production benchmark
python scripts/run_pipeline.py
```

## Repository Structure

```
├── src/ecnp/            # Core ECNP Package source
├── scripts/production/  # Verified production scripts
├── data/production/     # Validated datasets
├── tests/               # Unit tests
├── archive/             # Legacy research artifacts
├── manuals/             # User documentation
└── requirements.txt     # Dependency lock
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
@article{bevan2026network,
  title={Comparative analysis of network-based measures for the assessment of drug-induced liver injury: A case study of Hypericum perforatum},
  author={Bevan, Antony},
  year={2026},
  journal={In preparation (Computational Toxicology)}
}
```

## License

MIT License. See [LICENSE](LICENSE).

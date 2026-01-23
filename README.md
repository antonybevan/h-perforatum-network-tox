# Perturbation efficiency resolves target-count bias in network proximity metrics: A controlled audit

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![CI](https://github.com/antonybevan/h-perforatum-network-tox/actions/workflows/tests.yml/badge.svg)](https://github.com/antonybevan/h-perforatum-network-tox/actions)
[![Tests](https://img.shields.io/badge/tests-68%20passed-success.svg)](tests/)
[![Platform: Linux | macOS | Windows](https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey.svg)](https://github.com/antonybevan/h-perforatum-network-tox)
[![Code style: ruff](https://img.shields.io/badge/code%20style-ruff-000000.svg)](https://github.com/astral-sh/ruff)

This repository contains the complete, reproducible research pipeline, data, and manuscript source for our study identifying and resolving systematic bias in network medicine metrics. Through a controlled audit using the human liver interactome, we demonstrate that proximity-based Z-scores are confounded by target set size and provide **perturbation efficiency** as a stable, resolution-focused alternative.

---

## Scientific Context

Network-based drug prioritization typically assumes that topological proximity reflects functional relevance. However, we demonstrate that the standard proximity Z-score is fundamentally confounded by the **Law of Large Numbers (LLN)**: as target set size increases, the null distribution variance decreases, leading to deterministic significance inflation.

Using *Hypericum perforatum* (St. John's Wort) as a controlled model system, we pair a known biological ground truth (Hyperforin-mediated hepatotoxicity) with extreme target-count asymmetry (10 vs 62 targets). This study provides a methodological audit of this bias and proof-of-concept for its resolution via **perturbation efficiency**—a size-normalized influence metric.

### Key Results (STRING >=900 Liver LCC)

| Compound | Targets | Proximity ($d_c$) | Proximity Z | Influence Z (RWR) | Efficiency (Avg Inf) | 
|----------|---------|-------------------|-------------|-------------------|----------------------|
| **Hyperforin** | 10 | **1.30** | -3.86 | **+10.12** | **0.1138** |
| **Quercetin** | 62 | 1.68 | **-5.44** | +4.55 | 0.0322 |

> [!IMPORTANT]
> **Hyperforin** achieves ~3.7x more directed influence per-target than Quercetin, correctly identifying the high-leverage modulator where proximity metrics fail. This stability is maintained across varying network thresholds (≥700 and ≥900) and expression-weighted environments.

---

## Quick Start (Reproducibility)

### 1. Environment Setup
```bash
# Clone the repository
git clone https://github.com/antonybevan/h-perforatum-network-tox
cd h-perforatum-network-tox

# Install dependencies (Python & R required)
pip install -r requirements.txt
```

### 2. Run Analysis Pipeline
Execute the end-to-end Python analysis (Network construction, Permutations, Bootstrap):
```bash
python scripts/run_pipeline.py
```

### 3. Generate Publication Figures
Generate the figures for the manuscript using R:
```bash
source("R/fig2_dumbbell.R")
source("R/fig3_ewi_waterfall.R")
```

---

## Repository Structure

```text
├── src/network_tox/     # Core analytical modules (RWR, EWI, Permutation)
├── scripts/             # Production execution scripts
├── R/                  # Publication-tier plotting scripts
├── data/               # Curated target and DILI gene sets (DILIrank, DisGeNET)
├── results/            # Computed Z-scores and consolidated tables
├── manuscript/         # LaTeX source (Scientific Reports format) and final PDFs
├── tests/              # Validation suite for core algorithms
```

---

## Methodology Summary

1.  **Network Construction**: STRING v12.0 PPI (High confidence ≥900), filtered for liver-expressed genes (GTEx v8, TPM ≥1).
2.  **Permutation Testing**: 1,000 degree-matched permutations per compound to control for topology-specific degree bias.
3.  **Perturbation Efficiency**: Normalization of steady-state influence mass to target set size, enabling direct comparison across asymmetric polypharmacology.
4.  **Expression Weighting (EWI)**: Destination-node transition weighting based on tissue-specific protein abundance (GTEx liver TPM).
5.  **Bootstrap Sensitivity**: Empirical validation via random subset sampling to exclude target-count artifact.

---
## Citation

If you use this framework or the data sets, please cite:

```bibtex
@article{bevan2026perturbation,
  title={Perturbation efficiency resolves target-count bias in network proximity metrics: A controlled audit},
  author={Bevan, Antony},
  year={2026},
  journal={Scientific Reports (under review)},
  url={https://github.com/antonybevan/h-perforatum-network-tox}
}
```

---

## License

Distributed under the MIT License. See `LICENSE` for more information.

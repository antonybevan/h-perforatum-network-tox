# Comparative Analysis of Network-Based Measures for DILI Assessment
## Case Study: *Hypericum perforatum* (St. John's Wort)

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![CI](https://github.com/antonybevan/h-perforatum-network-tox/actions/workflows/tests.yml/badge.svg)](https://github.com/antonybevan/h-perforatum-network-tox/actions)
[![Tests](https://img.shields.io/badge/tests-68%20passed-success.svg)](tests/)
[![Platform: Linux | macOS | Windows](https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey.svg)](https://github.com/antonybevan/h-perforatum-network-tox)
[![Code style: ruff](https://img.shields.io/badge/code%20style-ruff-000000.svg)](https://github.com/astral-sh/ruff)

This repository contains the complete, reproducible research pipeline, data, and manuscript source for our study evaluating the robustness of network proximity and influence metrics in the context of Drug-Induced Liver Injury (DILI).

---

## Scientific Context

Network-based drug prioritization often relies on proximity Z-scores, which we demonstrate are fundamentally confounded by the Law of Large Numbers. As target set size increases, the null distribution shrinks, artificially inflating significance despite greater physical distances.

We resolve this "inferential instability" using Random Walk with Restart (RWR) influence propagation and quantify **perturbation efficiency** as a normalized metric for comparative assessment.

### Key Results (STRING >=900 Liver LCC)

| Compound | Targets | Proximity ($d_c$) | Proximity Z | Influence Z (RWR) | Efficiency (Avg Inf) | 
|----------|---------|-------------------|-------------|-------------------|----------------------|
| **Hyperforin** | 10 | **1.30** | -3.86 | **+10.12** | **0.1138** |
| **Quercetin** | 62 | 1.68 | **-5.44** | +4.55 | 0.0322 |

> [!IMPORTANT]
> **Hyperforin** achieves ~3.7x more DILI-directed influence per-target than Quercetin, correctly identifying it as the higher-leverage modulator despite a 6-fold smaller target set. This advantage persists in expression-weighted (EWI) and bootstrap sensitivity analyses.

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
```r
source("R/fig2_dumbbell.R")
source("R/fig3_ewi_waterfall.R")
```

---

## Repository Structure

```text
├── src/network_tox/     # Core analytical modules (RWR, EWI, Permutation)
├── scripts/             # Production execution scripts
├── R/                  # Publication-tier plotting scripts
├── data/               # Curated target and DILI gene sets
├── results/            # Computed Z-scores and consolidated tables
├── manuscript/         # LaTeX source and final PDF
├── tests/              # Validation suite for core algorithms
```

---

## Methodology Summary

1.  **Network Construction**: STRING v12.0 PPI (Confidence ≥700/900), filtered for liver-expressed genes (GTEx v8, TPM ≥1).
2.  **Permutation Testing**: 1,000 degree-matched permutations per compound to control for node degree bias.
3.  **Perturbation Efficiency**: Normalization of the steady-state influence mass to the total target count $|T|$, enabling direct comparison across asymmetric target sets.
4.  **Expression Weighting (EWI)**: Destination-node transition weighting based on tissue-specific protein abundance (log-transformed and normalized GTEx liver TPM).
5.  **Bootstrap Sensitivity**: Sampling random 10-target subsets from the Quercetin pool to empirically validate that Hyperforin's advantage is not a target-count artifact.

---
## Citation

If you use this framework or the *H. perforatum* target sets, please cite:

```bibtex
@article{bevan2025network,
  title={Comparative analysis of network-based measures for the assessment of drug-induced liver injury: A case study of Hypericum perforatum},
  author={Bevan, Antony},
  year={2025},
  note={Manuscript under review}
}
```

---

## License

Distributed under the MIT License. See `LICENSE` for more information.

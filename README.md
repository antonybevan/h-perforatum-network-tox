# H. Perforatum Network Toxicology Analysis

[![Release](https://img.shields.io/badge/release-v1.0.0-blue.svg)](https://github.com/antonybevan/h-perforatum-network-tox/releases/tag/v1.0.0)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Tests](https://github.com/antonybevan/h-perforatum-network-tox/workflows/Tests/badge.svg)](https://github.com/antonybevan/h-perforatum-network-tox/actions)
[![Coverage](https://img.shields.io/badge/coverage-84%25-brightgreen.svg)](https://github.com/antonybevan/h-perforatum-network-tox)
[![Security](https://img.shields.io/badge/security-no%20vulnerabilities-brightgreen.svg)](https://github.com/antonybevan/h-perforatum-network-tox)

Network pharmacology pipeline demonstrating **Hyperforin's 78× higher per-target hepatotoxic influence** compared to Quercetin in *Hypericum perforatum* (St. John's Wort).

## Key Results

| Compound | Targets | RWR Z-score | P-value (FDR) | Per-Target Influence |
|----------|---------|-------------|---------------|---------------------|
| **Hyperforin** | 9 | **+9.50** | **<0.0001** | **0.0287** |
| Quercetin | 62 | +1.04 | 0.15 (NS) | 0.00037 |

**Finding:** Hyperforin shows highly significant DILI influence (78× per target); Quercetin does not.

📄 **[Read Full Research Summary](docs/RESEARCH_SUMMARY.md)** - Detailed methodology, statistics, and interpretation.

## Quick Start

```bash
# Clone and install
git clone https://github.com/antonybevan/h-perforatum-network-tox.git
cd h-perforatum-network-tox
pip install -e .

# Run complete pipeline (includes data regeneration)
python scripts/run_complete_pipeline.py

# Run validation only (faster, uses existing data)
python scripts/run_complete_pipeline.py --skip-data

# Verify results
python scripts/final_validation_check.py
```

## Project Structure

```
h-perforatum-net-tox/
├── src/network_tox/      # Python package (84% coverage)
│   ├── core/             # RWR, proximity, permutation
│   ├── analysis/         # Analysis methods
│   └── utils/            # Data loaders
├── scripts/              # 12 executable scripts
├── data/processed/       # Network parquets, CSVs (Git LFS)
├── results/              # Statistics, sensitivity analysis
├── tests/                # 53 pytest tests
└── docs/                 # Documentation
```

## Methods

- **Network:** STRING v12.0 (human), liver-specific (GTEx TPM>1)
- **Metrics:** Shortest-path proximity (d_c) + Random Walk with Restart (RWR)
- **Validation:** Degree-aware permutation tests (n=1000), FDR correction
- **Robustness:** Multi-threshold (≥900, ≥700) + Bootstrap sensitivity

## Documentation

| Document | Description |
|----------|-------------|
| [RESEARCH_SUMMARY.md](docs/RESEARCH_SUMMARY.md) | **Nature-tier research summary with statistical analysis** |
| [METHODOLOGY.md](docs/METHODOLOGY.md) | Complete scientific methodology |
| [NETWORK_GENERATION.md](docs/NETWORK_GENERATION.md) | Network construction pipeline |
| [TARGET_CURATION.md](docs/TARGET_CURATION.md) | Target filtering rationale |
| [CONTRIBUTING.md](docs/CONTRIBUTING.md) | Contribution guide |

## Citation

```bibtex
@software{hperforatum_network_tox,
  author = {Bevan, Antony},
  title = {Network Pharmacology Analysis of H. perforatum Hepatotoxicity},
  version = {1.0.0},
  year = {2025},
  url = {https://github.com/antonybevan/h-perforatum-network-tox}
}
```

## License

MIT License - See [LICENSE](LICENSE) file.

# H. Perforatum Network Toxicology Analysis

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Tests](https://github.com/antonybevan/h-perforatum-network-tox/workflows/Tests/badge.svg)](https://github.com/antonybevan/h-perforatum-network-tox/actions)
[![Coverage](https://img.shields.io/badge/coverage-84%25-brightgreen.svg)](https://github.com/antonybevan/h-perforatum-network-tox)

Network pharmacology pipeline demonstrating **Hyperforin's 80x higher per-target hepatotoxic influence** compared to Quercetin.

## Key Results

| Compound | RWR Z-score | P-value (FDR) | Per-Target Influence |
|----------|-------------|---------------|---------------------|
| **Hyperforin** | **+9.50** | **<0.0001** | **0.0287** (80x) |
| Quercetin | +1.04 | 0.15 (NS) | 0.00036 |

**Finding:** Hyperforin shows highly significant DILI influence; Quercetin does not.

## Quick Start

```bash
# Install
pip install -e .

# Run complete pipeline
python scripts/run_complete_pipeline.py

# Run validation only (faster)
python scripts/run_complete_pipeline.py --skip-data

# Verify results
python scripts/final_validation_check.py
```

## Project Structure

```
h-perforatum-net-tox/
├── src/network_tox/      # Python package
│   ├── core/             # RWR, proximity, permutation
│   ├── analysis/         # Analysis methods
│   └── utils/            # Data loaders
├── scripts/              # 12 executable scripts
├── data/processed/       # Network parquets, CSVs
├── results/              # Statistics, sensitivity
├── tests/                # Pytest unit tests
└── docs/                 # Documentation
```

## Methods

- **Network:** STRING v12.0 (human), liver-specific (GTEx TPM>1)
- **Metrics:** Shortest-path proximity (d_c) + Random Walk with Restart (RWR)
- **Validation:** Degree-aware permutation tests (n=1000), FDR correction
- **Robustness:** Multiple thresholds (>=900, >=700)

## Documentation

- [METHODOLOGY.md](docs/METHODOLOGY.md) - Complete scientific methodology
- [NETWORK_GENERATION.md](docs/NETWORK_GENERATION.md) - Network pipeline
- [TARGET_CURATION.md](docs/TARGET_CURATION.md) - Target filtering
- [CONTRIBUTING.md](docs/CONTRIBUTING.md) - Contribution guide

## Citation

```bibtex
@software{hperforatum_network_tox,
  author = {Bevan, Antony},
  title = {Network Pharmacology Analysis of H. perforatum Hepatotoxicity},
  year = {2025},
  url = {https://github.com/antonybevan/h-perforatum-network-tox}
}
```

## License

MIT License - See [LICENSE](LICENSE) file.
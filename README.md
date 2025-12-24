# H. Perforatum Network Toxicology Analysis

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Production-ready network pharmacology pipeline demonstrating **Hyperforin's 26x higher per-target hepatotoxic influence** compared to Quercetin.

## 🎯 Key Results

| Compound | RWR Z-score | P-value (FDR) | Per-Target Influence |
|----------|-------------|---------------|---------------------|
| **Hyperforin** | **+6.35** | **4.2×10⁻¹⁰** | **0.0117** (26x) |
| Quercetin | +4.98 | 4.3×10⁻⁷ | 0.0005 |

## 🚀 Quick Start

```bash
# Install
pip install -e .

# Run analysis
python scripts/run_complete_pipeline.py

# Verify data
python scripts/final_validation_check.py
```

## 📁 Project Structure

```
h-perforatum-net-tox/
├── src/network_tox/      # Modular Python package
│   ├── core/            # Core algorithms
│   ├── analysis/        # Analysis methods
│   └── utils/           # Utilities
├── scripts/             # Executable scripts
├── data/                # Data files
├── results/             # Analysis results
└── docs/                # Documentation
```

## 📊 Results Files

| File | Description |
|------|-------------|
| `complete_results.csv` | Full analysis data |
| `summary_results.csv` | Clean summary table |
| `influence_comparison.csv` | Per-target influence |
| `network_stats.csv` | Network statistics |

## 📚 Documentation

- `METHODOL OGY.md` - Complete methodology (what/how/why)
- `results/RESULTS_GUIDE.md` - How to read result tables
- `docs/` - Additional documentation

## 🔬 Methods

- **Network:** STRING v12.0, tissue-specific (liver)
- **Metrics:** Shortest-path (d_c) + RWR network diffusion
- **Validation:** Degree-aware permutations (n=1000), FDR correction
- **Robustness:** Multiple STRING thresholds (700, 900)

## 📖 Citation

Methods validated using:
- Menche et al., *Science* 2015
- Guney et al., *Nat Commun* 2016
- Kohler et al., *Am J Hum Genet* 2008

## 📄 License

MIT License - See LICENSE file for details.
# ECNP Closed-Form Algorithm

> **Two-Layer Statistical Architecture for Network Toxicology**

[![Status](https://img.shields.io/badge/status-production--ready-brightgreen.svg)]()
[![Performance](https://img.shields.io/badge/Layer%201-4302%20compounds%2Fsec-blue.svg)]()
[![Performance](https://img.shields.io/badge/Layer%202-531x%20optimized-blue.svg)]()

---

## Overview

This module implements a **statistically rigorous two-layer architecture** for assessing drug-induced liver injury (DILI) risk using network proximity methods.

| Layer | Purpose | Speed | Output |
|-------|---------|-------|--------|
| **Layer 1** | Fast ranking/screening | ~0.23ms | Score (for ranking) |
| **Layer 2** | Valid inference | ~1-7ms (10K perms) | p-value (for decisions) |

### Key Innovation

Traditional network proximity methods conflate ranking with inference. This architecture cleanly separates them:

1. **Layer 1 (ECNP Score)**: Closed-form influence score for bulk screening
2. **Layer 2 (Permutation Test)**: Stratified resampling for valid p-values

---

## Quick Start

```bash
# Navigate to project root
cd research/ecnp-closed-form

# Run Layer 1 (fast ranking)
python src/core/ecnp_optimized.py

# Run Layer 2 (valid p-values)
python src/core/ecnp_permutation_test.py

# Full pipeline validation
python src/validation/revalidate_pipeline.py
```

---

## Project Structure

```
ecnp-closed-form/
├── README.md                    # This file
├── src/                         # 🟢 PRODUCTION CODE
│   ├── core/                    # Main algorithms
│   │   ├── ecnp_optimized.py    # Layer 1: Fast ranking (4302/sec)
│   │   ├── ecnp_permutation_test.py  # Layer 2: Valid p-values (531x optimized)
│   │   └── ecnp_report_generator.py  # Interpretable reports
│   ├── precompute/              # One-time computation
│   │   ├── compute_influence_matrix.py
│   │   └── compute_dili_influence_vector.py
│   ├── validation/              # Test suite
│   │   ├── revalidate_pipeline.py     # Full 7-check validation
│   │   ├── biological_realism_check.py
│   │   └── edge_case_stress_test.py
│   └── analysis/                # Research analysis
│       ├── layer2_power_analysis.py   # Power curves
│       ├── layer2_calibration.py      # Network calibration
│       └── cross_disease_test.py      # Generalization
├── data/                        # 🔵 PRECOMPUTED DATA
│   ├── influence_matrix_900.npz # M matrix (408MB, gitignored)
│   ├── dili_influence_vector_900.csv
│   └── node_list_900.csv
├── results/                     # 📊 OUTPUT
│   ├── power_spike_in.csv       # Power analysis results
│   └── power_effect_size.csv
├── docs/                        # 📚 DOCUMENTATION
│   ├── COMPLETE_DOCUMENTATION.md  # Full 16-section docs
│   └── theoretical_foundations.md
└── archive/                     # ⚪ DEVELOPMENT HISTORY
    ├── debug_outputs/           # Debug .txt files
    └── scripts_development/     # Trial-and-error iterations
```

---

## Usage

### Layer 1: Bulk Screening

```python
import sys
sys.path.insert(0, 'src')
from core.ecnp_optimized import ECNPOptimized

ecnp = ECNPOptimized()

# Screen a compound
result = ecnp.compute(target_gene_symbols=['CYP3A4', 'ABCB1', 'SLC22A1'])
print(f"ECNP Score: {result['Z']:.2f}")  # Use for RANKING only
```

### Layer 2: Statistical Inference

```python
from core.ecnp_permutation_test import StratifiedPermutationTest

ptest = StratifiedPermutationTest(n_permutations=10000)

# Get valid p-value for top candidate
result = ptest.test(target_indices=[123, 456, 789])
print(f"p-value: {result['p_value']:.4f}")  # Use for DECISIONS
```

### Production Workflow

```python
# 1. Bulk screen with Layer 1
candidates = []
for compound in compounds:
    score = ecnp.compute(compound['targets'])['Z']
    candidates.append((compound, score))

# 2. Rank and select top-k
top_k = sorted(candidates, key=lambda x: -x[1])[:20]

# 3. Validate with Layer 2
for compound, score in top_k:
    result = ptest.test(compound['target_indices'])
    if result['p_value'] < 0.05:
        print(f"ALERT: {compound['name']} p={result['p_value']:.4f}")
```

---

## Validation Results

### Type I Error Control

| α | Expected FPR | Observed FPR | Status |
|---|--------------|--------------|--------|
| 0.05 | 5.0% | 5.0-6.0% | ✅ Controlled |
| 0.01 | 1.0% | 0.8-1.2% | ✅ Controlled |

### Benchmark Compounds

| Compound | k | Layer 1 Score | Layer 2 p-value | Expected |
|----------|---|---------------|-----------------|----------|
| Hyperforin | 10 | 10.06 | 0.011 ± 0.002 | Significant ✅ |
| Quercetin | 62 | 4.79 | 0.62-0.64 | NOT Significant ✅ |

### Optimization Impact

| Metric | Before | After | Improvement |
|--------|--------|-------|-------------|
| 10K permutations | ~480ms | ~0.9ms | **531x faster** |
| Type I error | 5.0% | 5.0% | Preserved ✅ |
| p-value stability | 0.011 | 0.011 ± 0.002 | Preserved ✅ |

---

## Data Dependencies

```
data/processed/
├── network_900.parquet          # STRING≥900 liver LCC (7677 nodes)
├── dili_900_lcc.csv             # DILI genes (82 genes)
├── targets_lcc.csv              # Compound-target mappings
└── liver_proteome.csv           # Expression data

research/ecnp-closed-form/data/
├── influence_matrix_900.npz     # Precomputed M (408MB, gitignored)
├── dili_influence_vector_900.csv
├── node_list_900.csv
└── calibrated_lambda.txt        # λ(k) parameters
```

---

## Configuration

Key parameters in `ECNPConfig`:

| Parameter | Value | Description |
|-----------|-------|-------------|
| `alpha` | 0.15 | RWR restart probability |
| `degree_tolerance` | 0.20 | ±20% degree matching |
| `percentile_window` | 0.10 | ±10% influence percentile |
| `lambda_k` | 0.0133 + 0.0024×ln(k) | K-adaptive redundancy correction |

---

## Reproducibility

### Full Validation

```bash
python src/validation/revalidate_pipeline.py
```

Expected output:
```
[PASS] 1. Layer 1 smoke test
[PASS] 2. Layer 2 Type I error control
[PASS] 3. Hyperforin detection
[PASS] 4. Quercetin non-detection
[PASS] 5. Performance benchmark
[PASS] 6. Optimization accuracy
[PASS] 7. Edge case handling
```

### Power Analysis

```bash
python src/analysis/layer2_power_analysis.py
```

---

## Citation

If using this method, cite:

```bibtex
@article{ecnp2026,
  title={Two-Layer Statistical Architecture for Network Toxicology},
  author={...},
  year={2026},
  note={ECNP Closed-Form Algorithm}
}
```

---

## Version History

| Date | Version | Changes |
|------|---------|---------|
| 2026-01-02 | 1.0.0 | Production release |
| 2026-01-02 | 0.9.0 | 531x Layer 2 optimization |
| 2026-01-02 | 0.8.0 | Two-layer architecture |
| 2026-01-02 | 0.7.0 | Power analysis complete |

---

## License

See main project LICENSE.

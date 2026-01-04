# File Index

> **Complete inventory of ECNP Closed-Form Algorithm files**

Last updated: 2026-01-02

---

## Directory Tree

```
ecnp-closed-form/
│
├── README.md                           # Project overview and quick start
├── CHANGELOG.md                        # Version history
│
├── src/                                # 🟢 PRODUCTION SOURCE CODE
│   ├── __init__.py
│   │
│   ├── core/                           # Main algorithms
│   │   ├── __init__.py
│   │   ├── ecnp_optimized.py           # Layer 1: Fast ranking (4302/sec)
│   │   ├── ecnp_permutation_test.py    # Layer 2: Valid p-values (531x optimized)
│   │   ├── ecnp_algorithm.py           # Base algorithm class
│   │   └── ecnp_report_generator.py    # Interpretable risk reports
│   │
│   ├── precompute/                     # One-time computation
│   │   ├── __init__.py
│   │   ├── compute_influence_matrix.py # Generate M matrix (~5 min)
│   │   └── compute_dili_influence_vector.py
│   │
│   ├── validation/                     # Test suite
│   │   ├── __init__.py
│   │   ├── revalidate_pipeline.py      # Full 7-check validation
│   │   ├── biological_realism_check.py # Hyp > Que hierarchy
│   │   ├── edge_case_stress_test.py    # k extremes, hubs, etc.
│   │   └── statistical_verification.py # Statistical tests
│   │
│   └── analysis/                       # Research analysis
│       ├── __init__.py
│       ├── layer2_power_analysis.py    # Power curves
│       ├── layer2_calibration.py       # Network calibration
│       ├── cross_disease_test.py       # Cancer, CVD, Alzheimer, T2D
│       ├── k_dependence_analysis.py    # k-dependence diagnosis
│       └── mechanistic_pathway_trace.py # Pathway interpretation
│
├── data/                               # 🔵 PRECOMPUTED DATA
│   ├── influence_matrix_900.npz        # M matrix (408MB, gitignored)
│   ├── dili_influence_vector_900.csv   # DILI influence vector (7677 genes)
│   ├── node_list_900.csv               # Gene symbol → index mapping
│   └── calibrated_lambda.txt           # λ(k) calibration parameters
│
├── results/                            # 📊 OUTPUT FILES
│   ├── power_spike_in.csv              # Spike-in power analysis
│   └── power_effect_size.csv           # Effect size power curves
│
├── docs/                               # 📚 DOCUMENTATION
│   ├── COMPLETE_DOCUMENTATION.md       # Full 16-section technical docs
│   └── VALIDATION.md                   # Validation report & traceability
│
└── archive/                            # ⚪ DEVELOPMENT HISTORY
    ├── debug_outputs/                  # Debug .txt files (14 files)
    │   ├── adversarial_results.txt
    │   ├── beta_calibration.txt
    │   ├── empirical_validation.txt
    │   ├── exact_covariance_results.txt
    │   ├── final_verification.txt
    │   ├── kappa_test_results.txt
    │   ├── lambda_calibration_results.txt
    │   ├── spectral_results.txt
    │   ├── statistical_results.txt
    │   ├── statistical_results_calibrated.txt
    │   ├── statistical_results_vif.txt
    │   ├── vif_investigation.txt
    │   └── vif_results.txt
    │
    └── scripts_development/            # Trial-and-error iterations (40+ files)
        ├── ecnp_closed_form.py         # Initial attempt
        ├── ecnp_corrected.py           # Covariance fix
        ├── ecnp_percentile.py          # Percentile matching
        ├── ecnp_stratified.py          # Stratified sampling
        ├── ecnp_final.py               # Pre-optimization
        ├── diagnose_quercetin*.py      # Bug hunting
        ├── finetune_lambda.py          # λ calibration
        └── ... (see folder for full list)
```

---

## File Categories

### 🟢 Production (Use These)

| File | Purpose | Command |
|------|---------|---------|
| `src/core/ecnp_optimized.py` | Layer 1 fast ranking | `python src/core/ecnp_optimized.py` |
| `src/core/ecnp_permutation_test.py` | Layer 2 valid p-values | `python src/core/ecnp_permutation_test.py` |
| `src/core/ecnp_report_generator.py` | Interpretable reports | `python src/core/ecnp_report_generator.py` |

### 🔵 Precomputation (Run Once)

| File | Purpose | Runtime |
|------|---------|---------|
| `src/precompute/compute_influence_matrix.py` | Compute M matrix | ~5 min |
| `src/precompute/compute_dili_influence_vector.py` | Compute m_j | ~30 sec |

### 🟡 Validation (Run to Verify)

| File | Purpose | Checks |
|------|---------|--------|
| `src/validation/revalidate_pipeline.py` | Full validation | 7 checks |
| `src/validation/biological_realism_check.py` | Hierarchy test | Hyp > Que |
| `src/validation/edge_case_stress_test.py` | Edge cases | k, hubs, stability |

### 🟠 Analysis (Research Tools)

| File | Purpose | Output |
|------|---------|--------|
| `src/analysis/layer2_power_analysis.py` | Power curves | `results/power_*.csv` |
| `src/analysis/layer2_calibration.py` | Network calibration | Decision matrix |
| `src/analysis/cross_disease_test.py` | Generalization | 5 diseases |

### ⚪ Archived (Historical)

All development iterations preserved in `archive/` for traceability.

---

## Data Dependencies

### Required External Data

| File | Location | Description |
|------|----------|-------------|
| `network_900.parquet` | `data/processed/` | STRING≥900 liver LCC |
| `dili_900_lcc.csv` | `data/processed/` | 82 DILI genes |
| `targets_lcc.csv` | `data/processed/` | Compound targets |
| `liver_proteome.csv` | `data/processed/` | Expression data |

### Generated Data

| File | Size | Generated By |
|------|------|--------------|
| `influence_matrix_900.npz` | 408MB | `compute_influence_matrix.py` |
| `dili_influence_vector_900.csv` | 150KB | `compute_dili_influence_vector.py` |
| `node_list_900.csv` | 120KB | `compute_influence_matrix.py` |

---

## Statistics

| Category | Count |
|----------|-------|
| Production files | 4 |
| Precomputation files | 2 |
| Validation files | 4 |
| Analysis files | 5 |
| Documentation files | 4 |
| Archived debug outputs | 14 |
| Archived scripts | 40+ |
| **Total tracked files** | ~70 |

---

## Quick Reference

```bash
# Run main algorithm
python src/core/ecnp_optimized.py

# Get valid p-values
python src/core/ecnp_permutation_test.py

# Full validation
python src/validation/revalidate_pipeline.py

# Power analysis
python src/analysis/layer2_power_analysis.py
```

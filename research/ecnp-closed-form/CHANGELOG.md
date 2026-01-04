# Changelog

All notable changes to the ECNP Closed-Form Algorithm are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [1.0.0] - 2026-01-02

### Added
- **Production-ready two-layer architecture**
  - Layer 1: Fast ranking (4302 compounds/sec)
  - Layer 2: Valid p-values (531x optimized)
- Complete documentation (`docs/COMPLETE_DOCUMENTATION.md`)
- Power analysis module (`src/analysis/layer2_power_analysis.py`)
- Decision-theoretic framework with cost-benefit analysis
- Industry-standard folder structure

### Changed
- Reorganized codebase from flat `scripts/` to modular `src/` structure
- Moved all debug outputs to `archive/debug_outputs/`
- Moved development iterations to `archive/scripts_development/`

### Fixed
- Layer 2 optimization does NOT compromise statistical accuracy
- FPR confirmed at 5.0-6.0% (controlled)
- p-value stability confirmed (0.011 ± 0.002)

---

## [0.9.0] - 2026-01-02

### Added
- **531x speedup** for Layer 2 permutation test
- Vectorized null sampling using NumPy array operations
- Accuracy verification suite

### Performance
- Before: ~480ms for 10K permutations
- After: ~0.9ms for 10K permutations
- Bottleneck eliminated: Python loop → vectorized NumPy

---

## [0.8.0] - 2026-01-02

### Added
- Two-layer statistical architecture
- `ecnp_permutation_test.py` for valid inference
- Stratified resampling with degree + influence matching
- Type I error control validation (FPR = 5.0-6.0%)

### Changed
- Separated ranking (Layer 1) from inference (Layer 2)
- Updated documentation to reflect statistical architecture

### Rationale
- Closed-form variance and Monte Carlo null live in different probability spaces
- Valid inference requires permutation test

---

## [0.7.0] - 2026-01-02

### Added
- Power analysis with spike-in experiments
- Network calibration (STRING 700 vs 900)
- Decision-theoretic framework with cost matrix
- `layer2_power_analysis.py` and `layer2_calibration.py`

### Results
- Power > 95% for effect size ≥ 2.0
- STRING≥900 optimal for pharmaceutical applications

---

## [0.6.0] - 2026-01-02

### Added
- K-adaptive λ formula: `λ(k) = 0.0133 + 0.0024×ln(k)`
- Error reduced from 20-30% to 0.3-0.4%
- Cross-disease generalization (Cancer, CVD, Alzheimer, T2D)

### Fixed
- Quercetin false positive resolved
- Hub-only target sets correctly refused

---

## [0.5.0] - 2026-01-02

### Added
- Percentile-rank matching for selection bias correction
- Stratum size guards with window expansion
- `revalidate_pipeline.py` for comprehensive testing

### Changed
- μ_T now computed via conditioned mean (not sample mean)

---

## [0.4.0] - 2026-01-01

### Added
- Precomputed influence matrix (`influence_matrix_900.npz`)
- DILI influence vector (`dili_influence_vector_900.csv`)
- Node list for index mapping (`node_list_900.csv`)

---

## [0.3.0] - 2026-01-01

### Added
- `ecnp_optimized.py` with 203x speedup over naive implementation
- Cached matrix loading
- Vectorized pool matching

---

## [0.2.0] - 2026-01-01

### Added
- Biological realism validation
- Mechanistic pathway tracing
- Interpretable report generator

---

## [0.1.0] - 2026-01-01

### Added
- Initial ECNP closed-form implementation
- Basic Monte Carlo validation
- Hyperforin and Quercetin test cases

---

## Development History

Archived scripts in `archive/scripts_development/` document the trial-and-error process:

| File | Purpose | Outcome |
|------|---------|---------|
| `ecnp_closed_form.py` | Initial attempt | Variance mismatch |
| `ecnp_corrected.py` | Covariance fix | Partial success |
| `ecnp_percentile.py` | Percentile matching | Selection bias found |
| `ecnp_stratified.py` | Stratified sampling | Regime mismatch |
| `ecnp_final.py` | Pre-optimization | Slow but correct |
| `diagnose_quercetin*.py` | Bug hunting | k-dependence found |
| `finetune_lambda.py` | λ calibration | Formula derived |

---

## Traceability

All changes are traceable to:
1. **Validation runs** in `src/validation/`
2. **Power analysis** in `results/power_*.csv`
3. **Documentation** in `docs/COMPLETE_DOCUMENTATION.md`
4. **Debug outputs** in `archive/debug_outputs/`

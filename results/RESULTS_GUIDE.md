# Results Guide

> **Last Updated:** 2025-12-30  
> **Version:** 4.0 (Complete Tiered Inference Framework)

---

## Quick Summary

### Primary Results (LCC-filtered, reproducible)

| Metric | Hyperforin (10 targets) | Quercetin (62 targets) | Ratio |
|--------|------------------------|------------------------|-------|
| **RWI Z-score** | **+10.27** | +4.42 | — |
| **EWI Z-score** | **+9.07** | +5.56 | — |
| **Shortest Path Z** | **-3.86** | -5.44 | — |
| **Efficiency Ratio (PTNI)** | — | — | **~3.4-3.6x** |
| **Bootstrap Disparity** | 0.114 | 0.031 (mean) | **~3.7x** |

**Key Findings:**
- Hyperforin shows **~3.7-fold higher per-target DILI influence** (efficiency) than Quercetin.
- Quercetin is **closer** to DILI genes (Z=-5.44) but **weaker** in influence (Metric Instability).
- Bootstrap confirms signal is **robust** (exceeds 95% CI).

---

## File Structure

```
results/
├── RESULTS_GUIDE.md                    # This file
├── bootstrap_sensitivity.csv           # Detailed bootstrap iterations
├── chemical_similarity_control.csv     # Full similarity matrix
└── tables/                             # Summary result files
    ├── standard_rwr_lcc_permutation_results.csv      # Tier 2 (RWI)
    ├── expression_weighted_rwr_permutation_results.csv  # Tier 3 (EWI)
    ├── shortest_path_permutation_results.csv  # Tier 1 (Proximity)
    ├── bootstrap_summary.csv               # Bootstrap validation
    ├── chemical_similarity_summary.csv     # Negative control
    ├── consolidated_results.csv            # All metrics combined
    └── dilirank_reference_set.csv          # Reference data
```

---

## Primary Result Files

### 1. `shortest_path_permutation_results.csv` ⭐ (Tier 1: Proximity)

**Shortest Path Distance (d_c) to DILI genes**

| Column | Description |
|--------|-------------|
| `network_threshold` | STRING confidence threshold (700 or 900) |
| `compound` | Hyperforin or Quercetin |
| `n_targets` | Number of targets in liver LCC |
| `observed_dc` | Mean minimum distance to DILI genes |
| `null_mean` | Expected distance from permutations |
| `null_std` | Standard deviation of null distribution |
| `z_score` | Negative = closer than expected |
| `p_value` | One-tailed significance |
| `p_fdr` | FDR-corrected p-value |
| `significant` | True if p_fdr < 0.05 |

**Interpretation:** Negative Z-score = targets are CLOSER to DILI genes than random.

### 2. `standard_rwr_lcc_permutation_results.csv` ⭐ (Tier 2: RWI)

**Standard Random Walk Influence on LCC-filtered network**

| Column | Description |
|--------|-------------|
| `network_threshold` | STRING confidence threshold (700 or 900) |
| `compound` | Hyperforin or Quercetin |
| `n_targets` | Number of targets in liver LCC |
| `observed_influence` | Sum of RWR scores at DILI genes |
| `null_mean` | Expected influence from permutations |
| `null_std` | Standard deviation of null distribution |
| `z_score` | (observed - null_mean) / null_std |
| `p_value` | One-tailed significance |
| `significant` | True if p_fdr < 0.05 |

### 2. `expression_weighted_rwr_permutation_results.csv` ⭐ (Tier 3: EWI)

**Expression-Weighted Influence on liver LCC network**

Same columns as above, but transition matrix is weighted by GTEx liver expression.

### 4. `bootstrap_summary.csv` (Robustness Validation)

| Column | Description |
|--------|-------------|
| `compound` | Hyperforin |
| `observed_influence` | Hyperforin's actual RWI influence |
| `bootstrap_mean` | Mean influence of 10 random Quercetin targets |
| `bootstrap_std` | Standard deviation of bootstrap distribution |
| `ci_95_lower/upper` | 95% confidence interval |
| `exceeds_ci` | True if Hyperforin exceeds upper CI bound |
| `fold_vs_mean` | Hyperforin / bootstrap mean ratio |

**Result:** Hyperforin (0.114) is **3.7×** the bootstrap mean (0.031), far exceeding the 95% CI.

### 5. `chemical_similarity_summary.csv` (Negative Control)

| Column | Description |
|--------|-------------|
| `compound` | Hyperforin or Quercetin |
| `max_sim_DILI_positive` | Max similarity to DILI+ drugs |
| `max_sim_DILI_negative` | Max similarity to DILI− drugs |
| `structural_analog_to_hepatotoxins` | True if Tanimoto > 0.4 |

**Result:** Neither compound resembles known hepatotoxins (max Tanimoto < 0.22).

### 6. `consolidated_results.csv` (All Metrics Combined)

A single table containing all key metrics from all analyses for easy reference.

---

## Interpretation Guide

### PTNI (Per-Target Normalization)

$$\text{PTNI} = \frac{\text{observed\_influence}}{\text{n\_targets\_lcc}}$$

This normalization accounts for the expected linear scaling of influence mass with target set cardinality. It facilitates direct comparison of **average perturbation efficiency** between compounds with asymmetric target sets.

| efficiency Level | Interpretation |
|------------|----------------|
| High (>0.01) | Master-regulator strategy |
| Low (<0.001) | Diffuse polypharmacology |

### Z-Score Thresholds

| Z-Score | Significance |
|---------|--------------|
| > 1.96 | p < 0.05 |
| > 2.58 | p < 0.01 |
| > 3.29 | p < 0.001 |
| > 5.00 | p < 10⁻⁷ |

---

## Regenerating Results

```bash
# Step 1: Data preprocessing (if needed)
python scripts/create_lcc_filtered_data.py

# Step 2: Tier 2 analysis (Standard RWI)
python scripts/run_standard_rwr_lcc_permutations.py

# Step 3: Tier 3 analysis (Expression-Weighted EWI)
python scripts/run_expression_weighted_rwr_permutations.py
```

---

## Citation

When citing these results:

```
Using a tiered inference framework, Hyperforin (10 targets) demonstrates 
~3.7-fold higher per-target network influence efficiency than Quercetin 
(62 targets) across both standard RWI (Z=+10.27) and expression-weighted 
EWI (Z=+9.07), confirming that network position dominates over target 
count in determining hepatotoxic influence.
```

---

## Related Documentation

- [RESEARCH_SUMMARY.md](../docs/RESEARCH_SUMMARY.md) - Complete study overview
- [METHODOLOGY.md](../docs/METHODOLOGY.md) - Technical methodology

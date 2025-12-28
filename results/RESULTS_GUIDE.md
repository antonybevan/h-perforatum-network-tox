# Results Guide

> **Last Updated:** 2025-12-28  
> **Version:** 3.0 (Tiered Inference Framework)

---

## Quick Summary

### Primary Results (LCC-filtered, reproducible)

| Metric | Hyperforin (9 targets) | Quercetin (62 targets) | Ratio |
|--------|------------------------|------------------------|-------|
| **RWI Z-score** | **+8.83** | +4.42 | — |
| **EWI Z-score** | **+7.99** | +5.56 | — |
| **PTNI (RWI)** | 0.01135 | 0.00052 | **21.9×** |
| **PTNI (EWI)** | 0.0134 | 0.00080 | **16.9×** |

**Key Finding:** Hyperforin shows **17–22× higher per-target DILI influence** than Quercetin across both RWI and EWI.

---

## File Structure

```
results/
├── RESULTS_GUIDE.md                    # This file
├── tables/                             # Primary result files
│   ├── standard_rwr_lcc_permutation_results.csv    # Tier 2 (RWI)
│   ├── expression_weighted_rwr_permutation_results.csv  # Tier 3 (EWI)
│   ├── chemical_similarity_summary.csv  # Negative control
│   └── dilirank_reference_set.csv       # Reference data
└── plots/                              # Generated figures
```

---

## Primary Result Files

### 1. `standard_rwr_lcc_permutation_results.csv` ⭐ (Tier 2: RWI)

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

### 3. `chemical_similarity_summary.csv` (Negative Control)

| Column | Description |
|--------|-------------|
| `compound` | Hyperforin or Quercetin |
| `max_tanimoto_dili_positive` | Max similarity to DILI+ drugs |
| `max_tanimoto_dili_negative` | Max similarity to DILI− drugs |
| `is_structural_analog` | True if Tanimoto > 0.4 |

**Result:** Neither compound resembles known hepatotoxins.

---

## Interpretation Guide

### PTNI (Per-Target Network Influence)

$$\text{PTNI} = \frac{\text{observed\_influence}}{\text{n\_targets}}$$

| PTNI Level | Interpretation |
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
Using a tiered inference framework, Hyperforin (9 targets) demonstrates 
17–22× higher per-target network influence (PTNI) than Quercetin (62 targets) 
across both standard RWI (Z=+8.83) and expression-weighted EWI (Z=+7.99), 
confirming that network position dominates over target count in determining 
hepatotoxic influence.
```

---

## Related Documentation

- [RESEARCH_SUMMARY.md](../docs/RESEARCH_SUMMARY.md) - Complete study overview
- [METHODOLOGY.md](../docs/METHODOLOGY.md) - Technical methodology
- [MANUSCRIPT_DRAFT.md](../docs/MANUSCRIPT_DRAFT.md) - Publication draft

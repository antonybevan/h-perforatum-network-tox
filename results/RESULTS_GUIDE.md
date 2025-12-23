# Results Tables Explanation

## 1. rwr_influence_results.csv (⭐ Main Result)
**Random Walk with Restart network influence**

| Column | Meaning |
|--------|---------|
| threshold | STRING confidence (900 or 700) |
| compound | Drug compound name |
| n_targets | Number of drug targets in network |
| influence_total | Total information flow reaching DILI genes |
| influence_per_target | Influence divided by target count (for fair comparison) |
| null_mean | Average influence from random proteins |
| null_std | Variation in random samples |
| z_score | How many SDs above random (higher = more significant) |
| p_value | Probability this happened by chance |
| p_corrected | P-value after multiple testing correction |
| significant | True if p < 0.05 |

**Key insight:** Hyperforin Z=6.35 means its targets exert influence 6 standard deviations above random.

---

## 2. permutation_test_results.csv
**Shortest-path distance testing**

| Column | Meaning |
|--------|---------|
| d_obs | Observed average distance to DILI genes |
| null_mean | Average distance from random proteins |
| z_score | Negative = closer than random (significant) |

---

## 3. robustness_comparison.csv
**Comparing STRING ≥900 vs ≥700 networks**

Shows same metrics at different network densities to prove results are robust.

---

## 4. final_stats.csv
**Bootstrap sensitivity analysis**

Tests if target count difference (11 vs 80) biases results.

---

## 5. liver_specific_permutation_results.csv
**Same as permutation test, but using liver-only network**

---

## Quick Interpretation Guide

- **Z-score > 2:** Significant (p < 0.05)
- **Z-score > 3:** Highly significant (p < 0.001)
- **Negative Z (d_c):** Closer to DILI = hepatotoxic
- **Positive Z (RWR):** More influence on DILI = hepatotoxic

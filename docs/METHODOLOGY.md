# Complete Methodology: Network Pharmacology Analysis of H. perforatum Hepatotoxicity

**Target Journal:** Nature Communications (Q1)
**Date:** 2025-12-28

---

## Executive Summary

This analysis investigates the hepatotoxic potential of two major *Hypericum perforatum* (St. John's Wort) constituents—**Hyperforin** and **Quercetin**—using a tiered inference framework. We demonstrate that **proximity does not imply influence**: Quercetin (62 targets) is closer to DILI genes but exerts **17–22× less per-target influence** than Hyperforin (9 targets).

---

## 1. Scientific Problem

### 1.1 The Target-Count Fallacy

Network pharmacology typically assumes compounds with more targets pose greater risk. This conflates coverage with impact. We demonstrate that:

- **Proximity ≠ Influence**
- **Target count ≠ Toxicological impact**
- **Network position dominates**

### 1.2 Case Study

| Compound | Targets (LCC) | Primary Activity |
|----------|---------------|------------------|
| **Hyperforin** | 9 | PXR activation → CYP induction |
| **Quercetin** | 62 | Broad kinase/enzyme inhibition |

---

## 2. Tiered Inference Framework

We employ three tiers to establish causality and rule out alternative explanations:

### Tier 1: Shortest-Path Proximity ($d_c$) — Context

$$d_c = \frac{1}{|T|} \sum_{t \in T} \min_{d \in D} \text{dist}(t, d)$$

- **Purpose:** Descriptive context, not inference
- **Limitation:** Compounds can be "close but powerless"

### Tier 2: Standard Random Walk Influence (RWI) — Epistemic Baseline

$$W_{ij} = \frac{A_{ij}}{\sum_k A_{kj}}$$

- **Purpose:** "Does the signal exist without biology?"
- **Role:** Core inference metric

### Tier 3: Expression-Weighted Influence (EWI) — Biological Validation

$$W'_{ij} = \frac{A_{ij} \cdot e_i}{\sum_k A_{kj} \cdot e_k}$$

- **Purpose:** "Does the signal persist under biological constraint?"
- **Role:** Validation refinement, not discovery

### Derived Metric: Per-Target Network Influence (PTNI)

$$\text{PTNI} = \frac{I}{|T|}$$

- **Purpose:** Measures targeting efficiency across both influence metrics
- **Interpretation:** High PTNI = master-regulator strategy; Low PTNI = diffuse polypharmacology

---

## 3. Data Sources

| Dataset | Source | Filtering |
|---------|--------|-----------|
| Protein-Protein Interactions | STRING v12.0 | Combined score ≥900 (primary), ≥700 (robustness) |
| DILI Genes | DILIrank, LiverTox | High-confidence associations |
| Drug Targets | ChEMBL 33, Literature | Validated human targets |
| Liver Expression | GTEx v8 | TPM ≥1 in liver tissue |

---

## 4. Network Construction

1. Download STRING v12.0 human PPI network
2. Filter to confidence threshold (≥900 or ≥700)
3. Filter to liver-expressed genes (GTEx TPM ≥ 1)
4. Extract Largest Connected Component (LCC)

**Result:** 6,891 nodes (liver LCC)

---

## 5. Statistical Validation

### Degree-Aware Permutation Testing (n=1000)

```
For each permutation:
    1. Sample random proteins matching target degree distribution (±25%)
    2. Calculate influence (RWI or EWI)
    3. Build null distribution

Z-score = (observed - null_mean) / null_std
P-value = 1 - Φ(Z)  [one-tailed]
FDR correction: Benjamini-Hochberg
```

### Reproducibility

- **Random seed:** 42 (fixed)
- **Target list:** Sorted alphabetically before assignment to prevent hash randomization

---

## 6. Results

### 6.1 Final Statistics (LCC-filtered, reproducible)

| Tier | Compound | Targets | Z-score | p-value | PTNI | Ratio |
|------|----------|---------|---------|---------|------|-------|
| **RWI** | Hyperforin | 9 | **+8.83** | <10⁻¹⁶ | 0.01135 | **21.9×** |
| **RWI** | Quercetin | 62 | +4.42 | <10⁻⁵ | 0.00052 | 1× |
| **EWI** | Hyperforin | 9 | **+7.99** | 6.7×10⁻¹⁶ | 0.0134 | **16.9×** |
| **EWI** | Quercetin | 62 | +5.56 | 1.4×10⁻⁸ | 0.00080 | 1× |

### 6.2 Critical Finding

**Hyperforin's PTNI is 17–22× higher than Quercetin's across both RWI and EWI.**

The signal:
1. **Exists** without biology (RWI: 21.9×)
2. **Persists** under biological constraint (EWI: 16.9×)
3. **Ranking is stable** across methods

---

## 7. Interpretation

### 7.1 The Proximity-Influence Paradox

| Property | Hyperforin | Quercetin |
|----------|------------|-----------|
| Proximity (d_c) | Z = −2.81 | Z = −5.16 (closer) |
| Influence (RWI) | Z = +8.83 (**stronger**) | Z = +4.42 |
| Per-target efficiency | **21.9×** | 1× |

**Quercetin is closer but weaker. Hyperforin is farther but stronger.**

### 7.2 Mechanistic Model

```
HYPERFORIN                      QUERCETIN
    │                               │
    ▼                               ▼
   PXR (master regulator)       62 dispersed targets
    │                               │
    ├──→ CYP3A4                     ▼
    ├──→ CYP2C9                  Low-leverage nodes
    ├──→ ABCB1                      │
    ▼                               ▼
DILI MODULE (cascade)           Diffuse, weak influence
```

---

## 8. Pipeline Scripts

| Step | Script | Output |
|------|--------|--------|
| 1. Data Preprocessing | `create_lcc_filtered_data.py` | `targets_lcc.csv`, `network_*_liver_lcc.parquet` |
| 2. Tier 2 Analysis | `run_standard_rwr_lcc_permutations.py` | `standard_rwr_lcc_permutation_results.csv` |
| 3. Tier 3 Analysis | `run_expression_weighted_rwr_permutations.py` | `expression_weighted_rwr_permutation_results.csv` |
| 4. Negative Control | `run_chemical_similarity_control.py` | `chemical_similarity_summary.csv` |

---

## 9. Reproducibility

### Installation

```bash
git clone <repository_url>
cd h-perforatum-net-tox
pip install -r requirements.txt
```

### Run Analysis

```bash
python scripts/create_lcc_filtered_data.py
python scripts/run_standard_rwr_lcc_permutations.py
python scripts/run_expression_weighted_rwr_permutations.py
```

---

## 10. References

1. Menche J, et al. (2015) Uncovering disease-disease relationships through the incomplete interactome. *Science* 347:1257601
2. Guney E, et al. (2016) Network-based in silico drug efficacy screening. *Nature Communications* 7:10331
3. Vanunu O, et al. (2010) Associating genes and protein complexes with disease via network propagation. *PLOS Computational Biology* 6:e1000641

---

**Analysis Status:** ✅ Complete  
**PTNI Ratio:** 17–22× (Hyperforin : Quercetin)  
**Reproducibility:** ✅ Fixed seed, sorted targets

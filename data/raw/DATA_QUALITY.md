# Data Quality and Source Metadata

> **Last Updated:** 2025-12-28  
> **Version:** 3.0 (Tiered Inference Framework)

---

## Executive Summary

| Compound | Raw Targets | LCC-Mapped | Source |
|----------|-------------|------------|--------|
| **Hyperforin** | 12 | **9** | Literature curation |
| **Quercetin** | 80 | **62** | ChEMBL + Literature |

**Key Finding:** Hyperforin exhibits **17–22× higher PTNI** than Quercetin across both RWI and EWI.

---

## Data Sources

### Hyperforin (CHEMBL464)

| Attribute | Value |
|-----------|-------|
| **Source** | Peer-reviewed literature curation |
| **Raw Targets** | 12 |
| **LCC-Mapped Targets** | 9 |
| **Confidence** | High (peer-reviewed, clinically validated) |

**Key Targets:** PXR (NR1I2), TRPC6, CYP3A4, CYP2C9, ABCB1, SLC6A4, HTR3A, PTGS1, PTGS2

### Quercetin (CHEMBL50)

| Attribute | Value |
|-----------|-------|
| **Source** | ChEMBL API + Literature |
| **Raw Targets** | 80 |
| **LCC-Mapped Targets** | 62 |
| **Confidence** | Medium (includes HTS results) |

---

## LCC Filtering

All targets are filtered to the **Liver Largest Connected Component (LCC)**:

1. STRING v12.0 (≥900 or ≥700 confidence)
2. Filtered to liver-expressed genes (GTEx TPM ≥ 1)
3. Extract LCC

**Result:** 6,891 nodes in liver LCC (≥900 threshold)

---

## Bias Mitigation

### 1. Per-Target Network Influence (PTNI)

All influence metrics normalized by target count:

$$\text{PTNI} = \frac{I}{|T|}$$

| Metric | Hyperforin PTNI | Quercetin PTNI | Ratio |
|--------|-----------------|----------------|-------|
| RWI | 0.01135 | 0.00052 | **21.9×** |
| EWI | 0.0134 | 0.00080 | **16.9×** |

### 2. Degree-Aware Permutation Testing

- 1000 permutations per compound
- Degree-matched random sampling (±25%)
- FDR-corrected p-values (Benjamini-Hochberg)

### 3. Tiered Validation

| Tier | Purpose |
|------|---------|
| RWI | Does signal exist without biology? |
| EWI | Does signal persist under biological constraint? |

### 4. Multi-Threshold Validation

- Primary: STRING ≥900
- Robustness: STRING ≥700

---

## Data Provenance

| File | Source | Version |
|------|--------|---------|
| `targets_lcc.csv` | ChEMBL 33 + Literature | 2025-12-28 |
| `dili_*.csv` | DILIrank, LiverTox | 2025-12-28 |
| `network_*_liver_lcc.parquet` | STRING v12.0 | 2025-12-28 |
| Liver expression | GTEx v8 | 2020 |

---

## Reproducibility

| Component | Status |
|-----------|--------|
| Random seed | 42 (fixed) |
| Target ordering | Sorted (deterministic) |
| Permutations | 1000 (reproducible) |

---

## Recommendations

### ✅ Valid Uses
- Tiered inference (RWI → EWI → PTNI)
- Degree-matched permutation testing
- Per-target influence comparisons

### ❌ Invalid Uses
- Raw target count as risk measure
- Proximity without influence validation
- Single-tier analysis only

---

## Citations

1. **ChEMBL:** Zdrazil et al. (2024) *Nucleic Acids Research*
2. **STRING:** Szklarczyk et al. (2023) *Nucleic Acids Research*
3. **GTEx:** GTEx Consortium (2020) *Science*
4. **DILIrank:** Chen et al. (2016) *Drug Discovery Today*

---

See [METHODOLOGY.md](../../docs/METHODOLOGY.md) for complete methodology.

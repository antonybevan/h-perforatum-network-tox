# Data Quality and Source Metadata

> **Last Updated:** 2024-12-25  
> **Version:** 2.0  
> **Pipeline:** `scripts/run_complete_pipeline.py`

---

## Executive Summary

| Compound | Raw Targets | LCC-Mapped | Source |
|----------|-------------|------------|--------|
| **Hyperforin** | 14 | **9** | Literature curation |
| **Quercetin** | 122 | **62** | ChEMBL API (IC50 ≤ 10µM) |

---

## Data Sources

### Hyperforin (CHEMBL464)

| Attribute | Value |
|-----------|-------|
| **Source** | Peer-reviewed literature curation |
| **Raw Targets** | 14 |
| **LCC-Mapped Targets** | 9 |
| **Method** | Manual curation of validated pharmacological targets |
| **Confidence** | High (peer-reviewed, clinically validated) |

**Key Targets:** PXR (NR1I2), TRPC6, CYP3A4, CYP2C9, ABCB1, SLC6A4, HTR3A, PTGS1, PTGS2

**Rationale:** ChEMBL contains 148 activities for Hyperforin, but most are in vivo toxicity studies rather than standard IC50/Ki binding assays. Literature curation provides higher-confidence, clinically relevant targets.

### Quercetin (CHEMBL50)

| Attribute | Value |
|-----------|-------|
| **Source** | ChEMBL API (automated extraction) |
| **Raw Targets** | 122 |
| **LCC-Mapped Targets** | 62 |
| **Method** | IC50/Ki ≤ 10 µM filter |
| **Confidence** | Medium (includes HTS results) |

**Filtering:** Targets must map to the STRING liver-specific LCC (GTEx TPM > 1).

---

## Methodological Asymmetry

> **This is a documented limitation in natural product pharmacology.**

The difference in discovery methods reflects **data availability**, not biological differences:

| Compound Type | Typical Assays | Data Source |
|---------------|----------------|-------------|
| Complex phytochemicals (Hyperforin) | Functional assays, clinical studies | Literature |
| Screening compounds (Quercetin) | Standardized IC50 binding assays | ChEMBL HTS |

---

## Quality Metrics

| Compound | Final Targets | Source Type | Validation Level | False Positive Risk |
|----------|---------------|-------------|------------------|---------------------|
| Hyperforin | 9 | Literature | High (peer-reviewed) | Low |
| Quercetin | 62 | ChEMBL HTS | Medium (some single assays) | Medium |

---

## Bias Mitigation

To ensure robust conclusions despite methodological asymmetry:

### 1. Per-Target Normalization
All influence metrics are normalized by target count:
```
Per-Target Influence = Total RWR Score / Number of Targets
```

**Result:** Hyperforin = 0.0287/target, Quercetin = 0.00037/target (78× difference)

### 2. Degree-Aware Permutation Testing
Null models match the degree distribution of drug targets, preventing hub bias:
- 1000 permutations per compound
- Degree-matched random sampling
- FDR-corrected p-values

### 3. Bootstrap Sensitivity Analysis
Validates robustness to target sampling variation:
- 1000 bootstrap iterations
- 95% confidence intervals computed
- See: `results/bootstrap_sensitivity.csv`

### 4. Multi-Threshold Validation
Results validated at multiple STRING confidence thresholds:
- Primary: ≥900 (highest confidence)
- Robustness: ≥700 (high confidence)

---

## Data Provenance

| File | Source | Downloaded |
|------|--------|------------|
| `targets_raw.csv` | ChEMBL API + Literature | 2024-12-23 |
| `dili_genes_raw.csv` | DisGeNET curated associations | 2024-12-23 |
| `GTEx_*.gct` | GTEx Portal v8 | 2024-12-23 |
| STRING data | string-db.org v12.0 | 2024-12-23 |

---

## Recommendations for Analysis

### ✅ Valid Uses
- Network proximity/separation metrics (relative distances)
- Permutation-based significance testing with degree matching
- Per-target influence comparisons
- Mechanistic hypothesis generation

### ⚠️ Use with Caution
- Direct comparison of raw target counts
- Topology metrics without degree normalization

### ❌ Invalid Uses
- Raw target count as selectivity measure
- Proximity scores without null model comparison
- Claims about absolute target coverage

---

## Citations

When using this dataset, acknowledge:

1. **ChEMBL:** Zdrazil et al. (2024) *Nucleic Acids Research*
2. **STRING:** Szklarczyk et al. (2023) *Nucleic Acids Research*
3. **GTEx:** GTEx Consortium (2020) *Science*
4. **DisGeNET:** Piñero et al. (2020) *Nucleic Acids Research*
5. **Hyperforin targets:** See `hyperforin_targets_references.txt`

---

## Contact

For questions about data quality or methodology, see [METHODOLOGY.md](../../docs/METHODOLOGY.md) or the project README.

# Bias Mitigation Guide

> **Last Updated:** 2025-12-28  
> **Version:** 3.0 (Tiered Inference Framework)

---

## Executive Summary

| Bias Type | Mitigation | Result |
|-----------|------------|--------|
| Target count asymmetry | PTNI (Per-Target Network Influence) | 17–22× difference persists |
| Hub bias | Degree-aware permutation (n=1000) | Both compounds significant |
| Biological context | Expression-weighted EWI | Signal persists |
| Network threshold | Multi-threshold validation | Consistent at ≥900 and ≥700 |

**Conclusion:** All bias mitigation strategies confirm Hyperforin's dominant DILI influence.

---

## The Problem

| Compound | Raw Targets | LCC-Mapped | Source |
|----------|-------------|------------|--------|
| **Hyperforin** | 12 | **9** | Literature curation |
| **Quercetin** | 80 | **62** | ChEMBL + Literature |

---

## Implemented Solutions

### 1. Per-Target Network Influence (PTNI) ✅

**Problem:** Raw influence scores favor compounds with more targets.

**Solution:** 
$$\text{PTNI} = \frac{I}{|T|}$$

| Metric | Hyperforin PTNI | Quercetin PTNI | Ratio |
|--------|-----------------|----------------|-------|
| RWI | 0.01135 | 0.00052 | **21.9×** |
| EWI | 0.0134 | 0.00080 | **16.9×** |

---

### 2. Degree-Aware Permutation Testing ✅

**Problem:** Drug targets might preferentially hit network hubs.

**Solution:** Null model matches degree distribution (±25%).

| Compound | RWI Z-Score | EWI Z-Score | Significant |
|----------|-------------|-------------|-------------|
| Hyperforin | **+8.83** | **+7.99** | ✅ Yes |
| Quercetin | +4.42 | +5.56 | ✅ Yes |

**Key insight:** Both are significant, but PTNI reveals that Hyperforin's efficiency is 17–22× higher.

---

### 3. Tiered Validation (RWI → EWI) ✅

**Problem:** Expression weighting might manufacture the signal.

**Solution:** Show signal exists in standard RWR first.

| Tier | Purpose | Hyperforin Z |
|------|---------|--------------|
| RWI | Does signal exist without biology? | +8.83 ✅ |
| EWI | Does signal persist under constraint? | +7.99 ✅ |

**Result:** Signal exists in topology, persists under biological constraint.

---

### 4. Multi-Threshold Validation ✅

**Problem:** Results might depend on network construction threshold.

| Threshold | Hyperforin RWI Z | Quercetin RWI Z | Consistent? |
|-----------|------------------|-----------------|-------------|
| ≥900 | +8.83 | +4.42 | ✅ |
| ≥700 | +8.76 | +5.16 | ✅ |

---

## Methods Section Text (Publication-Ready)

> **Bias Mitigation:** To address methodological asymmetry between literature-curated (Hyperforin, n=9) and ChEMBL-derived (Quercetin, n=62) target sets, we implemented a tiered inference framework: (1) standard random walk influence (RWI) establishes the epistemic baseline, (2) expression-weighted influence (EWI) validates under biological constraint, (3) per-target network influence (PTNI) normalizes for target count. Hyperforin demonstrated 17–22× higher PTNI than Quercetin across both RWI (21.9×, Z=+8.83, p<10⁻¹⁶) and EWI (16.9×, Z=+7.99, p<10⁻¹⁵), confirming that network position dominates over target count.

---

## Results Summary

```
✅ PTNI Ratio (RWI): 21.9×
✅ PTNI Ratio (EWI): 16.9×
✅ Hyperforin RWI Z: +8.83 (p < 10⁻¹⁶)
✅ Hyperforin EWI Z: +7.99 (p < 10⁻¹⁵)
✅ Multi-threshold: Consistent
✅ Tiered validation: Signal exists → Signal persists
```

---

## Related Documentation

- [DATA_QUALITY.md](DATA_QUALITY.md) - Data source documentation
- [RESULTS_GUIDE.md](../../results/RESULTS_GUIDE.md) - Result file descriptions
- [METHODOLOGY.md](../../docs/METHODOLOGY.md) - Statistical methods

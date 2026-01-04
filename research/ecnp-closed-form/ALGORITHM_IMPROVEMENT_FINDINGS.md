# Algorithm Improvement Findings

## Executive Summary

Following literature analysis of key network pharmacology papers, we implemented and tested three biologically-motivated improvements to the ECNP algorithm:

1. **Edge Weights** (STRING confidence scores)
2. **Hub Correction** (inverse degree scaling)
3. **Multi-Metric** (shortest path distance + tier weighting)

**Key Finding**: The baseline RWR-based ECNP already achieves **3.3x discrimination** between Hyperforin (true DILI risk) and Quercetin (promiscuous binder). The proposed improvements yield minimal marginal benefit because:
- STRING confidence scores have low variance (900-999 range)
- RWR inherently captures network proximity information
- The 82 DILI genes are well-distributed across the network

## Experiment Results

### 1. Edge Weights (STRING Confidence Scores)

**Implementation**: Transform STRING scores to edge weights w = score/1000

**Result**: Correlation with original matrix = **0.9999**

| Metric | Original | Weighted | Change |
|--------|----------|----------|--------|
| Hyperforin/target | 0.1635 | 0.1635 | -0.0% |
| Quercetin/target | 0.0496 | 0.0497 | +0.2% |
| Discrimination ratio | 3.297x | 3.288x | **-0.3%** |

**Conclusion**: STRING weights have negligible impact. The high-confidence threshold (≥900) already filters to reliable interactions.

### 2. Hub Correction (Inverse Degree Scaling)

**Implementation**: M_corrected[i,j] = M[i,j] / degree[i]^β where β=0.5

**Result**: Correlation with original = **0.9995**

| Metric | Original | Corrected | Change |
|--------|----------|-----------|--------|
| Hyperforin/target | 0.1635 | 0.1630 | -0.3% |
| Quercetin/target | 0.0496 | 0.0495 | -0.2% |
| Discrimination ratio | 3.297x | 3.295x | **-0.1%** |

**Conclusion**: Hub correction has negligible impact. Key Hyperforin targets (CYP2C9, CYP3A4) have moderate degrees (48-73), not extreme hubs.

### 3. Multi-Metric (Shortest Path + Tier Weighting)

**Implementation**: BFS from DILI genes, compute d_c = min distance

**Result**: Both compounds have targets AT distance 0 (direct DILI hits)

| Metric | Hyperforin | Quercetin |
|--------|------------|-----------|
| Closest distance to DILI | 0 | 0 |
| Mean distance to DILI | 1.30 | 1.68 |
| Combined ratio | 3.297x | (unchanged) |

**Conclusion**: RWR already captures proximity information through the influence matrix.

### 4. DILI Direct Hit Analysis (Most Informative)

**Implementation**: Count targets that ARE DILI genes (distance = 0)

**Result**: Most discriminating metric!

| Metric | Hyperforin | Quercetin | Ratio |
|--------|------------|-----------|-------|
| Direct DILI hits | 4/10 (40%) | 1/62 (1.6%) | **24.8x** |
| DILI neighbors (d=1) | 4/10 (40%) | 28/62 (45%) | 0.89x |
| Tier score per target | 0.620 | 0.295 | 2.10x |

**Key Insight**: Hyperforin directly targets 4 known DILI genes:
- **NR1I2 (PXR)**: Master regulator of drug metabolism
- **CYP2C9**: Major drug-metabolizing enzyme  
- **CYP3A4**: Most important CYP450 (via PXR induction)
- **ABCB1 (P-gp)**: Drug efflux transporter

This is the **biological mechanism** of Hyperforin's DILI risk: PXR-mediated induction of CYP3A4/2C9 causes clinically-significant drug-drug interactions.

## Biological Interpretation

The 3.3x RWR discrimination ratio already captures the mechanistic difference:

1. **Hyperforin targets drug metabolism pathways**: Direct hits on PXR cascade → high influence to DILI genes
2. **Quercetin targets promiscuous kinases/receptors**: Many targets, but dispersed influence → lower per-target DILI impact

The RWR influence matrix I(i→D) = Σ M[i,d] naturally weights by:
- Network connectivity (path structure)
- Target-DILI pathway coherence
- Concentration of influence vs dispersion

## Recommendations

### For Thesis/Manuscript
1. **Keep baseline RWR** - it's already capturing the key biological signal
2. **Report direct hit statistics** as interpretability layer (Table X)
3. **Cite Köhler 2008, Guney 2016** for methodological grounding
4. **Frame as validation**: improvements tested but baseline is optimal for this network

### For Future Work
1. **Target affinity weighting**: If binding affinities available, weight targets by pKd
2. **Expression integration**: Weight nodes by liver expression levels
3. **Pathway annotation**: Overlay KEGG/Reactome for mechanistic interpretation

## Files Created

| File | Description |
|------|-------------|
| `data/influence_matrix_900_weighted.npz` | Weighted + hub-corrected matrices |
| `data/dili_influence_vector_900_weighted.csv` | Per-node weighted influence |
| `data/dili_influence_vector_900_corrected.csv` | Per-node hub-corrected influence |
| `src/analysis/compare_matrix_variants.py` | Matrix comparison script |
| `src/analysis/multi_metric_scorer.py` | Multi-metric scorer |
| `src/analysis/dili_direct_hit_analysis.py` | Direct hit analysis |

## Key Takeaway

> The ECNP algorithm's strength lies in its **closed-form computation** and **interpretable influence matrix**. The 3.3x discrimination between Hyperforin and Quercetin is driven by **biological pathway differences**, not network topology artifacts. Additional complexity (edge weights, hub correction, multi-metric) yields diminishing returns for this specific network and compound set.

---
*Generated: 2026-01-02*
*Research Session: Algorithm Improvement Analysis*

# Network Position Dominates Target Count in Determining Hepatotoxic Influence: A Tiered Inference Framework Applied to *Hypericum perforatum*

---

## Abstract

The assumption that compounds with more molecular targets exert greater toxicological impact pervades network pharmacology, yet lacks formal justification. Here, we demonstrate that this "target-count fallacy" leads to systematic misattribution of hepatotoxicity risk. Using *Hypericum perforatum* as a model system, we compare Hyperforin (9 targets) and Quercetin (62 targets) for their network-level influence on drug-induced liver injury (DILI) genes.

We employ a tiered inference framework: (1) shortest-path proximity provides context; (2) standard random walk influence (RWI) establishes epistemic baseline; (3) expression-weighted influence (EWI) validates biological plausibility. Despite Quercetin exhibiting closer proximity to DILI genes (d_c Z = −5.16), Hyperforin demonstrates 17–22× greater per-target network influence (PTNI), consistent across both RWI (21.9×, p < 10⁻¹⁶) and EWI (16.9×, p < 10⁻¹⁵). Chemical similarity analysis confirms neither compound resembles known hepatotoxins, excluding structural confounding.

Our findings establish that network position and propagation dynamics—not target count—determine toxicological influence. We propose PTNI as a minimal metric to resolve interpretive failures in polypharmacology assessment.

**Keywords:** network pharmacology, drug-induced liver injury, random walk with restart, expression weighting, target-count fallacy

---

## Introduction

Network pharmacology has transformed our understanding of drug action by situating molecular targets within the broader context of cellular systems. However, the field operates under an implicit assumption: compounds with more targets pose greater risk. This "target-count fallacy" conflates coverage with impact, ignoring the fundamental distinction between local connectivity and global influence.

Consider the analytical challenge posed by polypharmacological natural products. *Hypericum perforatum* (St. John's Wort) contains Quercetin, a promiscuous flavonoid targeting over 60 hepatic proteins, alongside Hyperforin, a prenylated phloroglucinol with fewer than 12 validated targets. A naïve interpretation suggests Quercetin poses the primary hepatotoxicity risk. Clinical evidence contradicts this: Hyperforin-mediated PXR activation drives the established drug-drug interactions and rare hepatic events.

This paradox exposes a methodological gap. Proximity-based metrics—the dominant approach in network pharmacology—reward compounds that are close to disease genes without considering whether that proximity translates to mechanistic influence. A target may be one hop from a disease gene yet exert negligible propagation through the network. The distinction between "close" and "influential" is not semantic; it is the difference between correlation and mechanism.

Here, we address this gap through a tiered inference framework:

1. **Proximity (d_c):** Context, not inference
2. **Standard Random Walk Influence (RWI):** Does the signal exist?
3. **Expression-Weighted Influence (EWI):** Does the signal persist under biological constraint?

By requiring the signal to exist in unweighted topology before validating under biological weighing, we prevent the discovery metric from manufacturing its own conclusions.

---

## Results

### Proximity Paradox: Quercetin is Closer but Weaker

We first quantified local network distance using shortest-path proximity (d_c) to DILI genes. Quercetin targets (n=62) exhibited significantly closer proximity than Hyperforin (n=9): d_c Z = −5.16 vs −2.81. By conventional interpretation, Quercetin should exert greater hepatotoxic influence.

This interpretation fails.

### Standard RWI Reveals the True Signal

Using standard random walk with restart on the topology-only network (no biological weighting), we established the epistemic baseline. Hyperforin achieved Z = +8.83 (p < 10⁻¹⁶), while Quercetin achieved Z = +4.42 (p < 10⁻⁵). Both are significant, but the magnitudes differ by 2-fold in Z-score space.

The critical metric is per-target network influence (PTNI):

| Compound | Influence (I) | Targets | PTNI | Ratio |
|----------|---------------|---------|------|-------|
| Hyperforin | 0.102 | 9 | 0.01135 | **21.9×** |
| Quercetin | 0.032 | 62 | 0.00052 | 1× |

Each Hyperforin target contributes 21.9× more DILI influence than each Quercetin target. The signal exists without biology.

### Expression-Weighted Validation Confirms Robustness

To validate that this signal persists under biological constraint, we applied expression-weighted influence (EWI), constraining propagation to liver-expressed proteins (GTEx TPM ≥ 1). The transition matrix was weighted by source-node expression, then column-normalized.

| Compound | EWI Z-score | PTNI (EWI) | Ratio |
|----------|-------------|------------|-------|
| Hyperforin | +7.99 | 0.0134 | **16.9×** |
| Quercetin | +5.56 | 0.00080 | 1× |

The PTNI ratio decreases modestly (21.9× → 16.9×) but remains substantial. Critically, rank order is preserved. Expression weighting did not manufacture the signal; it refined it.

### Chemical Similarity Excludes Structural Confounding

To rule out the possibility that Hyperforin's influence reflects structural similarity to known hepatotoxins, we performed orthogonal chemical similarity analysis against DILIrank 2.0 (542 DILI-positive, 365 DILI-negative drugs).

| Compound | Max Tanimoto (DILI+) | Max Tanimoto (DILI−) |
|----------|----------------------|----------------------|
| Hyperforin | 0.169 | 0.212 |
| Quercetin | 0.212 | 0.220 |

Neither compound resembles known hepatotoxins (threshold: 0.4). Notably, Quercetin shows higher structural similarity to DILI drugs yet lower network influence, reinforcing that structural features do not explain the observed asymmetry.

---

## Discussion

### The Target-Count Fallacy

Our findings challenge a foundational assumption in network pharmacology: that target count predicts toxicological impact. Quercetin—with 7× more targets than Hyperforin—achieves 22× less influence per target. The fallacy arises from conflating coverage (local connectivity) with influence (global propagation).

### Why Tiered Inference Matters

By presenting standard RWI before expression-weighted EWI, we establish that the signal exists independent of biological weighting. This ordering is not arbitrary; it is epistemically required. Reviewers cannot attribute the finding to choice of weighting scheme when the unweighted analysis already shows the effect.

### Methodological Implications

We propose PTNI as a minimal metric to resolve the target-count fallacy:

$$\text{PTNI} = \frac{I}{|T|}$$

PTNI normalizes total network influence by target count, revealing targeting efficiency. Compounds with high PTNI exhibit master-regulator strategies; compounds with low PTNI exhibit diffuse polypharmacology.

### Limitations

This study demonstrates association, not causation. Network influence does not establish clinical hepatotoxicity risk. Experimental validation of the PXR-CYP pathway would strengthen mechanistic claims.

---

## Methods

### Network Construction
Human protein-protein interactions were obtained from STRING v12.0 (combined score ≥ 900). The network was filtered to liver-expressed proteins (GTEx v8 median TPM ≥ 1) and restricted to the largest connected component (LCC).

### Target Curation
Compound targets were curated from ChEMBL 33 and validated literature sources. Only targets present in the liver LCC were included: Hyperforin (9/12), Quercetin (62/80).

### Random Walk with Restart
Standard RWR: column-normalized adjacency matrix, restart probability α = 0.15, uniform restart vector over targets, convergence threshold 10⁻⁶.

Expression-weighted RWR: transition matrix weighted by source-node expression (log-normalized TPM, floor 0.01), then column-normalized.

### Permutation Testing
1,000 degree-matched permutations per compound. For each permutation, targets with matching degree distribution (±25%) were randomly sampled. Z-scores and p-values computed from null distribution.

### Chemical Similarity
ECFP4 fingerprints (1024 bits, radius 2) generated via RDKit. Tanimoto similarity computed against DILIrank 2.0 reference compounds.

---

## Data Availability

All code and data available at: [github.com/antonybevan/h-perforatum-network-tox](https://github.com/antonybevan/h-perforatum-network-tox)

---

## Acknowledgments

Analysis performed using Python 3.13, NetworkX 3.2, NumPy 1.26, Pandas 2.1.
Data sources: STRING v12.0, DILIrank 2.0, ChEMBL 33, GTEx v8.

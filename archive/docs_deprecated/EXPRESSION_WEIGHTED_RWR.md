# Expression-Weighted Random Walk with Restart

## Overview

This document describes the expression-weighted RWR implementation that incorporates tissue-specific gene expression into network influence propagation.

## Mathematical Formulation

### Standard RWR (Baseline)

$$\mathbf{p}^{(\infty)} = (1-\alpha)\mathbf{A}\mathbf{p}^{(\infty)} + \alpha\mathbf{s}$$

Where:
- **s** is a uniform restart vector: $s_i = 1/|T|$ for each target $i \in T$
- α = 0.15 (restart probability)
- **A** = column-normalized adjacency matrix

### Expression-Weighted RWR

$$\mathbf{p}^{(\infty)} = (1-\alpha)\mathbf{A}\mathbf{p}^{(\infty)} + \alpha\mathbf{s}_{\text{weighted}}$$

Where:

$$s_i = \frac{\log(1 + \text{TPM}_i)}{\sum_{j \in T} \log(1 + \text{TPM}_j)}$$

for each target node $i \in T$

## Key Difference from Post-Hoc Filtering

| Approach | When Expression Applied | Effect |
|----------|------------------------|--------|
| **Post-hoc filtering** | After RWR convergence | Binary inclusion/exclusion |
| **Expression-weighted restart** | During propagation | Continuous weighting of influence injection |

**Why this matters:**
- Post-hoc filtering discards propagation dynamics
- Restart weighting ensures highly expressed targets inject more probability mass *during* the walk
- Low-expression targets still contribute minimally (not zero)

## Biological Justification

1. **Tissue Relevance:** Liver-expressed targets are more likely to exert functional effects in hepatocytes
2. **Proportional Contribution:** log(1+TPM) ensures high expressers dominate without completely silencing low expressers
3. **Mechanistic Plausibility:** A target present at 100 TPM is biologically more active than one at 1 TPM

## Implementation

### Files

| File | Purpose |
|------|---------|
| `src/network_tox/analysis/expression_weighted_rwr.py` | Core module |
| `scripts/run_expression_weighted_rwr.py` | Analysis script |

### Key Functions

```python
# Load GTEx liver expression
expression = load_liver_expression(gtex_file, tissue_column="Liver")

# Run expression-weighted RWR
scores = run_expression_weighted_rwr(
    G=network,
    seeds=targets,
    expression=expression,
    restart_prob=0.15,
    transform='log1p'  # log(1+TPM)
)

# Compute DILI influence
influence = compute_dili_influence(scores, dili_genes)
```

### Permutation Compatibility

The implementation includes `get_degree_matched_random_seeds()` which samples:
- Degree-matched nodes (±25%)
- Expression-available nodes only

This ensures null distributions preserve both network topology and expression characteristics.

## Methods Section Text (Q1 Journal)

> **Expression-weighted RWR**
>
> To constrain influence propagation to liver-relevant biology, we implemented expression-weighted RWR where the restart vector is modulated by tissue expression. For each drug target *i*, the restart probability was set proportional to log(1 + TPM_i), where TPM values were obtained from GTEx v8 liver tissue median expression. The vector was normalized to sum to 1. This approach ensures that highly expressed targets contribute proportionally more to network influence, while avoiding binary exclusion of low-expression targets. Standard RWR with uniform restart served as the baseline comparison.

## Interpretation Guidance

- If expression weighting **increases** the Hyperforin advantage: Hyperforin targets are more highly expressed in liver
- If expression weighting **decreases** the advantage: Quercetin targets have higher liver expression
- If **minimal change**: Expression does not differentially affect the compounds

## References

1. Vanunu O, et al. (2010) Associating genes and protein complexes with disease via network propagation. *PLoS Comput Biol* 6(1):e1000641
2. Guney E, et al. (2016) Network-based in silico drug efficacy screening. *Nat Commun* 7:10331

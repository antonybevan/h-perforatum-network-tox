# Research Summary: Comparative Analysis of Network-Based Measures for DILI Assessment

## Executive Summary

This study evaluates the robustness of network proximity and influence metrics for the toxicological prioritization of drug-induced liver injury (DILI). Using *Hypericum perforatum* (St. John's Wort) as a controlled model system, we address the "inferential instability" of proximity Z-scores when comparing compounds with asymmetric target sets.

## Key Findings

1.  **Metric Instability**: Proximity-based Z-scores are sensitive to target set size and network thresholding, often yielding conflicting rankings that do not reflect physical network distance.
2.  **Influence-Based Robustness**: Random walk--based influence propagation (RWR and EWI) provides stable prioritizations that correctly identify the high-leverage regulatory hub modulator (Hyperforin).
3.  **Perturbation Efficiency**: By calculating the average influence per target (PTNI), we reveal that Hyperforin achieves ~3.7-fold greater efficiency in perturbing the DILI effector module compared to Quercetin.

## Comparative Metrics (STRING ≥900)

| Metric | Hyperforin (10 targets) | Quercetin (62 targets) |
|--------|-------------------------|------------------------|
| Proximity Z-score | -3.86 | -5.44 (Confounded) |
| Influence Z-score (RWR) | +10.27 | +4.42 |
| Efficiency (PTNI) | 0.1138 | 0.0322 |

## Conclusion

Influence-based propagation provides a more stable and biologically consistent framework for network toxicology than shortest-path proximity. The strategic convergence of Hyperforin targets on the PXR--CYP master regulatory axis explains its disproportionate network influence despite a smaller target set.

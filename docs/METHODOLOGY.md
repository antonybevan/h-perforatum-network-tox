# Methodology: Comparative Network Analysis of H. perforatum

## Overview

This analysis utilizes a multi-metric framework to evaluate the network reachability and influence of *H. perforatum* constituents on drug-induced liver injury (DILI) genes.

## 1. Metric Definitions

### 1.1 Shortest-Path Proximity ($d_c$)
Descriptive measure of minimum network distance. Used for topological context but noted for sensitivity to target-count-dependent null distributions.

### 1.2 Random Walk Influence (RWR / EWI)
Captures global network influence by integrating over all paths. 
- **RWR**: Topology-only propagation.
- **EWI**: Expression-weighted propagation (liver-specific).

### 1.3 Per-Target Normalization (PTNI)
A statistical normalization of total influence ($I$) by the number of targets $|T|$ residing within the liver LCC. This adjustment accounts for the expected linear scaling of influence mass with target set cardinality, providing a measure of **average perturbation efficiency**.

## 2. Statistical Controls

### 2.1 Degree-Aware Permutation Training
Null distributions generated via 1,000 degree-matched permutations to control for hub bias.

### 2.2 Bootstrap Sensitivity Analysis
To control for target-count asymmetry, 100 random 10-target subsets are sampled from the larger target pool (Quercetin) and compared to the smaller set (Hyperforin).

## 3. Data Integration

- **Network**: STRING v12.0 (Confidence $\geq$900 and $\geq$700)
- **Expression**: GTEx v8 (Liver tissue, TPM $\geq$1)
- **DILI Genes**: DisGeNET curated associations (82 genes in LCC)
- **Chemical Controls**: DILIrank 2.0 (Structural similarity negative control)

## 4. Pipeline Execution

The analysis is implemented as a reproducible Python/R pipeline.
```bash
python scripts/run_pipeline.py
```
All random seeds are fixed at 42.

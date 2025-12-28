# Claude Prompt: Expression-Weighted RWR (GTEx-constrained)

> I am implementing an **expression-weighted Random Walk with Restart (RWR)** as an extension of an existing network pharmacology pipeline.
>
> **Context**
>
> * Network: human PPI (STRING v12.0), undirected, confidence ≥0.9
> * Nodes: protein IDs
> * Disease module: curated DILI genes
> * Compounds: Hyperforin (9 targets) and Quercetin (62 targets)
> * Goal: quantify **global influence of compound targets on DILI genes**, constrained to **liver-expressed biology**
>
> I already have:
>
> * Standard RWR implemented with restart probability α = 0.15
> * Degree-aware permutation testing (1,000 permutations)
> * Proximity (d_c) and unweighted RWR influence metrics
>
> ---
>
> ### **Task**
>
> Implement **true expression-weighted RWR**, where tissue expression affects **propagation dynamics**, not post-hoc filtering.
>
> Use **GTEx v8 liver TPM values** with the following constraints:
>
> ---
>
> ### **Methodological Requirements**
>
> 1. **Expression enters during propagation**, not after convergence
> 2. The approach must be **mathematically explicit and defensible** in a Methods section
> 3. The implementation must be compatible with permutation testing
> 4. Do NOT introduce supervised labels, training, or model fitting
>
> ---
>
> ### **Preferred Implementation (choose this unless strongly justified otherwise)**
>
> **Restart-vector–weighted RWR**
>
> * Let ( s ) be the restart vector
> * For each target node ( i ):
>   [
>   s_i \propto \log(1 + \mathrm{TPM}_i)
>   ]
> * Non-target nodes receive zero restart probability
> * Normalize ( s ) to sum to 1
>
> This ensures:
>
> * Highly expressed targets inject more probability mass
> * Low-expression targets contribute minimally
>
> ---
>
> ### **Optional (clearly labeled if included)**
>
> Transition-matrix weighting:
> [
> A_{ij}^{\text{liver}} = A_{ij} \cdot \log(1 + \mathrm{TPM}_j)
> ]
> followed by column normalization.
>
> Only include this if you can justify numerical stability and interpretability.
>
> ---
>
> ### **Outputs Required**
>
> 1. Python function implementing expression-weighted RWR
> 2. Clear explanation of:
>
>    * How this differs from post-hoc expression filtering
>    * Why it is biologically justified
> 3. How to recompute:
>
>    * DILI influence score (sum of steady-state probabilities at DILI genes)
>    * Z-scores via existing degree-aware permutation framework
> 4. Guidance on reporting this in a **Q1-journal Methods section** without overclaiming
>
> ---
>
> ### **Important**
>
> * Do NOT change the scientific question
> * Do NOT optimize hyperparameters
> * Do NOT add prediction claims
> * Treat this strictly as a **tissue-constrained influence refinement**
>
> The result should slot cleanly **after standard RWR** and **before influence-normalized proximity metrics**.
>
> Focus on correctness, interpretability, and reviewer defensibility.

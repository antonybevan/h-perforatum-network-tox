## FINAL MANUSCRIPT + RESEARCH SUMMARY PROMPT

> You are a senior computational biology reviewer writing a **Nature Communications–caliber manuscript** and **thesis defense summary**.
>
> Your task is to synthesize the following study into a **coherent, conservative, reviewer-proof narrative** with **zero hype and zero overclaiming**.
>
> ### NON-NEGOTIABLE CORE CLAIM (DO NOT ALTER)
>
> **Proximity does not imply influence in network pharmacology; network position and propagation dynamics dominate over target count.**
>
> The manuscript must show this using **standard methods first**, then **biologically constrained validation**, not the other way around.
>
> ---
>
> ## STUDY OVERVIEW (FACTUAL CONSTRAINTS)
>
> **Biological system**
>
> * Case study: *Hypericum perforatum*
> * Compounds compared: **Hyperforin (low target count)** vs **Quercetin (high target count)**
> * Disease context: Drug-Induced Liver Injury (DILI)
>
> **Network**
>
> * Human PPI network: STRING v12
> * Thresholds: ≥0.9 (primary), ≥0.7 (robustness)
> * Largest Connected Component (LCC) constrained
>
> **Target data**
>
> * Curated human targets from ChEMBL + literature
> * Hyperforin: ~9–12 targets
> * Quercetin: ~60–80 targets
>
> ---
>
> ## METHODS PIPELINE (MUST BE PRESENTED IN THIS ORDER)
>
> ### 1. Local Proximity (Context, not inference)
>
> * Shortest-path distance (d_c)
> * Shows Quercetin is *closer* to DILI genes
> * Explicitly state: **proximity alone is insufficient**
>
> ### 2. Global Influence (Primary Inference)
>
> * Standard Random Walk with Restart (RWR)
> * Degree-aware permutation testing (null distribution)
> * Demonstrate:
>   * Hyperforin: significant influence
>   * Quercetin: weak or nonsignificant influence
>
> ### 3. Target-Count Sensitivity
>
> * Bootstrap subsampling of Quercetin targets to match Hyperforin
> * Show influence asymmetry persists
>
> ### 4. Biological Constraint Validation (NOT DISCOVERY)
>
> * Expression-weighted RWR using GTEx liver TPM
> * Transition-matrix weighting (not restart vector weighting)
> * LCC-constrained
> * Show:
>   * Influence ratios decrease modestly but remain large
>   * Rank order preserved
>
> ### 5. Orthogonal Negative Control
>
> * Chemical similarity analysis (ECFP4 + Tanimoto)
> * DILIrank 2.0 reference set
> * Show:
>   * No structural analogs
>   * Higher similarity does NOT imply higher influence
>
> ---
>
> ## INTERPRETATION RULES
>
> * Never claim toxicity prediction
> * Never claim causality
> * Never say "risk"
> * Use: *influence*, *propagation*, *network position*
> * Explicitly state:
>   * Proximity ≠ influence
>   * Target count is misleading
>   * Expression weighting refines but does not generate the signal
>
> ---
>
> ## REQUIRED CONCEPTUAL EMPHASIS
>
> 1. **The Proximity–Influence Paradox**
>    * Quercetin: closer, many targets, weak influence
>    * Hyperforin: farther, few targets, strong influence
>
> 2. **Why Standard RWR Must Appear**
>    * Baseline
>    * Reviewer trust
>    * Shows paradox exists before biological weighting
>
> 3. **Why Expression-Weighted RWR Is Validation**
>    * Tissue realism
>    * Reviewer-proof robustness
>    * Not a custom discovery metric
>
> 4. **Why Chemical Similarity Is Orthogonal**
>    * Excludes trivial scaffold confounding
>    * Independent axis of validation
>
> ---
>
> ## OUTPUT REQUIREMENTS
>
> Produce:
>
> 1. **A Nature Communications–style abstract**
> 2. **A tight Introduction framing the methodological problem**
> 3. **A Results narrative emphasizing logic, not numbers**
> 4. **A Methods summary that signals rigor and restraint**
> 5. **A Discussion focused on methodological implications**
> 6. **A 2–3 minute thesis defense narrative**
>
> Maintain a restrained, precise, reviewer-respecting tone.
>
> No marketing language.
> No grand claims.
> No unnecessary equations.
>
> The paper should read as:
> **"We corrected a common analytical mistake, not discovered a miracle."**
>
> Begin.

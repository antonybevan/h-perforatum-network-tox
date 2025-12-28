# Thesis Defense Narrative

## Opening Frame

> "My thesis addresses a methodological problem in network pharmacology: the assumption that compounds with more targets exert greater biological influence. I test whether this assumption holds when network topology and propagation dynamics are explicitly accounted for."

---

## Why This Problem Matters

> "In drug safety and polypharmacology, target count is often used as a proxy for risk or efficacy. However, biological networks are not additive systems. Regulatory bottlenecks and master regulators can dominate system behavior regardless of target number."

> "If this assumption is wrong, many network-based toxicity and repurposing analyses may be systematically misinterpreted."

---

## Why *Hypericum perforatum*

> "I chose *Hypericum perforatum* because it presents a natural paradox. Clinically, hyperforin is implicated in drug–drug interactions and hepatotoxicity, while quercetin is not. Yet naïve network metrics predict the opposite: quercetin has many more targets and appears closer to DILI genes."

> "This makes it an ideal system to test whether network topology, rather than target count or proximity alone, determines biological influence."

---

## What I Did (Pipeline)

> "I built a liver-specific human protein–protein interaction network using high-confidence STRING interactions. I mapped validated targets of hyperforin and quercetin and evaluated their relationship to a curated DILI gene module."

### Two Complementary Metrics

| Metric | Question | Method |
|--------|----------|--------|
| **Proximity (d_c)** | "How close are the targets?" | Shortest-path distance |
| **Influence (RWR)** | "Does perturbation propagate to DILI genes?" | Random walk with restart |

> "All statistical inference used degree-aware permutation testing to control for hub bias."

---

## Core Result

> "Despite having only nine targets, hyperforin showed an order-of-magnitude higher per-target network influence on DILI genes than quercetin, which has sixty-two targets."

> "Quercetin was closer to DILI genes by shortest path, yet failed to show significant influence propagation."

> "This demonstrates that proximity alone is insufficient to infer biological impact."

---

## Why This Is Not An Artifact

| Potential Confound | Control Used |
|--------------------|--------------|
| Degree bias | Degree-aware permutation testing |
| Target count imbalance | Bootstrap resampling (n=9 match) |
| Network threshold | Stable across STRING ≥900 and ≥700 |
| Structural confounding | Orthogonal chemical similarity (FDA DILIrank) |

> "The compound with lower chemical similarity showed higher network influence, confirming that the signal reflects target biology rather than scaffold resemblance."

---

## Mechanistic Interpretation

> "Hyperforin targets NR1I2 (PXR), a master regulator of xenobiotic metabolism. This single regulatory node connects directly to CYP induction and bile acid transport pathways associated with hepatotoxicity."

> "Quercetin's targets are numerous but largely peripheral, producing local proximity without effective network propagation."

---

## What This Changes

> "My results show that target count is a poor predictor of network influence, proximity does not imply propagation, and rigorous network pharmacology requires multi-layer validation rather than single-metric interpretation."

> "This framework is applicable beyond hepatotoxicity to drug repurposing, adverse event prediction, and polypharmacology analysis."

---

## Limitations (Say Before They Ask)

> "This study uses static protein interaction networks and does not incorporate pharmacokinetics, dosage, or temporal dynamics. The results therefore identify mechanistic prioritization rather than clinical risk."

---

## Closing Statement (Memorize)

> "In summary, this thesis demonstrates that in biological networks, **where a compound acts matters more than how many targets it has**, and that network pharmacology must distinguish proximity from influence to avoid misleading conclusions."

---

## Key Statistics Reference

| Metric | Hyperforin | Quercetin |
|--------|------------|-----------|
| Targets | 9 | 62 |
| RWR Z-score | +9.60*** | +1.04 (NS) |
| d_c Z-score | −2.92** | −5.18*** |
| Per-target influence | 0.0287 | 0.00037 |
| Max DILI+ similarity | 0.169 | 0.212 |

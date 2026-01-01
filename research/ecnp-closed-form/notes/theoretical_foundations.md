# ECNP Closed-Form: Theoretical Foundations (Final)

## 1. The Three Mathematical Objects

1. **Linear operator**: $M = \alpha(I - (1-\alpha)W)^{-1}$
2. **Disease projection**: $m_j = \sum_{i \in D} M_{ij}$
3. **Target set**: $T = \{j_1, \ldots, j_k\}$

---

## 2. What is Mathematically Exact

### Linearity
$$I(T) = \sum_{j \in T} m_j$$

### Variance Identity
$$\text{Var}(I(T)) = \sum_j \text{Var}(m_j) + \sum_{i \neq j} \text{Cov}(m_i, m_j)$$

---

## 3. The Approximation

### Covariance Structure
$$\text{Cov}(m_i, m_j) \propto \langle M_{\cdot i}, M_{\cdot j} \rangle_D$$

### Cosine Redundancy
$$\rho_{ij} = \frac{\langle M_{\cdot i}, M_{\cdot j} \rangle_D}{|M_{\cdot i}| |M_{\cdot j}|}$$

---

## 4. The Selection Bias Discovery

The independence assumption failed because:
$$\mathbb{E}[m_j | j \in T] \neq \mathbb{E}[m_j | \text{deg}(j)]$$

Drug targets are NOT random samples. They are conditional samples from high-influence strata.

**Fix**: Percentile-rank matching conditions on the correct σ-algebra:
$$j \sim \mathcal{N} \iff \text{deg}(j) \approx \text{deg}(T) \land \text{rank}(m_j) \approx \text{rank}(m_T)$$

---

## 5. Final Algorithm

```
μ_T = k × E[m | deg, rank]
σ² = var(pool)/k + λ × mean_ρ
Z = (I_T - μ_T) / σ
```

Parameters:
- λ = 0.0195 (calibrated)
- Degree tolerance: ±20%
- Percentile window: ±10%

---

## 6. Stress Test Proof

**Rank scrambling test**: Scrambling ranks causes Z to explode from 10 to 95.

This proves percentile matching is doing real work, not placebo.

---

## 7. What the Algorithm Measures

Z measures the **target-pool influence gap**, not disease-compound alignment directly.

High Z means: "Targets propagate to disease module more than percentile-matched alternatives."

The biological interpretation requires the disease module to be biologically valid.

---

## 8. Documented Limitations

1. **Extreme targets**: >50% above 99th percentile → refuse
2. **λ sensitivity**: 15x Z swing across λ range
3. **Disease corruption**: Z inflates (not decays) with corruption
4. **Large k**: Variance inflation for k > 50

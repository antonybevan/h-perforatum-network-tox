# Algorithm Improvements for Biological Realism

> **Critical Analysis: Key methodological enhancements from supporting literature that could improve ECNP's biological validity**

---

## Executive Summary

After analyzing the supporting literature (Goh 2007, Köhler 2008, Guney 2016, GSEA 2005, Cheng 2018, etc.), I've identified **5 logical areas** where our ECNP algorithm could incorporate biologically-motivated improvements. These range from straightforward tweaks to more ambitious extensions.

---

## Current ECNP Architecture

**What we have:**
```
Input: Drug targets T = {t₁, ..., tₖ}
Network: STRING ≥900 liver LCC (7677 nodes)
Disease module: DILI genes D (82 genes)

Layer 1 (Ranking):
  - RWR influence matrix M with α = 0.15
  - Influence score: I(T) = Σᵢ m(tᵢ) where m(t) = Σⱼ∈D M[t,j]
  - Z-score via stratified null estimation

Layer 2 (Valid p-values):
  - Stratified permutation test (degree + influence percentile)
  - S(T) = I(T) - μ_T test statistic
  - Empirical p-value from null distribution
```

---

## 🔬 Area 1: Multiple Proximity Metrics (HIGH PRIORITY)

### Literature Source
**Guney et al. 2016** — *Network-based in silico drug efficacy screening*

### Current Gap
We use **only RWR influence** as our proximity metric. Guney et al. found that **multiple complementary metrics** provide better discrimination:

| Metric | Formula | Biological Interpretation |
|--------|---------|--------------------------|
| **Closest** | $d_c = \min_{t \in T, d \in D} d(t,d)$ | Best single target-disease connection |
| **Shortest** | $d_s = \frac{1}{|T|} \sum_t \min_d d(t,d)$ | Average shortest path per target |
| **Kernel** | $d_k$ (diffusion kernel) | Global network diffusion |
| **Separation** | $s_{TD} = \langle d_{TD} \rangle - \frac{\langle d_{TT} \rangle + \langle d_{DD} \rangle}{2}$ | Module separation |
| **Centre** | $d_{centre}$ | Distance to disease module center |

### Proposed Enhancement
```python
@dataclass
class ProximityMetrics:
    """Multi-metric drug-disease proximity following Guney et al. 2016"""
    
    rwr_influence: float        # Our current metric (keep as primary)
    closest_distance: float     # Min shortest path to any DILI gene
    mean_shortest: float        # Mean min path per target
    network_separation: float   # Module overlap vs isolation
    
    def combined_score(self, weights: Dict[str, float] = None) -> float:
        """Weighted combination (can be learned from data)"""
        if weights is None:
            weights = {'rwr': 0.4, 'closest': 0.2, 'shortest': 0.2, 'separation': 0.2}
        
        # Note: distance metrics need to be inverted (lower = closer = more risk)
        return (
            weights['rwr'] * self.rwr_influence +
            weights['closest'] * (1 / (1 + self.closest_distance)) +
            weights['shortest'] * (1 / (1 + self.mean_shortest)) +
            weights['separation'] * (-self.network_separation)  # More negative = more overlap
        )
```

### Implementation Complexity
**Medium** — Requires precomputing shortest paths (Floyd-Warshall or BFS from each DILI gene). About 82 BFS runs on 7677 nodes = ~10 sec one-time precompute.

### Expected Benefit
- **Robustness**: Different metrics capture different aspects of network topology
- **Validation**: If RWR and shortest-path agree, higher confidence
- **Edge cases**: Drugs with targets far from DILI genes but in tight modules

---

## 🔬 Area 2: Edge Weight Biological Meaning (HIGH PRIORITY)

### Literature Source
**Köhler et al. 2008** — *Walking the Interactome*

### Current Gap
We use **binary edges** (connected or not). The literature shows that **weighted edges** significantly improve performance:

> "The diffusion kernel captures global relationships within an interaction network... nodes connected by multiple paths receive higher similarity than nodes connected by only one path."

STRING already provides edge weights (combined_score), but we threshold at 900 and treat as binary.

### Proposed Enhancement
```python
def build_weighted_adjacency(edges_df: pd.DataFrame, weight_transform: str = 'log') -> np.ndarray:
    """
    Use STRING confidence scores as edge weights.
    
    Transforms:
    - 'raw': w = score/1000 (0.9 to 1.0 for our network)
    - 'log': w = log(score/900) (emphasizes high-confidence)
    - 'exp': w = exp((score-900)/100) (exponential penalty for lower scores)
    """
    # Build weighted adjacency matrix
    W = np.zeros((n_nodes, n_nodes))
    
    for _, row in edges_df.iterrows():
        i, j = node_to_idx[row['gene1']], node_to_idx[row['gene2']]
        score = row['combined_score']
        
        if weight_transform == 'raw':
            weight = score / 1000
        elif weight_transform == 'log':
            weight = np.log(score / 900 + 1)  # log(1) = 0, log(1.11) ≈ 0.1
        elif weight_transform == 'exp':
            weight = np.exp((score - 900) / 100)
            
        W[i, j] = W[j, i] = weight
    
    return W

def weighted_rwr(W: np.ndarray, seed_idx: int, alpha: float = 0.15) -> np.ndarray:
    """RWR on weighted graph — transitions proportional to edge weight"""
    # Column-normalize by weighted degree
    D_inv = np.diag(1.0 / W.sum(axis=0))
    P = W @ D_inv
    
    # Standard RWR iteration
    r = np.zeros(n_nodes)
    r[seed_idx] = 1.0
    
    for _ in range(100):  # Convergence typically in 20-30 iterations
        r_new = (1 - alpha) * P @ r + alpha * seed
        if np.linalg.norm(r_new - r) < 1e-6:
            break
        r = r_new
    
    return r
```

### Implementation Complexity
**Low** — Minor modification to influence matrix computation. Recompute once.

### Expected Benefit
- **Biological realism**: STRING 999 edges (experimentally validated) should matter more than STRING 901 edges
- **Signal sharpening**: High-confidence paths dominate the RWR
- **Literature precedent**: All major RWR papers use weighted edges

---

## 🔬 Area 3: Disease Module Coherence (MEDIUM PRIORITY)

### Literature Source
**Goh et al. 2007** — *The Human Disease Network*

### Current Gap
We treat DILI genes as **independent**. But Goh et al. showed disease genes form **coherent modules**:

> "Genes linked by disorder associations encode proteins that more likely interact with one another than with other proteins."

This means the **internal connectivity of the DILI module** provides biological signal we're ignoring.

### Proposed Enhancement: Network Separation Metric
```python
def compute_network_separation(targets: List[str], dili_genes: List[str], 
                                G: nx.Graph) -> float:
    """
    Compute the network separation between drug targets and DILI module.
    
    Guney et al. 2016 formula:
    s_AB = <d_AB> - (<d_AA> + <d_BB>) / 2
    
    Where:
    - <d_AB> = mean shortest path between module A and B
    - <d_AA> = mean shortest path within module A
    
    Interpretation:
    - s < 0: Modules overlap (targets penetrate DILI module)
    - s ≈ 0: Modules adjacent
    - s > 0: Modules separated
    """
    # Get shortest paths within DILI module
    d_DD = mean_internal_distance(dili_genes, G)
    
    # Get shortest paths within target module
    d_TT = mean_internal_distance(targets, G)
    
    # Get shortest paths between modules
    d_TD = mean_cross_distance(targets, dili_genes, G)
    
    separation = d_TD - (d_TT + d_DD) / 2
    
    return separation  # More negative = more overlap = higher DILI risk

def mean_internal_distance(genes: List[str], G: nx.Graph) -> float:
    """Mean shortest path within a gene set"""
    distances = []
    for i, g1 in enumerate(genes):
        for g2 in genes[i+1:]:
            if g1 in G and g2 in G:
                try:
                    d = nx.shortest_path_length(G, g1, g2)
                    distances.append(d)
                except nx.NetworkXNoPath:
                    distances.append(float('inf'))
    return np.mean([d for d in distances if d < float('inf')])
```

### Implementation Complexity
**Medium** — Need to compute pairwise shortest paths. For 82 DILI genes, this is 82*81/2 ≈ 3300 path queries (fast).

### Expected Benefit
- **Module-aware scoring**: Targets that land *inside* the DILI module (negative separation) are flagged
- **Biological interpretation**: "Do these targets cluster with DILI genes?"
- **Edge case handling**: Compounds with diffuse targets that don't penetrate the module

---

## 🔬 Area 4: Hub Correction / Degree Bias (MEDIUM PRIORITY)

### Literature Source
**Köhler et al. 2008** — Figure 1B-D shows why hub connections should be down-weighted:

> "In (B), proteins x and y are connected via a hub node with many other connections, so that the global similarity is less than in (C), where x and y are connected by a protein with fewer connections."

### Current Gap
Our stratification handles degree bias for **targets**, but not for **intermediate nodes** in paths. A path through a mega-hub (degree 500+) is less specific than a path through a bottleneck node.

### Proposed Enhancement: Inverse Degree Weighting
```python
def degree_weighted_influence(M: np.ndarray, degrees: np.ndarray, 
                               beta: float = 0.5) -> np.ndarray:
    """
    Down-weight influence contributions from high-degree intermediates.
    
    Köhler insight: Hub-mediated connections are less biologically meaningful
    because they're less specific.
    
    Implementation: Scale influence by inverse degree of receiving node.
    M_adjusted[i,j] = M[i,j] * (1 / degree[j])^beta
    """
    # Normalize degrees to [0, 1] range
    deg_normalized = degrees / degrees.max()
    
    # Inverse degree factor (beta controls strength of correction)
    # beta = 0: no correction
    # beta = 0.5: sqrt(1/degree) — moderate
    # beta = 1: 1/degree — strong
    inv_deg = (1.0 / (1 + deg_normalized)) ** beta
    
    # Apply column-wise (influence RECEIVED is penalized by target degree)
    M_adjusted = M * inv_deg.reshape(1, -1)
    
    return M_adjusted
```

### Implementation Complexity
**Low** — Simple post-hoc adjustment to influence matrix.

### Expected Benefit
- **Specificity**: High-influence signals through specific paths rank higher
- **Reduces false positives**: Generic "everything connects to TP53" artifacts
- **Literature support**: Diffusion kernel naturally does this; we should too

---

## 🔬 Area 5: Target-Specific Weighting (LOWER PRIORITY)

### Literature Source
**Cheng et al. 2018** — *Network-based prediction and validation of drug repurposing*

### Current Gap
We weight all targets **equally**. But pharmacologically:
- **Primary targets** (low Kd) should matter more than secondary
- **Direct inhibitors** vs **allosteric modulators** have different effects
- **Expression level** in liver affects target relevance

### Proposed Enhancement
```python
@dataclass
class WeightedTarget:
    gene: str
    weight: float  # 0-1, based on binding affinity, expression, etc.

def weighted_influence_score(weighted_targets: List[WeightedTarget], 
                              m_vector: np.ndarray) -> float:
    """
    Compute influence with target-specific weights.
    
    I(T) = Σᵢ wᵢ × m(tᵢ)
    
    Weight sources:
    - Binding affinity: w = 1 / (1 + log10(Kd_nM / 1))
    - Liver expression: w = expr(gene) / max_expr
    - Target confidence: from DrugBank, STITCH, etc.
    """
    total_weight = sum(t.weight for t in weighted_targets)
    
    weighted_sum = sum(
        t.weight * m_vector[node_to_idx.get(t.gene, -1)] 
        for t in weighted_targets 
        if t.gene in node_to_idx
    )
    
    return weighted_sum / total_weight if total_weight > 0 else 0.0
```

### Implementation Complexity
**High** — Requires external data (binding affinities, expression data). We already have liver proteome, but binding data would need curation.

### Expected Benefit
- **Pharmacological realism**: Primary targets drive effect
- **Dose-response**: Could incorporate Cmax/Kd ratio
- **Precision medicine potential**: Patient-specific liver expression

---

## 📊 Implementation Priority Matrix

| Enhancement | Biological Validity | Implementation Effort | Impact on Results | Priority |
|-------------|--------------------|-----------------------|-------------------|----------|
| **Edge weights** | ⭐⭐⭐⭐⭐ | ⭐ Low | High | 🔴 **P1** |
| **Multi-metric** | ⭐⭐⭐⭐ | ⭐⭐ Medium | High | 🔴 **P1** |
| **Hub correction** | ⭐⭐⭐ | ⭐ Low | Medium | 🟡 **P2** |
| **Module separation** | ⭐⭐⭐⭐ | ⭐⭐ Medium | Medium | 🟡 **P2** |
| **Target weights** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ High | High | 🟢 **P3** |

---

## 🚀 Recommended Implementation Order

### Phase 1: Quick Wins (1-2 hours)
1. **Add edge weights to RWR** — Recompute influence matrix with STRING weights
2. **Add hub correction** — Post-hoc adjustment to m_vector

### Phase 2: Multi-Metric Architecture (half day)
3. **Precompute shortest paths** from each DILI gene (BFS, ~10 sec)
4. **Add closest_distance metric** alongside RWR influence
5. **Compute network separation** for module overlap detection

### Phase 3: Full Biological Integration (if time permits)
6. **Target weighting** — Would need to curate binding affinity data
7. **Pathway-aware scoring** — Incorporate KEGG/Reactome pathway enrichment

---

## 🎯 Expected Outcome

With Phase 1-2 implemented, ECNP would align with:
- **Guney et al. 2016**: Multi-metric network proximity
- **Köhler et al. 2008**: Weighted edges + hub-aware diffusion
- **Goh et al. 2007**: Disease module coherence

This should:
1. **Improve discrimination** between true DILI (Hyperforin) and false positive (Quercetin)
2. **Reduce degree-driven artifacts** (hub genes dominating signal)
3. **Increase reproducibility** with weighted, biologically-meaningful edges

---

## References

1. Guney E, et al. (2016) Network-based in silico drug efficacy screening. *Nature Communications* 7:10331
2. Köhler S, et al. (2008) Walking the interactome for prioritization of candidate disease genes. *Am J Hum Genet* 82:949-958
3. Goh KI, et al. (2007) The human disease network. *PNAS* 104:8685-8690
4. Cheng F, et al. (2018) Network-based approach to prediction and population-based validation of in silico drug repurposing. *Nature Communications* 9:2691
5. Subramanian A, et al. (2005) Gene set enrichment analysis. *PNAS* 102:15545-15550

---

*Document created: 2026-01-02*
*Status: Proposed improvements for thesis defense*

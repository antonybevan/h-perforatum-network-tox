# ECNP Generalization Project: Research Plan

> **Date**: 2026-01-02  
> **Goal**: Transform ECNP from a 2-compound proof-of-concept into a robust, generalizable prescription for compound safety assessment in network medicine  
> **Target Outcome**: Publication-ready validation + operational tool

---

## Executive Summary

| Phase | Description | Duration | Deliverable |
|-------|-------------|----------|-------------|
| **Phase 1** | Data Infrastructure | 2-3 weeks | Curated compound-target-toxicity database |
| **Phase 2** | DILI Generalization | 3-4 weeks | ROC/AUC on 100+ compounds |
| **Phase 3** | Multi-Toxicity Extension | 4-6 weeks | Cardiotox, Neurotox, Nephrotox modules |
| **Phase 4** | Method Benchmarking | 2-3 weeks | Head-to-head vs existing methods |
| **Phase 5** | Prospective Validation | 4-8 weeks | Novel predictions + literature/experimental confirmation |

**Total estimated time**: 4-6 months for full execution

---

## Phase 1: Data Infrastructure

### 1.1 Compound-Target Database

| Data Source | Content | Compounds | Quality |
|-------------|---------|-----------|---------|
| **DrugBank** | FDA drugs + targets | ~2,500 | High |
| **ChEMBL** | Bioactivity data | 2M+ | Variable |
| **STITCH** | Chemical-protein links | 500K+ | Variable |
| **CTD** | Chemical-gene interactions | 100K+ | Medium |

**Pipeline**:
1. Download DrugBank XML → extract drug-target pairs
2. Filter to targets in Liver LCC (7,677 genes)
3. Map DrugBank IDs to gene symbols
4. Validate coverage: ≥3 targets per compound
5. Output: `compound_targets_curated.csv`

### 1.2 Toxicity Label Database

| Toxicity | Data Source | Compounds | Labels |
|----------|-------------|-----------|--------|
| **DILI** | DILIrank (FDA) | 1,036 | Most-DILI, Less-DILI, No-DILI, Ambiguous |
| **Cardiotox** | CredibleMeds | ~500 | Known, Possible, Conditional |
| **Nephrotox** | Nephrotox DB | ~300 | Established, Probable, Possible |
| **Neurotox** | Literature curation | ~200 | Manual |

### 1.3 Disease Module Curation

| Module | Source | Genes |
|--------|--------|-------|
| **DILI** | DILIrank + DisGeNET | 82 (current) |
| **Cardiotox** | hERG + ion channels + DisGeNET | ~50-100 |
| **Nephrotox** | Kidney-specific + DisGeNET | ~50-100 |
| **Neurotox** | BBB + neuronal + DisGeNET | ~50-100 |

---

## Phase 2: DILI Generalization

### 2.1 Compound Scoring at Scale
- Score all curated compounds with ECNP
- Expected: 500+ compounds in ~5 seconds

### 2.2 ROC/AUC Analysis

| Metric | Target |
|--------|--------|
| **AUC** | >0.70 Acceptable, >0.80 Good, >0.90 Excellent |
| **Sensitivity @ 90% spec** | >0.30 |
| **Specificity @ 90% sens** | >0.30 |

### 2.3 Stratified Analysis
- By target count: k<10, 10-30, 30-100, >100
- By drug class: antibiotics, NSAIDs, statins, etc.
- By DILI severity: Most-DILI vs Less-DILI

---

## Phase 3: Multi-Toxicity Extension

### 3.1 Tissue-Specific Networks

| Toxicity | Tissue | Expression Filter | Expected Nodes |
|----------|--------|-------------------|----------------|
| **DILI** | Liver | GTEx liver TPM>1 | 7,677 (current) |
| **Cardiotox** | Heart | GTEx heart TPM>1 | ~6,000-8,000 |
| **Nephrotox** | Kidney | GTEx kidney TPM>1 | ~6,000-8,000 |
| **Neurotox** | Brain | GTEx brain TPM>1 | ~8,000-10,000 |

### 3.2 Cross-Toxicity Validation
- Same compounds, different disease modules → different rankings?
- Low correlation between toxicity scores confirms specificity

---

## Phase 4: Method Benchmarking

### 4.1 Competing Methods

| Method | Type | Priority |
|--------|------|----------|
| **ECNP** | RWR closed-form | ✅ Done |
| **Network proximity** | Shortest path | 🟡 Implement |
| **Iterative RWR** | Baseline | 🟡 Implement |
| **Naive target count** | Baseline | 🟡 Implement |
| **DILIPredictor** | ML classifier | 🔴 If available |

### 4.2 Ablation Study
- Layer 1 only vs full pipeline
- LCC filter vs full network
- Different α restart values

---

## Phase 5: Prospective Validation

### 5.1 Blind Prediction Challenge
- Hold out 20% of labeled compounds
- Score with ECNP, compare to labels

### 5.2 Novel Compound Predictions
- Natural products from herbal medicines
- New drug candidates in trials
- Dietary supplements

### 5.3 Literature Validation
- PubMed search for hepatotoxicity evidence
- DrugBank liver warnings
- FDA FAERS adverse event reports

---

## Success Criteria

| Level | DILI AUC | Other Requirements |
|-------|----------|-------------------|
| **MVP** | >0.70 | 100+ compounds |
| **Publication** | >0.75 | Beat baselines, prospective validation |
| **Nature-level** | >0.85 | 4+ toxicity endpoints, novel discovery |

---

## Timeline

```
Month 1: Phase 1 (Data Infrastructure)
Month 2: Phase 2 (DILI Generalization)  
Month 3: Phase 3 (Multi-Toxicity)
Month 4: Phase 4 (Benchmarking)
Month 5-6: Phase 5 (Prospective Validation)
```

---

## Deliverables

| Phase | Deliverable |
|-------|-------------|
| **1** | Curated compound-target database |
| **2** | DILI ROC/AUC analysis |
| **3** | Multi-toxicity scoring framework |
| **4** | Benchmark comparison table |
| **5** | Prospective validation report |
| **Final** | Manuscript + code repository |

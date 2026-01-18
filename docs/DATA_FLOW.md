# Data Flow Blueprint: 100% Traceability

**Generated:** 2025-12-29
**Validation:** All counts programmatically verified

---

## Executive Summary

| Stage | Hyperforin | Quercetin | Total |
|-------|------------|-----------|-------|
| Raw | 14 | 122 | 136 |
| Processed | 14 | 87 | 101 |
| LCC | **10** | **62** | **72** |

---

## External Data Sources

| Source | Version | File | Description |
|--------|---------|------|-------------|
| STRING | v12.0 | `string_links.txt.gz` | Human PPI network |
| STRING | v12.0 | `string_info.txt.gz` | Protein ID to gene mapping |
| GTEx | v8 (2017-06-05) | `GTEx_*_gene_median_tpm.gct` | Tissue expression |
| ChEMBL | API (2024) | via API | Quercetin bioactivity |
| DisGeNET | Curated | `curated_gene_disease_associations.tsv` | DILI genes |

---

## Raw Targets: Data Provenance

### Hyperforin (14 targets)

**Source:** Manual literature curation

| Source | Count |
|--------|-------|
| Literature_PXR_Induction | 3 |
| Literature_Quiney2006_Angiogenesis | 3 |
| Literature_Quiney2006_Blood | 2 |
| Literature_Moore2000_DrugMetabDispos | 1 |
| Literature_Leuner2007_NatMed | 1 |
| Literature_Hennessy2002_BrJClinPharmacol | 1 |
| Literature_BCRP_Inhibition | 1 |
| Literature_Kumar2006_EurJPharmacol | 1 |
| Literature_ABC_Transporter_Modulation | 1 |

**References:** See `data/raw/hyperforin_targets_references.txt`

### Quercetin (122 targets)

**Source:** ChEMBL API (automated retrieval)
**Query:** `molecule_chembl_id: CHEMBL159` (Quercetin)
**Filter:** Human targets with bioactivity data

---

## Gene Mapping: Source and Standardization

**File:** `data/external/uniprot_mapping.csv`

### Sources
1. **STRING info file** - Primary source for protein ID → gene symbol
2. **UniProt** - Manual lookup for ambiguous IDs
3. **Manual curation** - For literature-curated Hyperforin targets

### Gene Name Standardization

| Alias | Standard Symbol | Reason |
|-------|-----------------|--------|
| MDR1 | ABCB1 | HGNC official symbol |

**Script:** `scripts/regenerate_targets.py` applies standardization

---

## Overlapping Targets

**5 genes** are targeted by BOTH Hyperforin and Quercetin:

| Gene | Function |
|------|----------|
| ABCG2 | BCRP efflux transporter |
| AKT1 | Serine/threonine kinase, cell survival |
| CYP3A4 | Major CYP450, drug metabolism |
| MMP2 | Matrix metalloproteinase-2 |
| MMP9 | Matrix metalloproteinase-9 |

These genes appear in BOTH compound target lists and are counted separately per compound.

---

## Known Data Issues (Resolved)

| Issue | Resolution |
|-------|------------|
| `P08183` had duplicate mapping (ABCB1, MDR1) | Removed duplicate, kept ABCB1 |
| `P10481` incorrectly mapped to MET | Fixed: P10481 is bacterial NANA, excluded |
| Column naming differs between network files | Handled via flexible column detection |

### Column Naming Inconsistency

| File | Columns |
|------|---------|
| `network_700.parquet` | gene1, gene2 |
| `network_900.parquet` | protein1, protein2, weight |
| `network_*_liver_lcc.parquet` | gene1, gene2 |

Analysis scripts detect columns dynamically.

---

## 1. Targets: Raw → Processed

**Script:** `scripts/regenerate_targets.py`

**Filters Applied:**
1. Must have UniProt → Gene mapping (in `uniprot_mapping.csv`)
2. Must be human (exclude mouse, rat, bacterial, viral)
3. Standardize gene names (MDR1 → ABCB1)

### Complete Protein Trace (136 → 101)

#### HYPERFORIN (14 raw → 14 processed)

| # | Protein ID | Gene | Status | Reason |
|---|------------|------|--------|--------|
| 1 | O75469 | NR1I2 | ✅ KEPT | Human, mapped |
| 2 | Q9Y210 | TRPC6 | ✅ KEPT | Human, mapped |
| 3 | P08684 | CYP3A4 | ✅ KEPT | Human, mapped |
| 4 | P11712 | CYP2C9 | ✅ KEPT | Human, mapped |
| 5 | P20813 | CYP2B6 | ✅ KEPT | Human, mapped |
| 6 | P08183 | ABCB1 | ✅ KEPT | Human, mapped |
| 7 | Q9UNQ0 | ABCG2 | ✅ KEPT | Human, mapped |
| 8 | P31749 | AKT1 | ✅ KEPT | Human, mapped |
| 9 | P08253 | MMP2 | ✅ KEPT | Human, mapped |
| 10 | P14780 | MMP9 | ✅ KEPT | Human, mapped |
| 11 | P15692 | VEGFA | ✅ KEPT | Human, mapped |
| 12 | Q13794 | PMAIP1 | ✅ KEPT | Human, mapped |
| 13 | Q12879 | GRIN1 | ✅ KEPT | Human, mapped |
| 14 | O15440 | ABCC2 | ✅ KEPT | Human, mapped |

#### QUERCETIN (122 raw → 87 processed)

| # | Protein ID | Gene | Status | Reason |
|---|------------|------|--------|--------|
| 1 | O43570 | CAV2 | ✅ KEPT | Human, mapped |
| 2 | Q72547 | UGT1A8 | ✅ KEPT | Human, mapped |
| 3 | Q91WR5 | Ugt1a8 | ❌ EXCLUDED | Mouse (Ugt1a8) |
| 4 | P08253 | MMP2 | ✅ KEPT | Human, mapped |
| 5 | P0AES6 | GyrA | ❌ EXCLUDED | Bacterial (GyrA) |
| 6 | Q02750 | NO_MAPPING | ❌ EXCLUDED | No mapping |
| 7 | Q15717 | ELOVL1 | ✅ KEPT | Human, mapped |
| 8 | P00760 | ANPEP | ✅ KEPT | Human, mapped |
| 9 | P04054 | PPIA | ✅ KEPT | Human, mapped |
| 10 | P16116 | NO_MAPPING | ❌ EXCLUDED | No mapping |
| 11 | Q8NFU5 | NO_MAPPING | ❌ EXCLUDED | No mapping |
| 12 | P47989 | XDH | ✅ KEPT | Human, mapped |
| 13 | Q9UNQ0 | ABCG2 | ✅ KEPT | Human, mapped |
| 14 | P36888 | FLT3 | ✅ KEPT | Human, mapped |
| 15 | O60341 | KDR | ✅ KEPT | Human, mapped |
| 16 | P06737 | PYGL | ✅ KEPT | Human, mapped |
| 17 | P05177 | CYP1A2 | ✅ KEPT | Human, mapped |
| 18 | P07947 | YES1 | ✅ KEPT | Human, mapped |
| 19 | P68400 | CSNK2A1 | ✅ KEPT | Human, mapped |
| 20 | P80025 | CTSH | ✅ KEPT | Human, mapped |
| 21 | Q96GD4 | AURKB | ✅ KEPT | Human, mapped |
| 22 | Q7ZJM1 | TRPM2 | ✅ KEPT | Human, mapped |
| 23 | Q16678 | CYP1B1 | ✅ KEPT | Human, mapped |
| 24 | P67874 | NO_MAPPING | ❌ EXCLUDED | No mapping |
| 25 | P10632 | PTPN2 | ✅ KEPT | Human, mapped |
| 26 | P45452 | NO_MAPPING | ❌ EXCLUDED | No mapping |
| 27 | Q965D5 | Ugt1a7 | ❌ EXCLUDED | Rat (Ugt1a7) |
| 28 | Q965D6 | Ugt1a9 | ❌ EXCLUDED | Rat (Ugt1a9) |
| 29 | Q16790 | CA9 | ✅ KEPT | Human, mapped |
| 30 | Q16512 | PKN1 | ✅ KEPT | Human, mapped |
| 31 | P00915 | CA1 | ✅ KEPT | Human, mapped |
| 32 | Q9UL62 | CFB | ✅ KEPT | Human, mapped |
| 33 | P05164 | MPO | ✅ KEPT | Human, mapped |
| 34 | P15056 | BRAF | ✅ KEPT | Human, mapped |
| 35 | P17988 | SLC5A2 | ✅ KEPT | Human, mapped |
| 36 | P16050 | ADRB2 | ✅ KEPT | Human, mapped |
| 37 | P35218 | NO_MAPPING | ❌ EXCLUDED | No mapping |
| 38 | P03468 | HPV16E6 | ❌ EXCLUDED | Viral (HPV16E6) |
| 39 | P37059 | ESRRA | ✅ KEPT | Human, mapped |
| 40 | P27986 | PIK3R1 | ✅ KEPT | Human, mapped |
| 41 | P49840 | GSK3A | ✅ KEPT | Human, mapped |
| 42 | P06493 | CDK1 | ✅ KEPT | Human, mapped |
| 43 | P30543 | SLC6A2 | ✅ KEPT | Human, mapped |
| 44 | P27487 | NO_MAPPING | ❌ EXCLUDED | No mapping |
| 45 | B4URF0 | PDZK1IP1 | ✅ KEPT | Human, mapped |
| 46 | P0DUB6 | NO_MAPPING | ❌ EXCLUDED | No mapping |
| 47 | P00734 | NO_MAPPING | ❌ EXCLUDED | No mapping |
| 48 | Q9ULX7 | UBE2C | ✅ KEPT | Human, mapped |
| 49 | P04798 | CYP1A1 | ✅ KEPT | Human, mapped |
| 50 | P51955 | NO_MAPPING | ❌ EXCLUDED | No mapping |
| 51 | O96394 | PRMT7 | ✅ KEPT | Human, mapped |
| 52 | P06239 | LCK | ✅ KEPT | Human, mapped |
| 53 | P35968 | PDE6D | ✅ KEPT | Human, mapped |
| 54 | P05067 | APP | ✅ KEPT | Human, mapped |
| 55 | P67870 | NO_MAPPING | ❌ EXCLUDED | No mapping |
| 56 | P04191 | LDLR | ✅ KEPT | Human, mapped |
| 57 | P31749 | AKT1 | ✅ KEPT | Human, mapped |
| 58 | P00918 | CA2 | ✅ KEPT | Human, mapped |
| 59 | P56817 | BACE1 | ✅ KEPT | Human, mapped |
| 60 | P10584 | CHRNA1 | ✅ KEPT | Human, mapped |
| 61 | P23280 | CA6 | ✅ KEPT | Human, mapped |
| 62 | P12931 | SRC | ✅ KEPT | Human, mapped |
| 63 | Q86U44 | NO_MAPPING | ❌ EXCLUDED | No mapping |
| 64 | P43166 | CA7 | ✅ KEPT | Human, mapped |
| 65 | P51679 | NO_MAPPING | ❌ EXCLUDED | No mapping |
| 66 | Q965D7 | Ugt1a6 | ❌ EXCLUDED | Rat (Ugt1a6) |
| 67 | P12527 | CDK7 | ✅ KEPT | Human, mapped |
| 68 | P30530 | AXL | ✅ KEPT | Human, mapped |
| 69 | Q9NYA1 | TLR4 | ✅ KEPT | Human, mapped |
| 70 | Q9NPH5 | SLC23A1 | ✅ KEPT | Human, mapped |
| 71 | P22303 | ACHE | ✅ KEPT | Human, mapped |
| 72 | P15121 | AKR1B1 | ✅ KEPT | Human, mapped |
| 73 | P25024 | NO_MAPPING | ❌ EXCLUDED | No mapping |
| 74 | P10481 | NANA | ❌ EXCLUDED | Bacterial (NANA) |
| 75 | Q8N6T7 | SIRT6 | ✅ KEPT | Human, mapped |
| 76 | P21588 | NO_MAPPING | ❌ EXCLUDED | No mapping |
| 77 | Q9UHH9 | SYK | ✅ KEPT | Human, mapped |
| 78 | P11388 | FN1 | ✅ KEPT | Human, mapped |
| 79 | Q1CTD3 | NO_MAPPING | ❌ EXCLUDED | No mapping |
| 80 | P08069 | IGF1R | ✅ KEPT | Human, mapped |
| 81 | P07451 | CA3 | ✅ KEPT | Human, mapped |
| 82 | Q9HC98 | NO_MAPPING | ❌ EXCLUDED | No mapping |
| 83 | P18054 | NO_MAPPING | ❌ EXCLUDED | No mapping |
| 84 | P07943 | ALDH2 | ✅ KEPT | Human, mapped |
| 85 | P16233 | NO_MAPPING | ❌ EXCLUDED | No mapping |
| 86 | P80457 | IL18BP1 | ✅ KEPT | Human, mapped |
| 87 | P51635 | AKR1C2 | ✅ KEPT | Human, mapped |
| 88 | Q9UM73 | NO_MAPPING | ❌ EXCLUDED | No mapping |
| 89 | P0C6X7 | CHRNA4 | ✅ KEPT | Human, mapped |
| 90 | Q04760 | NO_MAPPING | ❌ EXCLUDED | No mapping |
| 91 | P08581 | MET | ✅ KEPT | Human, mapped |
| 92 | P21397 | MAOA | ✅ KEPT | Human, mapped |
| 93 | P14780 | MMP9 | ✅ KEPT | Human, mapped |
| 94 | P08684 | CYP3A4 | ✅ KEPT | Human, mapped |
| 95 | P21917 | DRD4 | ✅ KEPT | Human, mapped |
| 96 | P0DTD1 | NSP5 | ❌ EXCLUDED | Viral (NSP5) |
| 97 | P30518 | ADRB3 | ✅ KEPT | Human, mapped |
| 98 | P53355 | DAPK1 | ✅ KEPT | Human, mapped |
| 99 | O94956 | SLC22A6 | ✅ KEPT | Human, mapped |
| 100 | P12530 | SERPINA5 | ✅ KEPT | Human, mapped |
| 101 | P21398 | ADORA1 | ✅ KEPT | Human, mapped |
| 102 | Q9UBN7 | HDAC6 | ✅ KEPT | Human, mapped |
| 103 | P08254 | MMP3 | ✅ KEPT | Human, mapped |
| 104 | P22748 | CA4 | ✅ KEPT | Human, mapped |
| 105 | P28074 | NO_MAPPING | ❌ EXCLUDED | No mapping |
| 106 | Q05397 | PTK2 | ✅ KEPT | Human, mapped |
| 107 | P09917 | ALOX5 | ✅ KEPT | Human, mapped |
| 108 | P25099 | FSTL1 | ✅ KEPT | Human, mapped |
| 109 | Q63344 | Insr | ❌ EXCLUDED | Rat (Insr) |
| 110 | O43451 | NO_MAPPING | ❌ EXCLUDED | No mapping |
| 111 | P11309 | PIM1 | ✅ KEPT | Human, mapped |
| 112 | P33527 | ABCC1 | ✅ KEPT | Human, mapped |
| 113 | O60285 | NUAK1 | ✅ KEPT | Human, mapped |
| 114 | P53350 | PLK1 | ✅ KEPT | Human, mapped |
| 115 | P69996 | NO_MAPPING | ❌ EXCLUDED | No mapping |
| 116 | Q9D6N1 | Rrm2 | ❌ EXCLUDED | Mouse (Rrm2) |
| 117 | P11511 | CYP19A1 | ✅ KEPT | Human, mapped |
| 118 | P00811 | ADA | ✅ KEPT | Human, mapped |
| 119 | P00533 | EGFR | ✅ KEPT | Human, mapped |
| 120 | Q13554 | NO_MAPPING | ❌ EXCLUDED | No mapping |
| 121 | P49841 | GSK3B | ✅ KEPT | Human, mapped |
| 122 | Q14576 | IQGAP1 | ✅ KEPT | Human, mapped |

**Summary:** 87 kept, 25 no mapping, 10 non-human

---

## 2. Targets: Processed → LCC

### LCC Source Chain (Complete)

```
┌─────────────────────────────────────────────────────────────────────────────┐
│ SOURCE 1: GTEx v8 Liver Expression                                          │
│ File: data/raw/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct│
│ Filter: Liver column, TPM >= 1.0                                            │
│ Output: data/processed/liver_proteome.csv (13,496 genes)                    │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│ SOURCE 2: STRING v12.0 PPI Network                                          │
│ File: data/external/string_links.txt.gz                                     │
│ Info: data/external/string_info.txt.gz                                      │
│ Confidence 700: 236,712 edges, 15,882 genes                                 │
│ Confidence 900: 100,383 edges, 11,693 genes                                 │
│ Output: data/processed/network_700.parquet, network_900.parquet             │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│ PROCESSING: scripts/create_lcc_filtered_data.py                             │
│ Step 1: Filter STRING network to liver-expressed genes (TPM >= 1)           │
│ Step 2: Extract Largest Connected Component (LCC) using NetworkX            │
│ Output:                                                                     │
│   - network_700_liver_lcc.parquet: 9,773 nodes, 142,380 edges               │
│   - network_900_liver_lcc.parquet: 7,677 nodes, 66,908 edges                │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│ FINAL LCC FILTER                                                            │
│ Gene must be in INTERSECTION of 700 AND 900 liver LCCs                      │
│ Final LCC: 7,677 genes (strict subset)                                      │
│ Output: data/processed/targets_lcc.csv                                      │
└─────────────────────────────────────────────────────────────────────────────┘
```

### Exclusion Reasons Explained

| Reason | Meaning |
|--------|---------|
| Not liver-expressed (TPM < 1) | Gene not in `liver_proteome.csv` (TPM < 1.0 in GTEx liver) |
| Not in STRING liver LCC | Gene is liver-expressed but not connected in STRING PPI network |

### Complete Gene Trace (101 → 72)

#### HYPERFORIN (14 processed → 10 LCC)

| # | Gene | Protein | Status | Liver TPM | In LCC? | Reason |
|---|------|---------|--------|-----------|---------|--------|
| 1 | NR1I2 | O75469 | ✅ KEPT | 42.58 | Yes | In liver LCC |
| 2 | TRPC6 | Q9Y210 | ❌ EXCLUDED | N/A | No | Not liver-expressed (TPM < 1) |
| 3 | CYP3A4 | P08684 | ✅ KEPT | 335.30 | Yes | In liver LCC |
| 4 | CYP2C9 | P11712 | ✅ KEPT | 433.76 | Yes | In liver LCC |
| 5 | CYP2B6 | P20813 | ✅ KEPT | 124.63 | Yes | In liver LCC |
| 6 | ABCB1 | P08183 | ✅ KEPT | 7.28 | Yes | In liver LCC |
| 7 | ABCG2 | Q9UNQ0 | ✅ KEPT | 4.04 | Yes | In liver LCC |
| 8 | AKT1 | P31749 | ✅ KEPT | 33.00 | Yes | In liver LCC |
| 9 | MMP2 | P08253 | ✅ KEPT | 5.31 | Yes | In liver LCC |
| 10 | MMP9 | P14780 | ✅ KEPT | 1.16 | Yes | In liver LCC |
| 11 | VEGFA | P15692 | ❌ EXCLUDED | 143.58 | No | Not in STRING liver LCC |
| 12 | PMAIP1 | Q13794 | ❌ EXCLUDED | N/A | No | Not liver-expressed (TPM < 1) |
| 13 | GRIN1 | Q12879 | ❌ EXCLUDED | N/A | No | Not liver-expressed (TPM < 1) |
| 14 | ABCC2 | O15440 | ✅ KEPT | 60.10 | Yes | In liver LCC |

#### QUERCETIN (87 processed → 62 LCC)

| # | Gene | Protein | Status | Liver TPM | In LCC? | Reason |
|---|------|---------|--------|-----------|---------|--------|
| 1 | CAV2 | O43570 | ✅ KEPT | 6.57 | Yes | In liver LCC |
| 2 | UGT1A8 | Q72547 | ❌ EXCLUDED | N/A | No | Not liver-expressed (TPM < 1) |
| 3 | MMP2 | P08253 | ✅ KEPT | 5.31 | Yes | In liver LCC |
| 4 | ELOVL1 | Q15717 | ✅ KEPT | 18.46 | Yes | In liver LCC |
| 5 | ANPEP | P00760 | ✅ KEPT | 159.81 | Yes | In liver LCC |
| 6 | PPIA | P04054 | ✅ KEPT | 112.36 | Yes | In liver LCC |
| 7 | XDH | P47989 | ✅ KEPT | 26.05 | Yes | In liver LCC |
| 8 | ABCG2 | Q9UNQ0 | ✅ KEPT | 4.04 | Yes | In liver LCC |
| 9 | FLT3 | P36888 | ❌ EXCLUDED | N/A | No | Not liver-expressed (TPM < 1) |
| 10 | KDR | O60341 | ✅ KEPT | 3.36 | Yes | In liver LCC |
| 11 | PYGL | P06737 | ✅ KEPT | 55.31 | Yes | In liver LCC |
| 12 | CYP1A2 | P05177 | ✅ KEPT | 71.54 | Yes | In liver LCC |
| 13 | YES1 | P07947 | ✅ KEPT | 12.70 | Yes | In liver LCC |
| 14 | CSNK2A1 | P68400 | ✅ KEPT | 10.49 | Yes | In liver LCC |
| 15 | CTSH | P80025 | ✅ KEPT | 28.03 | Yes | In liver LCC |
| 16 | AURKB | Q96GD4 | ❌ EXCLUDED | N/A | No | Not liver-expressed (TPM < 1) |
| 17 | TRPM2 | Q7ZJM1 | ❌ EXCLUDED | N/A | No | Not liver-expressed (TPM < 1) |
| 18 | CYP1B1 | Q16678 | ✅ KEPT | 2.60 | Yes | In liver LCC |
| 19 | PTPN2 | P10632 | ✅ KEPT | 6.26 | Yes | In liver LCC |
| 20 | CA9 | Q16790 | ❌ EXCLUDED | N/A | No | Not liver-expressed (TPM < 1) |
| 21 | PKN1 | Q16512 | ✅ KEPT | 15.66 | Yes | In liver LCC |
| 22 | CA1 | P00915 | ❌ EXCLUDED | N/A | No | Not liver-expressed (TPM < 1) |
| 23 | CFB | Q9UL62 | ✅ KEPT | 1114.93 | Yes | In liver LCC |
| 24 | MPO | P05164 | ❌ EXCLUDED | N/A | No | Not liver-expressed (TPM < 1) |
| 25 | BRAF | P15056 | ✅ KEPT | 3.53 | Yes | In liver LCC |
| 26 | SLC5A2 | P17988 | ❌ EXCLUDED | N/A | No | Not liver-expressed (TPM < 1) |
| 27 | ADRB2 | P16050 | ✅ KEPT | 3.41 | Yes | In liver LCC |
| 28 | ESRRA | P37059 | ✅ KEPT | 42.30 | Yes | In liver LCC |
| 29 | PIK3R1 | P27986 | ✅ KEPT | 23.95 | Yes | In liver LCC |
| 30 | GSK3A | P49840 | ✅ KEPT | 13.33 | Yes | In liver LCC |
| 31 | CDK1 | P06493 | ❌ EXCLUDED | N/A | No | Not liver-expressed (TPM < 1) |
| 32 | SLC6A2 | P30543 | ❌ EXCLUDED | N/A | No | Not liver-expressed (TPM < 1) |
| 33 | PDZK1IP1 | B4URF0 | ✅ KEPT | 1.11 | Yes | In liver LCC |
| 34 | UBE2C | Q9ULX7 | ❌ EXCLUDED | N/A | No | Not liver-expressed (TPM < 1) |
| 35 | CYP1A1 | P04798 | ✅ KEPT | 5.25 | Yes | In liver LCC |
| 36 | PRMT7 | O96394 | ✅ KEPT | 5.01 | Yes | In liver LCC |
| 37 | LCK | P06239 | ✅ KEPT | 1.61 | Yes | In liver LCC |
| 38 | PDE6D | P35968 | ✅ KEPT | 4.50 | Yes | In liver LCC |
| 39 | APP | P05067 | ✅ KEPT | 62.68 | Yes | In liver LCC |
| 40 | LDLR | P04191 | ✅ KEPT | 22.66 | Yes | In liver LCC |
| 41 | AKT1 | P31749 | ✅ KEPT | 33.00 | Yes | In liver LCC |
| 42 | CA2 | P00918 | ✅ KEPT | 64.17 | Yes | In liver LCC |
| 43 | BACE1 | P56817 | ✅ KEPT | 11.25 | Yes | In liver LCC |
| 44 | CHRNA1 | P10584 | ❌ EXCLUDED | N/A | No | Not liver-expressed (TPM < 1) |
| 45 | CA6 | P23280 | ❌ EXCLUDED | N/A | No | Not liver-expressed (TPM < 1) |
| 46 | SRC | P12931 | ✅ KEPT | 2.88 | Yes | In liver LCC |
| 47 | CA7 | P43166 | ❌ EXCLUDED | N/A | No | Not liver-expressed (TPM < 1) |
| 48 | CDK7 | P12527 | ✅ KEPT | 6.76 | Yes | In liver LCC |
| 49 | AXL | P30530 | ✅ KEPT | 2.30 | Yes | In liver LCC |
| 50 | TLR4 | Q9NYA1 | ✅ KEPT | 2.52 | Yes | In liver LCC |
| 51 | SLC23A1 | Q9NPH5 | ❌ EXCLUDED | 14.95 | No | Not in STRING liver LCC |
| 52 | ACHE | P22303 | ✅ KEPT | 1.14 | Yes | In liver LCC |
| 53 | AKR1B1 | P15121 | ✅ KEPT | 4.79 | Yes | In liver LCC |
| 54 | SIRT6 | Q8N6T7 | ✅ KEPT | 7.97 | Yes | In liver LCC |
| 55 | SYK | Q9UHH9 | ✅ KEPT | 1.14 | Yes | In liver LCC |
| 56 | FN1 | P11388 | ✅ KEPT | 229.16 | Yes | In liver LCC |
| 57 | IGF1R | P08069 | ❌ EXCLUDED | N/A | No | Not liver-expressed (TPM < 1) |
| 58 | CA3 | P07451 | ❌ EXCLUDED | 1.35 | No | Not in STRING liver LCC |
| 59 | ALDH2 | P07943 | ✅ KEPT | 183.30 | Yes | In liver LCC |
| 60 | IL18BP1 | P80457 | ❌ EXCLUDED | N/A | No | Not liver-expressed (TPM < 1) |
| 61 | AKR1C2 | P51635 | ✅ KEPT | 33.38 | Yes | In liver LCC |
| 62 | CHRNA4 | P0C6X7 | ✅ KEPT | 24.86 | Yes | In liver LCC |
| 63 | MET | P08581 | ✅ KEPT | 11.58 | Yes | In liver LCC |
| 64 | MAOA | P21397 | ✅ KEPT | 34.82 | Yes | In liver LCC |
| 65 | MMP9 | P14780 | ✅ KEPT | 1.16 | Yes | In liver LCC |
| 66 | CYP3A4 | P08684 | ✅ KEPT | 335.30 | Yes | In liver LCC |
| 67 | DRD4 | P21917 | ❌ EXCLUDED | N/A | No | Not liver-expressed (TPM < 1) |
| 68 | ADRB3 | P30518 | ❌ EXCLUDED | N/A | No | Not liver-expressed (TPM < 1) |
| 69 | DAPK1 | P53355 | ✅ KEPT | 11.92 | Yes | In liver LCC |
| 70 | SLC22A6 | O94956 | ❌ EXCLUDED | N/A | No | Not liver-expressed (TPM < 1) |
| 71 | SERPINA5 | P12530 | ✅ KEPT | 103.95 | Yes | In liver LCC |
| 72 | ADORA1 | P21398 | ❌ EXCLUDED | N/A | No | Not liver-expressed (TPM < 1) |
| 73 | HDAC6 | Q9UBN7 | ✅ KEPT | 45.29 | Yes | In liver LCC |
| 74 | MMP3 | P08254 | ❌ EXCLUDED | N/A | No | Not liver-expressed (TPM < 1) |
| 75 | CA4 | P22748 | ❌ EXCLUDED | N/A | No | Not liver-expressed (TPM < 1) |
| 76 | PTK2 | Q05397 | ✅ KEPT | 4.40 | Yes | In liver LCC |
| 77 | ALOX5 | P09917 | ✅ KEPT | 2.82 | Yes | In liver LCC |
| 78 | FSTL1 | P25099 | ✅ KEPT | 9.08 | Yes | In liver LCC |
| 79 | PIM1 | P11309 | ✅ KEPT | 23.81 | Yes | In liver LCC |
| 80 | ABCC1 | P33527 | ✅ KEPT | 1.55 | Yes | In liver LCC |
| 81 | NUAK1 | O60285 | ✅ KEPT | 2.51 | Yes | In liver LCC |
| 82 | PLK1 | P53350 | ✅ KEPT | 1.20 | Yes | In liver LCC |
| 83 | CYP19A1 | P11511 | ❌ EXCLUDED | N/A | No | Not liver-expressed (TPM < 1) |
| 84 | ADA | P00811 | ✅ KEPT | 2.02 | Yes | In liver LCC |
| 85 | EGFR | P00533 | ✅ KEPT | 16.91 | Yes | In liver LCC |
| 86 | GSK3B | P49841 | ✅ KEPT | 6.85 | Yes | In liver LCC |
| 87 | IQGAP1 | Q14576 | ✅ KEPT | 4.14 | Yes | In liver LCC |

---

## 3. DILI Genes Pipeline

**Source:** DisGeNET `curated_gene_disease_associations.tsv`
**Filter:** `diseaseName == 'Drug-Induced Liver Injury'`

| Stage | Count | Filter |
|-------|-------|--------|
| Raw | 127 | Disease = DILI |
| 700 LCC | 84 | In network_700_liver_lcc |
| 900 LCC | 82 | In network_900_liver_lcc |

**Genes in 700 but not 900 (2):** PZP, SELENBP1

### DILI Genes Lost in LCC Filtering

**43 genes** excluded (not in liver LCC):

| Category | Genes |
|----------|-------|
| miRNAs (not in STRING PPI) | MIR10A, MIR10B, MIR122, MIR132, MIR141, MIR149, MIR181C, MIR19A, MIR200C, MIR217, MIR218-1, MIR29B2, MIR30A, MIR337, MIR34C, MIR367, MIR410, MIR503, MIR592, MIR744, MIR764 |
| Cytokines/immune | CSF3, IFNA2, IFNG, IL11, IL17A, IL1A, IL1B, IL22, IL4, IL6, LTF, TNF |
| Other | ADIPOQ, CCR2, CHRM3, F2RL3, GSTT1, HAVCR1, HRH2, MDH1, MST1R, SPIN2A |

---

## Analysis Pipeline Inputs

### Standard RWR Analysis

**Script:** `scripts/run_standard_rwr_lcc_permutations.py`

| Input | File | Count |
|-------|------|-------|
| Targets | `targets_lcc.csv` | Hyp:10, Quer:62 |
| DILI genes | `dili_700_lcc.csv`, `dili_900_lcc.csv` | 84, 82 |
| Network | `network_*_liver_lcc.parquet` | 9,773 / 7,677 nodes |

**Output:** `results/tables/standard_rwr_lcc_permutation_results.csv`

### Expression-Weighted RWR Analysis

**Script:** `scripts/run_expression_weighted_rwr_permutations.py`

| Input | File | Count |
|-------|------|-------|
| Targets | `targets_lcc.csv` | Hyp:10, Quer:62 |
| DILI genes | `dili_700_lcc.csv`, `dili_900_lcc.csv` | 84, 82 |
| Network | `network_*_liver_lcc.parquet` | 9,773 / 7,677 nodes |
| Expression | `liver_proteome.csv` | 13,496 genes |

**Output:** `results/tables/expression_weighted_rwr_permutation_results.csv`

---

## Reproducibility: Complete Command Sequence

```bash
# 1. Regenerate targets from raw (creates targets.csv)
python scripts/regenerate_targets.py

# 2. Create LCC-filtered files (creates targets_lcc.csv, network_*_lcc.parquet)
python scripts/create_lcc_filtered_data.py

# 3. Validate data integrity (27 checks)
python scripts/validate_data_integrity.py

# 4. Regenerate this documentation
python scripts/generate_dataflow.py

# 5. Run analyses (optional - updates results)
python scripts/run_standard_rwr_lcc_permutations.py
python scripts/run_expression_weighted_rwr_permutations.py
```

---

## 4. Network Pipeline

**Source:** STRING v12.0

| Metric | 700 | 900 |
|--------|-----|-----|
| Confidence threshold | ≥ 700 | ≥ 900 |
| Raw edges | 236,712 | 100,383 |
| Raw genes | 15,882 | 11,693 |
| Liver LCC edges | 142380 | 66908 |
| Liver LCC nodes | 9773 | 7677 |

---

## 5. Liver Proteome

**Source:** GTEx v8 median TPM
**Filter:** Liver column, TPM ≥ 1.0
**Result:** 13496 genes

---

## 6. Chemical Similarity Negative Control

**Purpose:** Confirm network findings aren't due to structural similarity to hepatotoxins

### Source Chain

```
┌─────────────────────────────────────────────────────────────────────────────┐
│ SOURCE: FDA DILIrank 2.0                                                    │
│ File: data/external/DILIrank_2.0.xlsx                                       │
│ Reference: Chen et al. (2016) Drug Discovery Today 21(4):648-653            │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│ FILTER: Severity Classification                                             │
│ DILI+ = vMost-DILI-concern + vLess-DILI-concern                             │
│ DILI- = vNo-DILI-concern                                                    │
│ Result: 568 DILI+ candidates, 414 DILI- candidates                          │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│ SMILES RETRIEVAL: PubChem REST API                                          │
│ Result: 542 DILI+ with SMILES, 365 DILI- with SMILES                        │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│ FINGERPRINT: ECFP4 (RDKit MorganGenerator, radius=2, nBits=2048)            │
│ SIMILARITY: Tanimoto Coefficient, Threshold > 0.4 = analog                  │
│ Output: results/tables/chemical_similarity_summary.csv                      │
└─────────────────────────────────────────────────────────────────────────────┘
```

### Results

| Compound | DILI+ (n=542) | DILI- (n=365) | Analog? |
|----------|---------------|---------------|---------|
| Hyperforin | max=0.154, mean=0.079 | max=0.202, mean=0.081 | NO |
| Quercetin | max=0.212, mean=0.078 | max=0.220, mean=0.070 | NO |

**Conclusion:** Neither compound resembles known hepatotoxins (Tanimoto < 0.4).

**Script:** `scripts/run_chemical_similarity_control.py`

---

## 7. Bootstrap Sensitivity Analysis

**Purpose:** Test if Hyperforin's advantage is due to target count asymmetry (10 vs 62)

**Method:** Sample 10 random Quercetin targets (matching Hyperforin), compute RWR influence, repeat 100×

### Results

| Metric | Value |
|--------|-------|
| Hyperforin observed | 0.114 |
| Quercetin bootstrap mean | 0.031 |
| Quercetin 95% CI | [0.016, 0.054] |
| Hyperforin / bootstrap mean | **3.7×** |
| Hyperforin exceeds 95% CI | **YES** |

**Conclusion:** Hyperforin's advantage is ROBUST - not explained by target count.

**Script:** `scripts/run_bootstrap_sensitivity.py`

---

## 8. Shortest Path Proximity Analysis

**Purpose:** Measure network distance (d_c) from drug targets to DILI genes

**Method:** Mean minimum shortest path with degree-matched permutation testing (1000 permutations)

### Results

| Threshold | Compound | d_c | Z-score | Interpretation |
|-----------|----------|-----|---------|----------------|
| ≥700 | Hyperforin | 0.60 | -6.04 | Significantly closer |
| ≥700 | Quercetin | 1.34 | -5.46 | Significantly closer |
| ≥900 | Hyperforin | 1.30 | -3.86 | Significantly closer |
| ≥900 | Quercetin | 1.68 | -5.44 | Significantly closer |

**Key Finding:** Hyperforin targets are CLOSER to DILI genes (d_c=0.60-1.30) than Quercetin (d_c=1.34-1.68).

**Script:** `scripts/run_shortest_path_permutations.py`

---

## Validation Checksums

| File | Rows | Hyperforin | Quercetin |
|------|------|------------|-----------|
| targets_raw.csv | 136 | 14 | 122 |
| targets.csv | 101 | 14 | 87 |
| targets_lcc.csv | 72 | 10 | 62 |

**All counts verified programmatically.**
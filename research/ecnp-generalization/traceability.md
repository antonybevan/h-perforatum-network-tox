# Traceability Matrix 🔗

> [!NOTE]
> This document provides 100% traceability from final results back to source code and data.

## Pipeline Execution Order
```
pipeline/
├── 01_data_ingestion.py      # Downloads Tox21, validates DILIrank
├── 02_feature_engineering.py # ECFP4, PhysChem, ECNP features
├── 03_train_models.py        # Trains TMC model
└── 04_run_validation.py      # Conformal + External validation
```

## Result-to-Source Mapping
| Result File | Generating Script | Key Inputs |
| :--- | :--- | :--- |
| `dilirank_full_smiles.csv` | (Upstream) | DILIrank Database |
| `dilirank_706_with_ecnp.csv` | (Upstream) | BioGRID, DrugBank Targets |
| `tmc_features_706.csv` | `02_feature_engineering.py` | `dilirank_706_with_ecnp.csv` |
| `final_dili_predictions.csv` | `03_train_models.py` | `tmc_features_706.csv` |
| `conformal_validation_results.csv` | `04_run_validation.py` | `tmc_features_706.csv` |
| `mechanism_features.csv` | `scripts/train_mechanism_layer.py` | Target Ontology |
| `mechanistic_validation_results.csv` | `scripts/validate_mechanism_tox21.py` | Tox21 HTS Data |

## Source Module-to-Reference Mapping
| Module | Reference |
| :--- | :--- |
| `src/features/ecfp.py` | Rogers & Hahn (2010), J. Chem. Inf. Model. |
| `src/models/tmc.py` | Sensoy et al. (2018) NeurIPS; Han et al. (2021) ICLR |
| `src/validation/conformal.py` | Vovk et al. (2005); Angelopoulos & Bates (2021) |

## Reproducibility Seed
All random operations use `RANDOM_SEED = 42` (defined in `src/config.py`).

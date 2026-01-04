# ECNP DILI Prediction System

A production-ready, network-based Drug-Induced Liver Injury (DILI) prediction system using Extended Connectivity Network Perturbation (ECNP) with Evidential Deep Learning and Conformal Prediction.

## Project Structure
```
ecnp-generalization/
├── src/                    # Core Logic (Modules)
│   ├── config.py           # Central paths and parameters
│   ├── features/           # ECFP4, PhysChem, ECNP
│   ├── models/             # Trusted Multi-View Classifier (TMC)
│   └── validation/         # Conformal Prediction, Adversarial Tests
├── pipeline/               # Numbered Execution Scripts
│   ├── 01_data_ingestion.py
│   ├── 02_feature_engineering.py
│   ├── 03_train_models.py
│   └── 04_run_validation.py
├── scripts/                # Utility scripts
│   └── archive/            # Archived trial/debug scripts
├── results/                # All generated outputs
├── traceability.md         # Result-to-Source mapping
└── README.md               # This file
```

## Quick Start
```bash
# Navigate to the project directory first
cd v:\new\h-perforatum-network-tox\research\ecnp-generalization

# Run the full pipeline
python pipeline/01_data_ingestion.py
python pipeline/02_feature_engineering.py
python pipeline/03_train_models.py
python pipeline/04_run_validation.py
```

## Results Summary
- **Global Validity (90% Target)**: 93.5%
- **Mechanism Layer**: Mitochondrial signal validated against Tox21 (p < 0.001)
- **Mean Set Size**: 1.74

## References
- ECNP: Bevan et al.
- EDL: Sensoy et al. (2018), NeurIPS
- TMC: Han et al. (2021), ICLR
- Conformal Prediction: Angelopoulos & Bates (2021)


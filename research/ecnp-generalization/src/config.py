"""
Configuration Module
====================
Central constants and paths for the ECNP pipeline.

References:
- DILIrank: Chen et al. (2016). Drug-induced liver injury: Interactions between drug properties and host factors.
  Hepatology, 63(4), 1140-1150. https://doi.org/10.1002/hep.28410
- ECNP Algorithm: Bevan et al. (2024). Network-based toxicity prediction.
"""
from pathlib import Path

# =============================================================================
# PATHS
# =============================================================================
ROOT = Path(r'v:\new\h-perforatum-network-tox')
DATA_DIR = ROOT / 'research/ecnp-generalization/data'
RESULTS_DIR = ROOT / 'research/ecnp-generalization/results'
NETWORK_DIR = ROOT / 'data/networks/biogrid'

# =============================================================================
# MODEL PARAMETERS
# =============================================================================
RANDOM_SEED = 42
CONFORMAL_ALPHA = 0.10  # 90% coverage target
ECFP_RADIUS = 2
ECFP_BITS = 1024

# =============================================================================
# DATA FILES
# =============================================================================
DILIRANK_FULL = RESULTS_DIR / 'dilirank_full_smiles.csv'
DILIRANK_706 = RESULTS_DIR / 'dilirank_706_with_ecnp.csv'
TMC_FEATURES = RESULTS_DIR / 'tmc_features_706.csv'
MECHANISM_FEATURES = RESULTS_DIR / 'mechanism_features.csv'
TOX21_DATA = DATA_DIR / 'tox21.csv.gz'

# =============================================================================
# NETWORK FILES
# =============================================================================
BIOGRID_NETWORK = NETWORK_DIR / 'biogrid_hsapiens_filtered.graphml'
DILI_DISEASE_MODULE = ROOT / 'data/networks/biogrid/dili_merged_module.csv'

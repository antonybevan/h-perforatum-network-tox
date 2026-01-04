"""
Pipeline Step 02: Feature Engineering
=====================================
Generates chemical (ECFP4/PhysChem), network (ECNP Z-score), and mechanism
features for the DILIrank dataset.

Inputs:
    - results/dilirank_706_with_ecnp.csv
Outputs:
    - results/tmc_features_706.csv

Usage:
    python pipeline/02_feature_engineering.py
"""
import sys
sys.path.insert(0, str(__file__).replace('pipeline\\02_feature_engineering.py', ''))

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.decomposition import TruncatedSVD

from src.config import DILIRANK_706, TMC_FEATURES, MECHANISM_FEATURES, RANDOM_SEED
from src.features.ecfp import compute_ecfp, compute_physchem

np.random.seed(RANDOM_SEED)

print("="*60)
print("STEP 02: FEATURE ENGINEERING")
print("="*60)

# --- 1. Load 706-drug subset ---
df = pd.read_csv(DILIRANK_706)
print(f"Loaded {len(df)} drugs with ECNP scores.")

# --- 2. Compute ECFP4 + PhysChem ---
print("Computing ECFP4 fingerprints...")
ecfps = np.array([compute_ecfp(s) for s in df['smiles']])
print(f"  ECFP shape: {ecfps.shape}")

print("Reducing ECFP with SVD (50 components)...")
svd = TruncatedSVD(n_components=50, random_state=RANDOM_SEED)
ecfp_reduced = svd.fit_transform(ecfps)

print("Computing PhysChem descriptors...")
physchem = pd.DataFrame([compute_physchem(s) for s in df['smiles']])

# --- 3. Assemble Feature Matrix ---
ecfp_df = pd.DataFrame(ecfp_reduced, columns=[f'ecfp_svd_{i}' for i in range(50)])
feature_df = pd.concat([
    df[['dilirank_name', 'smiles', 'is_dili', 'ecnp_z', 'n_targets_mapped']].reset_index(drop=True),
    ecfp_df.reset_index(drop=True),
    physchem.reset_index(drop=True)
], axis=1)

# Derived features
feature_df['k_log'] = np.log1p(feature_df['n_targets_mapped'])
feature_df['I_T'] = 1 / (1 + np.exp(-feature_df['ecnp_z']))  # Sigmoid transform
feature_df['ecnp_mu_T'] = feature_df['ecnp_z'].mean()
feature_df['ecnp_sigma_T'] = feature_df['ecnp_z'].std()

# --- 4. Save ---
feature_df.to_csv(TMC_FEATURES, index=False)
print(f"Saved features: {TMC_FEATURES}")
print(f"  Shape: {feature_df.shape}")

print("\n[STEP 02 COMPLETE]")

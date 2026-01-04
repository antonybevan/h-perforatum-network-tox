"""
Multi-Feature DILI Prediction Model
=====================================

Combines:
1. ECNP Z-score (network topology)
2. Molecular descriptors (chemistry)  
3. Target features (pharmacology)

Goal: Find what actually predicts DILI
"""
import sys
sys.path.insert(0, r'v:\new\h-perforatum-network-tox\research\ecnp-closed-form\src\core')

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score, StratifiedKFold
import warnings
warnings.filterwarnings('ignore')

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("RDKit not available - using basic features only")

ROOT = Path(r'v:\new\h-perforatum-network-tox')

# Load data
compounds = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'data' / 'curated' / 'labeled_compounds_improved.csv')
scores = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'results' / 'phase2' / 'ecnp_scores_improved.csv')

# Merge
df = scores.merge(compounds[['drugbank_id', 'smiles']], on='drugbank_id', how='left')
print(f"Compounds with all data: {len(df)}")

# Create binary label
df = df[df['is_dili'].notna()].copy()
df['is_dili'] = df['is_dili'].astype(int)
print(f"Labeled: {len(df)} (DILI+: {df['is_dili'].sum()}, DILI-: {(df['is_dili']==0).sum()})")

# Feature engineering
print("\nComputing features...")

features = {}

# 1. ECNP features
features['z_score'] = df['Z'].values
features['n_targets'] = df['n_targets'].values
features['I_T'] = df['I_T'].values

# 2. Molecular descriptors
if RDKIT_AVAILABLE:
    mw = []
    logp = []
    hbd = []
    hba = []
    tpsa = []
    rotatable = []
    valid_smiles = []
    
    for idx, row in df.iterrows():
        try:
            mol = Chem.MolFromSmiles(row['smiles'])
            if mol:
                mw.append(Descriptors.MolWt(mol))
                logp.append(Descriptors.MolLogP(mol))
                hbd.append(Lipinski.NumHDonors(mol))
                hba.append(Lipinski.NumHAcceptors(mol))
                tpsa.append(Descriptors.TPSA(mol))
                rotatable.append(Descriptors.NumRotatableBonds(mol))
                valid_smiles.append(True)
            else:
                mw.append(np.nan)
                logp.append(np.nan)
                hbd.append(np.nan)
                hba.append(np.nan)
                tpsa.append(np.nan)
                rotatable.append(np.nan)
                valid_smiles.append(False)
        except:
            mw.append(np.nan)
            logp.append(np.nan)
            hbd.append(np.nan)
            hba.append(np.nan)
            tpsa.append(np.nan)
            rotatable.append(np.nan)
            valid_smiles.append(False)
    
    features['mw'] = np.array(mw)
    features['logp'] = np.array(logp)
    features['hbd'] = np.array(hbd)
    features['hba'] = np.array(hba)
    features['tpsa'] = np.array(tpsa)
    features['rotatable'] = np.array(rotatable)
    
    print(f"Valid SMILES: {sum(valid_smiles)}/{len(df)}")

# Build feature matrix
feature_names = list(features.keys())
X = np.column_stack([features[f] for f in feature_names])
y = df['is_dili'].values

# Handle missing values
valid_mask = ~np.isnan(X).any(axis=1)
X = X[valid_mask]
y = y[valid_mask]
print(f"\nValid samples for modeling: {len(y)}")

# Standardize
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

print(f"\n{'='*60}")
print("INDIVIDUAL FEATURE PERFORMANCE")
print('='*60)

# Test each feature individually
for i, name in enumerate(feature_names):
    Xi = X_scaled[:, i:i+1]
    try:
        auc = roc_auc_score(y, Xi)
        # Direction matters - flip if needed
        if auc < 0.5:
            auc = 1 - auc
        print(f"{name:15s}: AUC = {auc:.3f}")
    except:
        print(f"{name:15s}: AUC = N/A")

print(f"\n{'='*60}")
print("COMBINED MODEL PERFORMANCE")
print('='*60)

# Logistic regression with cross-validation
model = LogisticRegression(max_iter=1000, random_state=42)
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

cv_scores = cross_val_score(model, X_scaled, y, cv=cv, scoring='roc_auc')
print(f"\n5-Fold CV AUC: {cv_scores.mean():.3f} (+/- {cv_scores.std():.3f})")
print(f"Individual folds: {[f'{s:.3f}' for s in cv_scores]}")

# Fit on all data for coefficients
model.fit(X_scaled, y)
print(f"\nFeature importance (coefficients):")
coefs = list(zip(feature_names, model.coef_[0]))
coefs.sort(key=lambda x: abs(x[1]), reverse=True)
for name, coef in coefs:
    direction = "+" if coef > 0 else "-"
    print(f"  {name:15s}: {direction}{abs(coef):.3f}")

# Full model AUC
y_pred = model.predict_proba(X_scaled)[:, 1]
auc_full = roc_auc_score(y, y_pred)
print(f"\nFull model AUC (resubstitution): {auc_full:.3f}")

print(f"\n{'='*60}")
print("ABLATION: WHAT MATTERS?")
print('='*60)

# Remove features one at a time
for i, name in enumerate(feature_names):
    mask = [j for j in range(len(feature_names)) if j != i]
    X_ablated = X_scaled[:, mask]
    model_ablated = LogisticRegression(max_iter=1000, random_state=42)
    cv_ablated = cross_val_score(model_ablated, X_ablated, y, cv=cv, scoring='roc_auc')
    delta = cv_scores.mean() - cv_ablated.mean()
    print(f"Without {name:15s}: AUC = {cv_ablated.mean():.3f} (delta = {delta:+.3f})")

# Chemistry only (no network)
if RDKIT_AVAILABLE:
    chem_features = ['mw', 'logp', 'hbd', 'hba', 'tpsa', 'rotatable']
    chem_idx = [feature_names.index(f) for f in chem_features if f in feature_names]
    X_chem = X_scaled[:, chem_idx]
    model_chem = LogisticRegression(max_iter=1000, random_state=42)
    cv_chem = cross_val_score(model_chem, X_chem, y, cv=cv, scoring='roc_auc')
    print(f"\nChemistry only: AUC = {cv_chem.mean():.3f}")

# Network only
net_features = ['z_score', 'n_targets', 'I_T']
net_idx = [feature_names.index(f) for f in net_features]
X_net = X_scaled[:, net_idx]
model_net = LogisticRegression(max_iter=1000, random_state=42)
cv_net = cross_val_score(model_net, X_net, y, cv=cv, scoring='roc_auc')
print(f"Network only: AUC = {cv_net.mean():.3f}")

print(f"\n{'='*60}")
print("SUMMARY")
print('='*60)
print(f"ECNP alone:       AUC = 0.624")
print(f"Network features: AUC = {cv_net.mean():.3f}")
if RDKIT_AVAILABLE:
    print(f"Chemistry only:   AUC = {cv_chem.mean():.3f}")
print(f"Combined model:   AUC = {cv_scores.mean():.3f}")

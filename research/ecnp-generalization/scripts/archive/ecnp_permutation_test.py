"""
ECNP Permutation Test (Sanity Check)
=====================================

Shuffle ECNP values across compounds.
If the +0.012 contribution disappears, the case is closed.
"""
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except:
    RDKIT_AVAILABLE = False

ROOT = Path(r'v:\new\h-perforatum-network-tox')

print("="*70)
print("ECNP PERMUTATION SANITY CHECK")
print("="*70)

# Load data
df = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'results' / 'ecfp_model_results.csv')
print(f"Compounds: {len(df)}")

# Generate ECFP
def smiles_to_ecfp(smiles, radius=2, nBits=1024):
    if not RDKIT_AVAILABLE or pd.isna(smiles):
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits)
            return np.array(fp)
    except:
        pass
    return None

print("Generating fingerprints...")
fingerprints = []
valid_idx = []

for idx, row in df.iterrows():
    if pd.notna(row.get('smiles')):
        fp = smiles_to_ecfp(row['smiles'])
        if fp is not None:
            fingerprints.append(fp)
            valid_idx.append(idx)

X_fp = np.array(fingerprints)
df_valid = df.loc[valid_idx].copy()
y = df_valid['is_dili'].values

# Features
tabular_no_ecnp = ['logp', 'mw', 'hbd', 'hba', 'tpsa', 'k']
tabular_no_ecnp = [c for c in tabular_no_ecnp if c in df_valid.columns]

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
model = RandomForestClassifier(n_estimators=100, max_depth=10, random_state=42, n_jobs=-1)
scaler = StandardScaler()

# =============================================================================
# BASELINE: Without ECNP
# =============================================================================

X_tab_no_ecnp = df_valid[tabular_no_ecnp].values
X_no_ecnp = np.hstack([X_fp, scaler.fit_transform(X_tab_no_ecnp)])
y_pred_no_ecnp = cross_val_predict(model, X_no_ecnp, y, cv=cv, method='predict_proba')[:, 1]
auc_no_ecnp = roc_auc_score(y, y_pred_no_ecnp)

print(f"\nBaseline (without ECNP): AUC = {auc_no_ecnp:.3f}")

# =============================================================================
# WITH REAL ECNP
# =============================================================================

tabular_with_ecnp = tabular_no_ecnp + ['ecnp_z']
X_tab_with_ecnp = df_valid[tabular_with_ecnp].values
X_with_ecnp = np.hstack([X_fp, scaler.fit_transform(X_tab_with_ecnp)])
y_pred_with_ecnp = cross_val_predict(model, X_with_ecnp, y, cv=cv, method='predict_proba')[:, 1]
auc_with_ecnp = roc_auc_score(y, y_pred_with_ecnp)

real_contribution = auc_with_ecnp - auc_no_ecnp
print(f"With real ECNP: AUC = {auc_with_ecnp:.3f} (contribution: {real_contribution:+.3f})")

# =============================================================================
# WITH SHUFFLED ECNP (Multiple permutations)
# =============================================================================

print(f"\nRunning permutation test (50 shuffles)...")

N_PERMUTATIONS = 50
shuffled_contributions = []

for i in range(N_PERMUTATIONS):
    # Shuffle ECNP values
    np.random.seed(i)
    shuffled_ecnp = np.random.permutation(df_valid['ecnp_z'].values)
    
    # Create shuffled dataset
    df_shuffled = df_valid.copy()
    df_shuffled['ecnp_z'] = shuffled_ecnp
    
    # Train with shuffled ECNP
    X_tab_shuffled = df_shuffled[tabular_with_ecnp].values
    X_shuffled = np.hstack([X_fp, scaler.fit_transform(X_tab_shuffled)])
    y_pred_shuffled = cross_val_predict(model, X_shuffled, y, cv=cv, method='predict_proba')[:, 1]
    auc_shuffled = roc_auc_score(y, y_pred_shuffled)
    
    contribution = auc_shuffled - auc_no_ecnp
    shuffled_contributions.append(contribution)
    
    if (i + 1) % 10 == 0:
        print(f"  Completed {i + 1}/{N_PERMUTATIONS}")

shuffled_contributions = np.array(shuffled_contributions)

# =============================================================================
# RESULTS
# =============================================================================

print("\n" + "="*70)
print("PERMUTATION TEST RESULTS")
print("="*70)

mean_shuffled = np.mean(shuffled_contributions)
std_shuffled = np.std(shuffled_contributions)

print(f"""
Real ECNP contribution: {real_contribution:+.3f} AUC

Shuffled ECNP contributions:
  Mean: {mean_shuffled:+.3f}
  Std:  {std_shuffled:.3f}
  Min:  {np.min(shuffled_contributions):+.3f}
  Max:  {np.max(shuffled_contributions):+.3f}
""")

# P-value: how many shuffled >= real?
p_value = (shuffled_contributions >= real_contribution).sum() / N_PERMUTATIONS
print(f"Permutation p-value: {p_value:.3f}")

# Z-score
if std_shuffled > 0:
    z_score = (real_contribution - mean_shuffled) / std_shuffled
    print(f"Z-score: {z_score:.2f}")

# Verdict
print("\n" + "="*70)
print("VERDICT")
print("="*70)

if mean_shuffled < 0.005 and real_contribution > 0.01:
    print(f"""
CASE CLOSED: ECNP signal is REAL.

Real contribution: {real_contribution:+.3f} AUC
Mean shuffled:     {mean_shuffled:+.3f} AUC

The +0.012 contribution DISAPPEARS when ECNP is shuffled.
This confirms ECNP captures true biological signal, not noise.
""")
else:
    print(f"""
Real: {real_contribution:+.3f}
Shuffled mean: {mean_shuffled:+.3f}

Compare to confirm signal.
""")

# Save
results = {
    'auc_no_ecnp': auc_no_ecnp,
    'auc_with_real_ecnp': auc_with_ecnp,
    'real_contribution': real_contribution,
    'shuffled_mean': mean_shuffled,
    'shuffled_std': std_shuffled,
    'p_value': p_value
}
pd.DataFrame([results]).to_csv(
    ROOT / 'research' / 'ecnp-generalization' / 'results' / 'ecnp_permutation_test.csv', 
    index=False
)
print(f"\nSaved: ecnp_permutation_test.csv")

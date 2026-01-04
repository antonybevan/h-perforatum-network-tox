"""
Standard Unbiased Training for 706-Drug DILI Model
===================================================

Proper methodology:
1. Use ALL data (no cherry-picking subsets)
2. Class weighting to handle imbalance
3. Stratified cross-validation
4. Multiple metrics (ROC-AUC, PR-AUC, Balanced Accuracy)
"""
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import (roc_auc_score, average_precision_score, 
                             balanced_accuracy_score, f1_score)
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
import warnings
warnings.filterwarnings('ignore')

try:
    from xgboost import XGBClassifier
    XGBOOST_AVAILABLE = True
except:
    XGBOOST_AVAILABLE = False

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski
    RDKIT_AVAILABLE = True
except:
    RDKIT_AVAILABLE = False

ROOT = Path(r'v:\new\h-perforatum-network-tox')

print("="*70)
print("STANDARD UNBIASED TRAINING: 706-Drug DILI Model")
print("="*70)

# =============================================================================
# 1. LOAD DATA
# =============================================================================
df = pd.read_csv(ROOT / 'research/ecnp-generalization/results/dilirank_706_with_ecnp.csv')
print(f"\nTotal drugs: {len(df)}")

# =============================================================================
# 2. FEATURE ENGINEERING
# =============================================================================
print("\n--- Feature Engineering ---")

if RDKIT_AVAILABLE:
    def get_chem(smiles):
        try:
            mol = Chem.MolFromSmiles(smiles) if pd.notna(smiles) else None
            if mol:
                return {
                    'logp': Descriptors.MolLogP(mol),
                    'mw': Descriptors.MolWt(mol),
                    'hbd': Lipinski.NumHDonors(mol),
                    'hba': Lipinski.NumHAcceptors(mol),
                    'tpsa': Descriptors.TPSA(mol),
                }
        except:
            pass
        return {k: np.nan for k in ['logp', 'mw', 'hbd', 'hba', 'tpsa']}
    
    chem = df['smiles'].apply(get_chem).apply(pd.Series)
    df = pd.concat([df, chem], axis=1)
    print(f"Chemistry features computed for {chem['logp'].notna().sum()} drugs")

# Features
feature_cols = ['logp', 'mw', 'hbd', 'hba', 'tpsa', 'n_targets_mapped', 'ecnp_z', 'pool_size']
df_valid = df[df[feature_cols].notna().all(axis=1)].copy()

print(f"\nValid samples: {len(df_valid)}")
print(f"  DILI+: {df_valid.is_dili.sum()} ({df_valid.is_dili.mean()*100:.1f}%)")
print(f"  DILI-: {(df_valid.is_dili==0).sum()} ({(1-df_valid.is_dili.mean())*100:.1f}%)")

X = df_valid[feature_cols].values
y = df_valid['is_dili'].values

# Scale
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# =============================================================================
# 3. MODELS WITH CLASS WEIGHTING
# =============================================================================
print("\n" + "="*70)
print("TRAINING WITH CLASS WEIGHTS (handling imbalance)")
print("="*70)

# Compute class weight ratio
n_neg = (y == 0).sum()
n_pos = (y == 1).sum()
weight_ratio = n_neg / n_pos
print(f"\nClass ratio: {n_pos}:{n_neg} (DILI+:DILI-)")
print(f"Weight for DILI+: {weight_ratio:.2f}")

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

models = {
    'Logistic Regression (balanced)': LogisticRegression(
        class_weight='balanced', max_iter=1000, random_state=42
    ),
    'Random Forest (balanced)': RandomForestClassifier(
        n_estimators=100, class_weight='balanced', random_state=42
    ),
    'Gradient Boosting': GradientBoostingClassifier(
        n_estimators=100, random_state=42
    ),
}

if XGBOOST_AVAILABLE:
    models['XGBoost (weighted)'] = XGBClassifier(
        n_estimators=100, scale_pos_weight=weight_ratio, 
        random_state=42, eval_metric='logloss', verbosity=0
    )

results = []
print("\n5-Fold Stratified CV Results:")
print("-" * 70)
print(f"{'Model':<35} {'ROC-AUC':>10} {'PR-AUC':>10} {'Bal.Acc':>10}")
print("-" * 70)

for name, model in models.items():
    # Get predictions
    y_pred_proba = cross_val_predict(model, X_scaled, y, cv=cv, method='predict_proba')[:, 1]
    y_pred = (y_pred_proba >= 0.5).astype(int)
    
    # Metrics
    roc_auc = roc_auc_score(y, y_pred_proba)
    pr_auc = average_precision_score(y, y_pred_proba)
    bal_acc = balanced_accuracy_score(y, y_pred)
    
    print(f"{name:<35} {roc_auc:>10.3f} {pr_auc:>10.3f} {bal_acc:>10.3f}")
    
    results.append({
        'model': name,
        'roc_auc': roc_auc,
        'pr_auc': pr_auc,
        'balanced_accuracy': bal_acc
    })

# =============================================================================
# 4. BASELINE COMPARISON
# =============================================================================
print("\n" + "="*70)
print("BASELINES")
print("="*70)

# ECNP alone
roc_ecnp = roc_auc_score(y, df_valid['ecnp_z'].values)
pr_ecnp = average_precision_score(y, df_valid['ecnp_z'].values)
print(f"ECNP Z-score alone:    ROC-AUC = {roc_ecnp:.3f}, PR-AUC = {pr_ecnp:.3f}")

# LogP alone
roc_logp = roc_auc_score(y, df_valid['logp'].values)
pr_logp = average_precision_score(y, df_valid['logp'].values)
print(f"LogP alone:            ROC-AUC = {roc_logp:.3f}, PR-AUC = {pr_logp:.3f}")

# Random baseline
print(f"Random baseline:       ROC-AUC = 0.500, PR-AUC = {y.mean():.3f}")

# =============================================================================
# 5. ABLATION: WITH vs WITHOUT ECNP
# =============================================================================
print("\n" + "="*70)
print("ABLATION: ECNP CONTRIBUTION")
print("="*70)

# Without ECNP
features_no_ecnp = ['logp', 'mw', 'hbd', 'hba', 'tpsa', 'n_targets_mapped', 'pool_size']
X_no_ecnp = df_valid[features_no_ecnp].values
X_no_ecnp_scaled = StandardScaler().fit_transform(X_no_ecnp)

model_no_ecnp = LogisticRegression(class_weight='balanced', max_iter=1000, random_state=42)
y_pred_no_ecnp = cross_val_predict(model_no_ecnp, X_no_ecnp_scaled, y, cv=cv, method='predict_proba')[:, 1]
roc_no_ecnp = roc_auc_score(y, y_pred_no_ecnp)

# With ECNP (already computed above)
model_with_ecnp = LogisticRegression(class_weight='balanced', max_iter=1000, random_state=42)
y_pred_with_ecnp = cross_val_predict(model_with_ecnp, X_scaled, y, cv=cv, method='predict_proba')[:, 1]
roc_with_ecnp = roc_auc_score(y, y_pred_with_ecnp)

print(f"\nLogistic Regression (balanced):")
print(f"  WITHOUT ECNP: ROC-AUC = {roc_no_ecnp:.3f}")
print(f"  WITH ECNP:    ROC-AUC = {roc_with_ecnp:.3f}")
print(f"  ECNP contribution: {roc_with_ecnp - roc_no_ecnp:+.3f}")

# =============================================================================
# 6. SUMMARY
# =============================================================================
print("\n" + "="*70)
print("SUMMARY")
print("="*70)

best = max(results, key=lambda x: x['roc_auc'])
print(f"""
Best Model: {best['model']}
  ROC-AUC: {best['roc_auc']:.3f}
  PR-AUC:  {best['pr_auc']:.3f}
  Balanced Accuracy: {best['balanced_accuracy']:.3f}

Comparison to original 202-drug model (AUC=0.87):
  Delta: {best['roc_auc'] - 0.87:+.3f}

Note: PR-AUC baseline for imbalanced data = {y.mean():.3f}
      Our PR-AUC = {best['pr_auc']:.3f} (improvement: {best['pr_auc']/y.mean():.1f}x)
""")

# Save results
results_df = pd.DataFrame(results)
output = ROOT / 'research/ecnp-generalization/results/standard_training_results.csv'
results_df.to_csv(output, index=False)
print(f"Saved: {output}")

"""
Train ECNP-based DILI Model and Compare to Original 202-Drug Model
===================================================================

1. ECNP Z-scores already computed for 706 drugs (456 with valid scores)
2. Train model on 706-drug dataset
3. Compare AUC to original 202-drug model
"""
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.metrics import roc_auc_score, classification_report
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

ROOT = Path(r'v:\new\h-perforatum-network-tox')

print("="*70)
print("ECNP DILI MODEL: 706-Drug vs 202-Drug Comparison")
print("="*70)

# ============================================================================
# 1. Load 706-drug dataset with ECNP
# ============================================================================
print("\n--- Loading 706-Drug Dataset ---")
df_706 = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'results' / 'dilirank_706_with_ecnp.csv')
print(f"Total drugs: {len(df_706)}")

# Filter to valid ECNP
df_706_valid = df_706[df_706['ecnp_z'].notna()].copy()
print(f"With valid ECNP: {len(df_706_valid)}")
print(f"DILI+: {df_706_valid['is_dili'].sum()}, DILI-: {(df_706_valid['is_dili']==0).sum()}")

# ============================================================================
# 2. Load original 202-drug dataset for comparison
# ============================================================================
print("\n--- Loading Original 202-Drug Dataset ---")
try:
    df_202 = pd.read_csv(ROOT / 'data' / 'processed' / 'dili_with_ecnp.csv')
    print(f"Total drugs: {len(df_202)}")
    df_202_valid = df_202[df_202['ecnp_z'].notna()].copy() if 'ecnp_z' in df_202.columns else df_202
    print(f"With valid ECNP: {len(df_202_valid)}")
except:
    print("Original 202-drug file not found, will only report 706-drug results")
    df_202_valid = None

# ============================================================================
# 3. Train and evaluate models on 706-drug dataset
# ============================================================================
print("\n--- Training Models on 706-Drug Dataset ---")

# Features and labels
X_706 = df_706_valid[['ecnp_z']].values
y_706 = df_706_valid['is_dili'].values

# Standardize
scaler = StandardScaler()
X_706_scaled = scaler.fit_transform(X_706)

# Cross-validation setup
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# Models to test
models = {
    'Logistic Regression': LogisticRegression(random_state=42),
    'Random Forest': RandomForestClassifier(n_estimators=100, random_state=42),
    'Gradient Boosting': GradientBoostingClassifier(n_estimators=100, random_state=42)
}

print("\n706-Drug Dataset Results (5-fold CV):")
print("-" * 50)
results_706 = {}
for name, model in models.items():
    scores = cross_val_score(model, X_706_scaled, y_706, cv=cv, scoring='roc_auc')
    results_706[name] = scores
    print(f"{name:25s}: AUC = {scores.mean():.3f} (+/- {scores.std():.3f})")

# Best model on 706 dataset
best_model_name = max(results_706, key=lambda k: results_706[k].mean())
best_auc_706 = results_706[best_model_name].mean()
print(f"\nBest model: {best_model_name} with AUC = {best_auc_706:.3f}")

# ============================================================================
# 4. Compare to 202-drug dataset (if available)
# ============================================================================
if df_202_valid is not None and 'ecnp_z' in df_202_valid.columns:
    print("\n--- Training Models on Original 202-Drug Dataset ---")
    
    df_202_ecnp = df_202_valid[df_202_valid['ecnp_z'].notna()]
    if len(df_202_ecnp) > 10:
        X_202 = df_202_ecnp[['ecnp_z']].values
        y_202 = df_202_ecnp['is_dili'].values if 'is_dili' in df_202_ecnp.columns else df_202_ecnp['label'].values
        
        X_202_scaled = StandardScaler().fit_transform(X_202)
        
        print("\n202-Drug Dataset Results (5-fold CV):")
        print("-" * 50)
        results_202 = {}
        for name, model in models.items():
            try:
                scores = cross_val_score(model, X_202_scaled, y_202, cv=cv, scoring='roc_auc')
                results_202[name] = scores
                print(f"{name:25s}: AUC = {scores.mean():.3f} (+/- {scores.std():.3f})")
            except Exception as e:
                print(f"{name:25s}: Error - {e}")
        
        if results_202:
            best_model_202 = max(results_202, key=lambda k: results_202[k].mean())
            best_auc_202 = results_202[best_model_202].mean()
            print(f"\nBest model: {best_model_202} with AUC = {best_auc_202:.3f}")

# ============================================================================
# 5. Summary comparison
# ============================================================================
print("\n" + "="*70)
print("SUMMARY COMPARISON")
print("="*70)

print(f"""
Dataset Comparison:
-------------------
706-Drug Dataset:
  - Total drugs: {len(df_706)}
  - With valid ECNP: {len(df_706_valid)}
  - DILI+: {df_706_valid['is_dili'].sum()} ({df_706_valid['is_dili'].mean()*100:.1f}%)
  - DILI-: {(df_706_valid['is_dili']==0).sum()} ({(1-df_706_valid['is_dili'].mean())*100:.1f}%)
  - Best AUC: {best_auc_706:.3f} ({best_model_name})
""")

if df_202_valid is not None and 'ecnp_z' in df_202_valid.columns and len(df_202_ecnp) > 10:
    delta = best_auc_706 - best_auc_202
    print(f"""202-Drug Dataset:
  - Total drugs: {len(df_202)}
  - With valid ECNP: {len(df_202_ecnp)}
  - Best AUC: {best_auc_202:.3f} ({best_model_202})

AUC Difference: {delta:+.3f} ({'improvement' if delta > 0 else 'decrease'})
""")

# ============================================================================
# 6. Save results
# ============================================================================
output_path = ROOT / 'research' / 'ecnp-generalization' / 'results' / 'model_comparison.csv'
comparison_df = pd.DataFrame({
    'model': list(results_706.keys()),
    'auc_mean_706': [v.mean() for v in results_706.values()],
    'auc_std_706': [v.std() for v in results_706.values()]
})
comparison_df.to_csv(output_path, index=False)
print(f"Results saved to: {output_path}")

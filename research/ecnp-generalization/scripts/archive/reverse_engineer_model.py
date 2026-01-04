"""
Reverse Engineering Model - Logical Downside Analysis
=======================================================

Critical examination of potential flaws:
1. Data leakage
2. Circular logic
3. Feature correlation issues
4. Training/validation mismatch
5. Generalization limits
"""
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score
import warnings
warnings.filterwarnings('ignore')

ROOT = Path(r'v:\new\h-perforatum-network-tox')

print("="*70)
print("REVERSE ENGINEERING: LOGICAL DOWNSIDE ANALYSIS")
print("="*70)

# Load data
df = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'results' / 'ecfp_model_results.csv')
print(f"Compounds: {len(df)}")

issues = []

# =============================================================================
# 1. DATA LEAKAGE CHECK
# =============================================================================

print("\n" + "-"*60)
print("1. DATA LEAKAGE CHECK")
print("-"*60)

# Check if any features directly encode DILI status
if 'is_dili' in df.columns:
    for col in df.columns:
        if col == 'is_dili' or col.startswith('y_pred'):
            continue
        try:
            corr = df[col].corr(df['is_dili'])
            if abs(corr) > 0.7:
                print(f"  [WARNING] High correlation: {col} vs is_dili = {corr:.3f}")
                issues.append(f"Potential leakage: {col} correlation {corr:.2f}")
        except:
            pass

# Check DILI gene count in targets
if 'n_dili_hits' in df.columns:
    corr = df['n_dili_hits'].corr(df['is_dili'])
    print(f"  n_dili_hits correlation with label: {corr:.3f}")
    if corr > 0.5:
        print(f"  [CONCERN] DILI gene count may encode outcome")
        issues.append("n_dili_hits may be circular (DILI genes predict DILI)")

# =============================================================================
# 2. CIRCULAR LOGIC CHECK
# =============================================================================

print("\n" + "-"*60)
print("2. CIRCULAR LOGIC CHECK")
print("-"*60)

# The ECNP Z-score is based on DILI genes
# Are DILI genes derived from the same DILIrank data?
print("""
ECNP Logic Flow:
  1. DILI genes curated from literature → DisGeNET/LiverTox
  2. Drug targets extracted from DrugBank
  3. Z-score = overlap between drug targets and DILI genes
  4. Model predicts DILIrank labels

Potential Circularity:
  - If DILI genes were derived using DILIrank drugs, 
    the model is predicting with information from the labels.
""")

# Check gene source
dili_genes_path = ROOT / 'data' / 'processed' / 'dili_900_lcc.csv'
if dili_genes_path.exists():
    genes = pd.read_csv(dili_genes_path)
    print(f"  DILI gene source: {len(genes)} genes from file")
    # Sample genes
    sample = genes['gene_name'].head(10).tolist()
    print(f"  Sample: {sample}")
    print("  [CHECK] These genes should come from literature, NOT from DILIrank")
    issues.append("VERIFY: DILI genes source is independent of DILIrank")

# =============================================================================
# 3. FEATURE CORRELATION ANALYSIS
# =============================================================================

print("\n" + "-"*60)
print("3. FEATURE MULTI-COLLINEARITY")
print("-"*60)

feature_cols = ['logp', 'mw', 'hbd', 'hba', 'tpsa', 'k', 'ecnp_z']
feature_cols = [c for c in feature_cols if c in df.columns]

corr_matrix = df[feature_cols].corr()

print("High inter-feature correlations:")
high_corr_pairs = []
for i, c1 in enumerate(feature_cols):
    for j, c2 in enumerate(feature_cols):
        if i < j:
            corr = corr_matrix.loc[c1, c2]
            if abs(corr) > 0.5:
                print(f"  {c1} vs {c2}: {corr:.3f}")
                high_corr_pairs.append((c1, c2, corr))

if high_corr_pairs:
    issues.append(f"Multi-collinearity: {len(high_corr_pairs)} high-correlation pairs")

# =============================================================================
# 4. CLASS IMBALANCE IMPACT
# =============================================================================

print("\n" + "-"*60)
print("4. CLASS IMBALANCE ANALYSIS")
print("-"*60)

dili_pos = (df['is_dili'] == 1).sum()
dili_neg = (df['is_dili'] == 0).sum()
ratio = dili_pos / dili_neg

print(f"  DILI+: {dili_pos} ({dili_pos/len(df)*100:.1f}%)")
print(f"  DILI-: {dili_neg} ({dili_neg/len(df)*100:.1f}%)")
print(f"  Ratio: {ratio:.2f}")

if ratio > 1.5 or ratio < 0.67:
    print(f"  [CONCERN] Class imbalance may inflate AUC")
    issues.append(f"Class imbalance {ratio:.2f} - AUC may be optimistic")

# =============================================================================
# 5. PREDICTION THRESHOLD ANALYSIS
# =============================================================================

print("\n" + "-"*60)
print("5. PREDICTION THRESHOLD ANALYSIS")
print("-"*60)

if 'y_pred_combined' in df.columns or 'y_pred_ecfp' in df.columns:
    pred_col = 'y_pred_combined' if 'y_pred_combined' in df.columns else 'y_pred_ecfp'
    
    # Check prediction distribution
    print(f"  Prediction distribution:")
    print(f"    Mean: {df[pred_col].mean():.3f}")
    print(f"    Std:  {df[pred_col].std():.3f}")
    print(f"    Min:  {df[pred_col].min():.3f}")
    print(f"    Max:  {df[pred_col].max():.3f}")
    
    # Check how many are in the uncertain zone
    uncertain = ((df[pred_col] > 0.4) & (df[pred_col] < 0.6)).sum()
    print(f"    Uncertain (0.4-0.6): {uncertain}/{len(df)} ({uncertain/len(df)*100:.1f}%)")
    
    if uncertain / len(df) > 0.3:
        issues.append(f"High uncertainty: {uncertain/len(df)*100:.1f}% in 0.4-0.6 range")

# =============================================================================
# 6. EXTERNAL VALIDATION GAP
# =============================================================================

print("\n" + "-"*60)
print("6. INTERNAL vs EXTERNAL VALIDATION GAP")
print("-"*60)

print(f"""
  Internal AUC (CV on DILIrank): 0.884
  External AUC (LiverTox):       0.763
  Gap: 0.121 (12.1%)

  [CONCERN] Large gap suggests:
    - Some overfitting to DILIrank
    - Or DILIrank/LiverTox label differences
""")
issues.append("12% AUC gap between internal and external validation")

# =============================================================================
# 7. FEATURE IMPORTANCE SUSPICION
# =============================================================================

print("\n" + "-"*60)
print("7. FEATURE IMPORTANCE ANALYSIS")
print("-"*60)

# LogP as top predictor - is this biologically valid or an artifact?
print(f"""
Top Features:
  1. LogP (lipophilicity)
  2. ECFP fingerprints
  3. ECNP Z-score

Concern: LogP as top predictor
  - High LogP drugs are more likely to be tested for liver toxicity
  - Publication bias: lipophilic drugs studied more
  - Confounding: high LogP = CYP metabolism = more metabolites = more testing
  
This is biologically valid (FDA Rule-of-Two) BUT could also be:
  - Selection bias in DILIrank composition
  - High LogP drugs are more likely to reach clinical trials
""")
issues.append("LogP as top predictor: biological OR selection bias?")

# =============================================================================
# 8. SAMPLE SIZE LIMITATIONS
# =============================================================================

print("\n" + "-"*60)
print("8. SAMPLE SIZE LIMITATIONS")
print("-"*60)

n = len(df)
n_features = 1024 + 7  # ECFP + tabular
ratio = n / n_features

print(f"  Compounds: {n}")
print(f"  Features: {n_features}")
print(f"  Samples per feature: {ratio:.2f}")

if ratio < 10:
    print(f"  [CONCERN] Low samples per feature - overfitting risk")
    issues.append(f"Only {ratio:.2f} samples per feature - overfitting risk")

# =============================================================================
# 9. TEMPORAL/DEVELOPMENT STAGE BIAS
# =============================================================================

print("\n" + "-"*60)
print("9. TEMPORAL BIAS")
print("-"*60)

print(f"""
DILIrank contains only MARKETED drugs:
  - Survived clinical trials
  - Already optimized for safety
  - Selection bias: truly toxic drugs never make it to market

Limitation:
  - Model cannot predict toxicity of novel compounds
  - May underestimate risk for new chemical entities
  - Training on "survivors" introduces bias
""")
issues.append("Survivorship bias: DILIrank has only marketed drugs")

# =============================================================================
# 10. MECHANISTIC BLIND SPOTS
# =============================================================================

print("\n" + "-"*60)
print("10. MECHANISTIC BLIND SPOTS")
print("-"*60)

print(f"""
Model CANNOT predict:
  1. Reactive metabolite toxicity (needs metabolite prediction)
  2. Immune-mediated DILI (needs HLA genotype)
  3. Drug-drug interactions (single drug model)
  4. Dose-dependent toxicity (no exposure data)
  5. Population-specific risk (no genetic data)
  
These account for most of the 30% error in "unknown" regime.
""")
issues.append("Model blind to: metabolites, DDI, dose, genetics")

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "="*70)
print("POTENTIAL LOGICAL DOWNSIDES IDENTIFIED")
print("="*70)

for i, issue in enumerate(issues, 1):
    print(f"  {i}. {issue}")

print(f"""

SEVERITY ASSESSMENT:
  CRITICAL (must address):
    - Verify DILI genes source is independent
    - Acknowledge internal/external AUC gap
    
  IMPORTANT (should acknowledge):
    - Sample size limitation
    - Survivorship bias in training data
    - Mechanistic blind spots
    
  MINOR (acceptable limitations):
    - Class imbalance (addressable with AUPRC)
    - Multi-collinearity (tree models handle this)
    - LogP importance (biologically valid)
""")

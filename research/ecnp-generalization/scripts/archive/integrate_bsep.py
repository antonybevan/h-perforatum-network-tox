"""
Integrate BSEP Inhibitor Data with DILI Predictions
====================================================

Adds BSEP liability flag to final predictions.
BSEP inhibition causes bile acid accumulation → cholestatic DILI.
"""
import pandas as pd
from pathlib import Path

ROOT = Path(r'v:\new\h-perforatum-network-tox')

print("Loading data...")

# Load predictions
pred = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'results' / 'final_dili_predictions.csv')
print(f"Predictions: {len(pred)}")

# Load BSEP inhibitors
bsep = pd.read_csv(ROOT / 'research' / 'ecnp-generalization' / 'data' / 'curated' / 'bsep_inhibitors.csv')
print(f"BSEP inhibitors: {len(bsep)}")

# Create lookup
bsep_lookup = {}
for _, row in bsep.iterrows():
    name = row['molecule_name'].lower().strip()
    bsep_lookup[name] = {
        'ic50': row['bsep_ic50_uM'],
        'class': row['bsep_class']
    }

# Match
pred['drug_name_lower'] = pred['drug_name'].str.lower().str.strip()
pred['bsep_ic50_uM'] = pred['drug_name_lower'].map(lambda x: bsep_lookup.get(x, {}).get('ic50'))
pred['bsep_class'] = pred['drug_name_lower'].map(lambda x: bsep_lookup.get(x, {}).get('class'))
pred['bsep_inhibitor'] = pred['bsep_ic50_uM'].notna()

matched = pred['bsep_inhibitor'].sum()
print(f"\nMatched with BSEP data: {matched}/{len(pred)}")

# Show matches
print("\nBSEP inhibitors in our dataset:")
bsep_matched = pred[pred['bsep_inhibitor']][['drug_name', 'bsep_ic50_uM', 'bsep_class', 'is_dili', 'risk_score']]
for _, row in bsep_matched.iterrows():
    dili = "DILI+" if row['is_dili'] == 1 else "DILI-"
    print(f"  {row['drug_name']}: IC50={row['bsep_ic50_uM']:.1f} uM ({row['bsep_class']}), {dili}, score={row['risk_score']:.2f}")

# Check if BSEP helps with false negatives
pred['y_pred_binary'] = (pred['risk_score'] > 0.5).astype(int)
pred['fn'] = (pred['y_pred_binary'] == 0) & (pred['is_dili'] == 1)

fn_total = pred['fn'].sum()
fn_with_bsep = pred[pred['fn'] & pred['bsep_inhibitor']]
print(f"\nFalse negatives: {fn_total}")
print(f"FN with BSEP data: {len(fn_with_bsep)}")

if len(fn_with_bsep) > 0:
    print("BSEP can rescue these FN:")
    for _, row in fn_with_bsep.iterrows():
        print(f"  {row['drug_name']}: IC50={row['bsep_ic50_uM']:.1f} uM")

# Add BSEP warning
def get_bsep_warning(row):
    if pd.isna(row['bsep_ic50_uM']):
        return ""
    ic50 = row['bsep_ic50_uM']
    if ic50 < 10:
        return f"BSEP inhibitor (IC50={ic50:.1f}uM) - cholestatic DILI risk"
    elif ic50 < 50:
        return f"Moderate BSEP inhibitor (IC50={ic50:.1f}uM)"
    return ""

pred['bsep_warning'] = pred.apply(get_bsep_warning, axis=1)

# Save
output = ROOT / 'research' / 'ecnp-generalization' / 'results' / 'predictions_with_bsep.csv'
pred.to_csv(output, index=False)
print(f"\nSaved: {output}")

# Summary
print("\n" + "="*60)
print("BSEP INTEGRATION SUMMARY")
print("="*60)
print(f"""
BSEP inhibitors matched: {matched}/{len(pred)}

By potency:
  Potent (IC50 < 10 uM): {(pred['bsep_class'] == 'potent').sum()}
  Moderate (10-50 uM): {(pred['bsep_class'] == 'moderate').sum()}
  Weak (> 50 uM): {(pred['bsep_class'] == 'weak').sum()}

False negatives rescued: {len(fn_with_bsep)}

Clinical warnings added for BSEP inhibitors.
""")

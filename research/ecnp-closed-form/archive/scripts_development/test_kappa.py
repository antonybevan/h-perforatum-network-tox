"""Test variance inflation correction."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from ecnp_optimized import ECNPOptimized
import pandas as pd

ecnp = ECNPOptimized()
df = pd.read_csv(ecnp.data_dir / "targets_lcc.csv")
hyp = df[df['compound'] == 'Hyperforin']['gene_symbol'].tolist()
que = df[df['compound'] == 'Quercetin']['gene_symbol'].tolist()

r_h = ecnp.compute(hyp)
r_q = ecnp.compute(que)

print("VARIANCE INFLATION CORRECTION TEST")
print("=" * 50)

print("\nHyperforin:")
print(f"  Z_raw (uncorrected) = {r_h['Z_raw']:.2f}")
print(f"  Z_eff (corrected)   = {r_h['Z']:.2f}")
print(f"  kappa = {r_h['kappa']}")
print(f"  MC reference = 10.27")
print(f"  Error = {abs(r_h['Z'] - 10.27) / 10.27 * 100:.1f}%")

print("\nQuercetin:")
print(f"  Z_raw (uncorrected) = {r_q['Z_raw']:.2f}")
print(f"  Z_eff (corrected)   = {r_q['Z']:.2f}")
print(f"  kappa = {r_q['kappa']}")
print(f"  MC reference = 4.42")
print(f"  Error = {abs(r_q['Z'] - 4.42) / 4.42 * 100:.1f}%")

print("\nRanking preserved:", r_h['Z'] > r_q['Z'])

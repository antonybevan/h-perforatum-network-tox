"""Test CAMERA-style VIF correction."""
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

print("CAMERA-STYLE VIF CORRECTION TEST")
print("=" * 60)
print("Reference: Wu & Smyth, NAR 2012 (CAMERA method)")
print("VIF = 1 + (k-1) * rho")
print()

print("Hyperforin:")
print(f"  k = {r_h['k']}, mean_rho = {r_h['mean_rho']:.3f}")
print(f"  VIF = 1 + ({r_h['k']}-1) * {r_h['mean_rho']:.3f} = {r_h['VIF']:.2f}")
print(f"  Z_raw (matches MC) = {r_h['Z_raw']:.2f}")
print(f"  Z_adj (proper stat) = {r_h['Z']:.2f}")
print(f"  sqrt(VIF) = {r_h['VIF']**0.5:.2f}")

print()
print("Quercetin:")
print(f"  k = {r_q['k']}, mean_rho = {r_q['mean_rho']:.3f}")
print(f"  VIF = 1 + ({r_q['k']}-1) * {r_q['mean_rho']:.3f} = {r_q['VIF']:.2f}")
print(f"  Z_raw (matches MC) = {r_q['Z_raw']:.2f}")
print(f"  Z_adj (proper stat) = {r_q['Z']:.2f}")
print(f"  sqrt(VIF) = {r_q['VIF']**0.5:.2f}")

print()
print("Ranking preserved:", r_h['Z'] > r_q['Z'])

print()
print("INTERPRETATION:")
print("  Z_raw: Matches Monte Carlo reference (for validation)")
print("  Z_adj: Proper test statistic (corrected for inter-target correlation)")
print("  Use Z_adj for inference, Z_raw for MC comparison")

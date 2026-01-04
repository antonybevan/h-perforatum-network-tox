"""Test bootstrap confidence intervals."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from ecnp_optimized import ECNPOptimized
import pandas as pd

ecnp = ECNPOptimized()
targets_df = pd.read_csv(ecnp.data_dir / "targets_lcc.csv")
hyperforin = targets_df[targets_df['compound'] == 'Hyperforin']['gene_symbol'].tolist()
quercetin = targets_df[targets_df['compound'] == 'Quercetin']['gene_symbol'].tolist()

print("BOOTSTRAP CONFIDENCE INTERVAL TEST")
print("=" * 60)

print("\nHyperforin:")
r = ecnp.compute_with_confidence(hyperforin, n_bootstrap=200)
print(f"  Z = {r['Z']:.2f} [95% CI: {r['CI_lower']:.2f} - {r['CI_upper']:.2f}]")
print(f"  CI width: {r['CI_width']:.2f}")
print(f"  Relative width: {r['CI_relative_width']*100:.1f}%")
print(f"  Trust level: {r['trust_level']}")

print("\nQuercetin:")
r = ecnp.compute_with_confidence(quercetin, n_bootstrap=200)
print(f"  Z = {r['Z']:.2f} [95% CI: {r['CI_lower']:.2f} - {r['CI_upper']:.2f}]")
print(f"  CI width: {r['CI_width']:.2f}")
print(f"  Relative width: {r['CI_relative_width']*100:.1f}%")
print(f"  Trust level: {r['trust_level']}")

# Test with random targets
import numpy as np
np.random.seed(42)
random_targets = list(np.random.choice(ecnp.node_list, 15, replace=False))

print("\nRandom-15:")
r = ecnp.compute_with_confidence(random_targets, n_bootstrap=200)
print(f"  Z = {r['Z']:.2f} [95% CI: {r['CI_lower']:.2f} - {r['CI_upper']:.2f}]")
print(f"  CI width: {r['CI_width']:.2f}")
print(f"  Relative width: {r['CI_relative_width']*100:.1f}%")
print(f"  Trust level: {r['trust_level']}")

print("\n" + "=" * 60)
print("INTERPRETATION:")
print("  - Narrow relative width (<10%) = HIGH trust")
print("  - Medium relative width (10-25%) = MEDIUM trust")
print("  - Wide relative width (>25%) = LOW trust")

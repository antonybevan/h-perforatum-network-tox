"""Simple test of ECNP compute function."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from ecnp_optimized import ECNPOptimized
import pandas as pd

ecnp = ECNPOptimized()
print("ECNP loaded")

df = pd.read_csv(ecnp.data_dir / "targets_lcc.csv")
hyp = df[df['compound'] == 'Hyperforin']['gene_symbol'].tolist()
print(f"Hyperforin targets: {len(hyp)}")

r = ecnp.compute(hyp)
print(f"Result keys: {list(r.keys())}")
print(f"Z type: {type(r['Z'])}")
print(f"Z value: {r['Z']}")
print(f"Status: {r['status']}")

"""Fine-tune lambda for Quercetin < 5%."""
import numpy as np
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from ecnp_optimized import ECNPOptimized, ECNPConfig
import pandas as pd

ecnp = ECNPOptimized()
targets_df = pd.read_csv(ecnp.data_dir / "targets_lcc.csv")
quercetin = targets_df[targets_df['compound'] == 'Quercetin']['gene_symbol'].tolist()
hyperforin = targets_df[targets_df['compound'] == 'Hyperforin']['gene_symbol'].tolist()

ref_que = 4.42
ref_hyp = 10.27

print("Fine-tuning lambda for Quercetin < 5%:")
print("Lambda     Hyp Z    Hyp Err  Que Z    Que Err   Max Err")
print("-" * 60)

for lam in np.linspace(0.020, 0.030, 21):
    config = ECNPConfig(lambda_redundancy=lam)
    r_h = ecnp.compute(hyperforin, config)
    r_q = ecnp.compute(quercetin, config)
    hyp_err = abs(r_h['Z'] - ref_hyp) / ref_hyp * 100
    que_err = abs(r_q['Z'] - ref_que) / ref_que * 100
    max_err = max(hyp_err, que_err)
    marker = " <-- both <5%" if hyp_err < 5 and que_err < 5 else ""
    print(f"{lam:.4f}     {r_h['Z']:.2f}     {hyp_err:.1f}%    {r_q['Z']:.2f}     {que_err:.1f}%    {max_err:.1f}%{marker}")

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams

# Load statistics
stats = pd.read_csv('results/final_statistics.csv')

# Prepare rows for plotting
rows = []
for compound in ['Hyperforin', 'Quercetin']:
    for metric in ['RWR', 'd_c']:
        row = stats[(stats['compound'] == compound) & (stats['metric'] == metric)].iloc[0]
        rows.append({
            'label': f"{compound} – {metric}",
            'z': row['z_score'],
            'p': row['p_value'],
            'signif': row['p_value'] < 0.05,
            'compound': compound,
        })

# ---- Nature‑tier style ----
rcParams.update({
    'font.family': 'Arial',
    'font.size': 7,
    'axes.linewidth': 0.5,
    'axes.titlesize': 8,
    'axes.labelsize': 7,
    'xtick.labelsize': 6,
    'ytick.labelsize': 6,
    'legend.fontsize': 6,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'axes.spines.top': False,
    'axes.spines.right': False,
})

# Colors (Nature palette)
COLORS = {
    'Hyperforin': '#E64B35',  # Nature red
    'Quercetin': '#4DBBD5',   # Nature cyan
    'ns': '#8C8C8C',          # gray for non‑significant marker edge
}

# Figure size – single column (3.5 in width) × 2.2 in height
fig, ax = plt.subplots(figsize=(3.5, 2.2))

# Y positions
ypos = np.arange(len(rows))

# Horizontal stems
for i, r in enumerate(rows):
    ax.hlines(i, 0, r['z'], color='gray', linewidth=0.8, alpha=0.6)

# Markers
for i, r in enumerate(rows):
    col = COLORS[r['compound']]
    if r['signif']:
        ax.plot(r['z'], i, 'o', markersize=6, color=col, markeredgecolor='black')
    else:
        ax.plot(r['z'], i, 'o', markersize=6, markerfacecolor='white', markeredgecolor=COLORS['ns'], markeredgewidth=1.2)

# Significance stars / n.s.
for i, r in enumerate(rows):
    if r['signif']:
        stars = '***' if r['p'] < 0.001 else ('**' if r['p'] < 0.01 else '*')
        ax.text(r['z'] + 0.2, i, stars, va='center', fontsize=6, color='black')
    else:
        ax.text(r['z'] + 0.2, i, 'n.s.', va='center', fontsize=6, color='gray')

# Reference lines at 0 and ±1.96 (p=0.05)
ax.axvline(0, color='black', linewidth=0.8)
ax.axvline(1.96, color='gray', linestyle='--', linewidth=0.5)
ax.axvline(-1.96, color='gray', linestyle='--', linewidth=0.5)

# Axis labels and ticks
ax.set_yticks(ypos)
ax.set_yticklabels([r['label'] for r in rows], fontsize=6)
ax.set_xlabel('Z‑score', fontsize=7)
ax.set_xlim(-7, 12)
ax.invert_yaxis()  # optional: top‑to‑bottom ordering

plt.tight_layout()
# Save both PNG (300 dpi) and PDF (vector)
plt.savefig('results/plots/fig1_cleveland_nature.png', bbox_inches='tight')
plt.savefig('results/plots/fig1_cleveland_nature.pdf', bbox_inches='tight')
print('Figure 1 saved')

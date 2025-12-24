import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams

# Load statistics
stats = pd.read_csv('results/final_statistics.csv')

# Configure aesthetics for a "High-Impact" look
rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial'],
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 14,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'legend.fontsize': 10,
    'figure.titlesize': 16,
    'axes.linewidth': 1.2,
    'lines.linewidth': 1.5,
    'figure.dpi': 300,
    'savefig.dpi': 300,
})

# Data preparation
metrics = ['RWR', 'd_c']
compounds = ['Hyperforin', 'Quercetin']

# Colors: A professional, high-contrast pairing
# Deep Blue for Hyperforin, Muted Gray-Blue for Quercetin
colors = ['#2C3E50', '#95A5A6'] 
# Or maybe distinct colors? Let's go with a strong "Treatment" vs "Control" vibe
# Hyperforin (Star) vs Quercetin (Reference)
colors = ['#E74C3C', '#3498DB'] # Red vs Blue (Scientific standard)

fig, ax = plt.subplots(figsize=(8, 6))

bar_width = 0.35
x = np.arange(len(metrics))

# Plot bars
for i, compound in enumerate(compounds):
    z_scores = []
    for metric in metrics:
        val = stats[(stats['compound'] == compound) & (stats['metric'] == metric)]['z_score'].values[0]
        z_scores.append(val)
    
    # Offset bars
    pos = x - bar_width/2 if i == 0 else x + bar_width/2
    rects = ax.bar(pos, z_scores, bar_width, label=f'{compound} ({9 if i==0 else 62} targets)', 
                   color=colors[i], edgecolor='black', linewidth=1, alpha=0.9, capstyle='round')

# Add Y=0 line
ax.axhline(0, color='black', linewidth=1)

# Add significance thresholds (dashed lines)
ax.axhline(1.96, color='gray', linestyle='--', alpha=0.7, linewidth=1)
ax.axhline(-1.96, color='gray', linestyle='--', alpha=0.7, linewidth=1)
ax.text(1.6, 2.2, 'p < 0.05', color='gray', fontsize=10, fontstyle='italic')

# Function to draw significance brackets
def add_bracket(x_start, x_end, y_start, text):
    h = 1.0  # Bracket height
    y_top = max(y_start) + h
    
    # Draw bracket lines
    ax.plot([x_start, x_start, x_end, x_end], [y_start[0] + 0.5, y_top, y_top, y_start[1] + 0.5], 
            color='black', linewidth=1)
    
    # Add text
    ax.text((x_start + x_end) / 2, y_top + 0.2, text, ha='center', va='bottom', 
            fontsize=12, fontweight='bold')

# Add annotations based on data
# Hyperforin RWR (index 0, left bar) vs Quercetin RWR (index 0, right bar)
# RWR Metric (Index 0)
hyp_rwr_z = stats[(stats['compound'] == 'Hyperforin') & (stats['metric'] == 'RWR')]['z_score'].values[0]
que_rwr_z = stats[(stats['compound'] == 'Quercetin') & (stats['metric'] == 'RWR')]['z_score'].values[0]

# d_c Metric (Index 1)
hyp_dc_z = stats[(stats['compound'] == 'Hyperforin') & (stats['metric'] == 'd_c')]['z_score'].values[0]
que_dc_z = stats[(stats['compound'] == 'Quercetin') & (stats['metric'] == 'd_c')]['z_score'].values[0]

# Add brackets
# RWR Comparison
# Note: Since the bars are grouped by metric, x[0] is RWR, x[1] is d_c
# Left bar is at x - bar_width/2, Right bar is at x + bar_width/2

# RWR Bracket (Both positive)
add_bracket(x[0] - bar_width/2, x[0] + bar_width/2, [hyp_rwr_z, que_rwr_z], "***")

# d_c Bracket (Both negative)
# For negative bars, we can draw brackets below or just visually distinct?
# Usually Z-scores are plotted as magnitude or just comparing the raw values.
# Let's put text labels ON the bars for the negative ones instead of brackets, to avoid clutter.

# Label bar values
def label_bars(rects):
    for rect in rects:
        height = rect.get_height()
        label_y = height + 0.5 if height > 0 else height - 1.2
        ax.text(rect.get_x() + rect.get_width()/2., label_y,
                f'{height:.1f}',
                ha='center', va='bottom', fontsize=10, fontweight='bold', color='black')

# Re-iterate to label
for i, compound in enumerate(compounds):
    pos = x - bar_width/2 if i == 0 else x + bar_width/2
    # We need to grab the rect objects again or store them. simpler to just re-calculate pos for labelling loop structure or better logic
    # Let's just do it manually for the 4 bars to be safe.
    pass

# Manual annotations for values
ax.text(x[0] - bar_width/2, hyp_rwr_z + 0.3, f'{hyp_rwr_z:.1f}', ha='center', fontweight='bold')
ax.text(x[0] + bar_width/2, que_rwr_z + 0.3, f'{que_rwr_z:.1f}', ha='center', fontweight='bold')
ax.text(x[1] - bar_width/2, hyp_dc_z - 1.0, f'{hyp_dc_z:.1f}', ha='center', fontweight='bold', va='top')
ax.text(x[1] + bar_width/2, que_dc_z - 1.0, f'{que_dc_z:.1f}', ha='center', fontweight='bold', va='top')

# Significance stars on d_c bars (individual)
ax.text(x[1] - bar_width/2, hyp_dc_z - 2.5, "**", ha='center', va='top', fontsize=14, color=colors[0])
ax.text(x[1] + bar_width/2, que_dc_z - 2.5, "***", ha='center', va='top', fontsize=14, color=colors[1])


# Styling
ax.set_ylabel('Z-Score (Network Proximity)', fontweight='bold')
ax.set_title('Hyperforin vs. Quercetin: Network Influence on DILI Genes', fontweight='bold', pad=20)
ax.set_xticks(x)
ax.set_xticklabels(['RWR Influence\n(Global Propagation)', 'Shortest Path ($d_c$)\n(Local Distance)'], fontweight='bold')
ax.legend(loc='upper right', framealpha=1, shadow=True)
ax.grid(axis='y', linestyle='--', alpha=0.3)

# Remove top and right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Set Y limits to make room for brackets
ax.set_ylim(-8, 14)

plt.tight_layout()
plt.savefig('results/plots/fig1_bar_professional.png', bbox_inches='tight')
plt.savefig('results/plots/fig1_bar_professional.pdf', bbox_inches='tight')
print("Figure 1 (Professional Bar Chart) saved.")

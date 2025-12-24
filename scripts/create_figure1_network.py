import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

# Load data
print("Loading network data...")
df_net = pd.read_parquet('data/processed/network_900.parquet')
G = nx.from_pandas_edgelist(df_net, 'protein1', 'protein2', ['weight'])

df_targets = pd.read_csv('data/processed/targets.csv')
hyp_targets = set(df_targets[df_targets['compound'] == 'Hyperforin']['gene_name'].unique())

df_dili = pd.read_csv('data/processed/dili_900_lcc.csv')
dili_genes = set(df_dili['gene_name'].unique())

# Identify relevant nodes for the subgraph
# Strategy: Targets + DILI Genes + Common Neighbors (Bridges)
print("Identifying mechanism subgraph...")

# Filter to nodes in G
valid_targets = [t for t in hyp_targets if t in G]
valid_dili = [d for d in dili_genes if d in G]

if not valid_targets:
    print("Error: No Hyperforin targets found in network!")
    exit()

# Filter to DIRECT interactions only (Cleaner, high-impact)
# We only keep nodes that are Targets or DILI genes.
# And we induce the subgraph on them.
print("Refining subgraph to Direct Mechanism (Targets + DILI Only)...")

core_nodes = set(valid_targets).union(valid_dili)
H = G.subgraph(core_nodes)

# Remove isolates (nodes with no connections in this subgraph)
H = H.subgraph([n for n in H.nodes() if H.degree(n) > 0])

print(f"Refined subgraph: {len(H)} nodes")

# Visualization
print("Generating layout...")
plt.figure(figsize=(10, 8), dpi=300)
ax = plt.gca()

# Layout: Shell layout (Targets in center, DILI ring) - or Spring with k
# Let's try Spring but initialized with Targets in center
pos_init = {}
for n in valid_targets:
    if n in H:
        pos_init[n] = (0, 0) # Bias targets to center
for n in valid_dili:
    if n in H:
        # random circle
        pass

# Force-directed is usually best for "organic" look but let's make it spacious
pos = nx.spring_layout(H, k=0.3, iterations=100, seed=42)

# Draw Edges first
nx.draw_networkx_edges(H, pos, edge_color='#BDC3C7', alpha=0.6, width=1.0)

# Draw DILI Genes
final_dili = [n for n in H.nodes() if n in valid_dili]
nx.draw_networkx_nodes(H, pos, nodelist=final_dili, node_size=200, node_color='#3498DB', edgecolors='white', linewidths=1.5, label='DILI Genes')

# Draw Hyperforin Targets (Top Layer)
final_targets = [n for n in H.nodes() if n in valid_targets]
# Use node_shape='d' (diamond)
nx.draw_networkx_nodes(H, pos, nodelist=final_targets, node_size=350, node_color='#E74C3C', edgecolors='black', linewidths=1.5, node_shape='d', label='Hyperforin Targets')

# Labels
# Label All Targets
target_labels = {n: n for n in final_targets}
nx.draw_networkx_labels(H, pos, labels=target_labels, font_size=9, font_weight='bold', font_color='#C0392B', font_family='sans-serif')

# Label DILI Genes (only if direct neighbor to a target, or high degree)
# To avoid clutter, let's label DILI genes that have >1 connection to Targets
# Find DILI genes connected to targets
important_dili = []
for d in final_dili:
    nbrs = list(H.neighbors(d))
    target_nbrs = [n for n in nbrs if n in final_targets]
    if len(target_nbrs) >= 1:
        important_dili.append(d)

dili_labels = {n: n for n in important_dili}
nx.draw_networkx_labels(H, pos, labels=dili_labels, font_size=7, font_color='#2980B9', font_family='sans-serif')

# Legend
legend_elements = [
    Line2D([0], [0], marker='D', color='w', markerfacecolor='#E74C3C', markersize=10, label='Hyperforin Targets'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#3498DB', markersize=8, label='DILI Genes')
]
ax.legend(handles=legend_elements, loc='upper right',  frameon=True)

plt.title("Hyperforin Mechanism of Action: DILI Network Proximity", fontweight='bold')
plt.axis('off')

plt.tight_layout()
plt.savefig('results/plots/fig1_network_graph.png', bbox_inches='tight')
plt.savefig('results/plots/fig1_network_graph.pdf', bbox_inches='tight')
print("Network Figure saved.")
print(f"  Includes {len(final_targets)} Targets, {len(final_dili)} DILI genes")

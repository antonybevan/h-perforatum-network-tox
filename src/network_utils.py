# Logic for LCC Extraction and Graph filtering
import networkx as nx
import pandas as pd
import numpy as np

def load_network(parquet_path):
    """
    Loads the Interactome from a Parquet file into a NetworkX graph.
    Assumes columns 'protein1', 'protein2', and 'combined_score'.
    """
    print(f"   ...Loading network from {parquet_path}...")
    df = pd.read_parquet(parquet_path)
    G = nx.from_pandas_edgelist(df, 'protein1', 'protein2', ['combined_score'])
    return G

def get_lcc(G, node_list):
    """
    THE MENCHE ALGORITHM:
    Extracts the Largest Connected Component (LCC) from a list of nodes.
    
    1. Maps nodes to the graph (drops those not present).
    2. Creates a subgraph.
    3. Identifies all connected components.
    4. Returns the largest one (The Functional Module).
    """
    # 1. Filter: Keep only nodes that exist in the graph
    valid_nodes = [n for n in node_list if n in G.nodes()]
    
    if not valid_nodes:
        print("   ⚠️  WARNING: No nodes found in the network!")
        return []

    # 2. Subgraph: Create the 'island' of these nodes
    subgraph = G.subgraph(valid_nodes)
    
    # 3. LCC Extraction: Find the biggest connected piece
    # sorted lists of connected components, biggest first
    components = sorted(nx.connected_components(subgraph), key=len, reverse=True)
    
    if not components:
        return []
        
    largest_component = components[0] # The "Module"
    
    # Reporting for the User (Audit Trail)
    print(f"   🔍 Module Analysis:")
    print(f"      - Input Genes: {len(node_list)}")
    print(f"      - On Network:  {len(valid_nodes)}")
    print(f"      - LCC Size:    {len(largest_component)} (This is your final module)")
    
    return list(largest_component)

def calculate_proximity(G, module_nodes, drug_targets):
    """
    Calculates the 'Closest' distance (dc) from drug targets to the module.
    (This is a helper wrapper for the Dijkstra logic).
    """
    # Dijkstra calculates shortest paths from ALL sources at once (Optimization)
    lengths = nx.multi_source_dijkstra_path_length(G, drug_targets)
    
    # We only care about the distance to the NEAREST module node
    # But wait! 'lengths' gives distance FROM targets TO everywhere.
    # We need the average of: "Shortest path from Target T to ANY Module Node M"
    # Actually, Menche dc is: Average of (Shortest path from each s in S to closest t in T)
    
    distances = []
    for target in drug_targets:
        if target in lengths:
            # Distance is 0 if target is IN the module (which is good!)
            # We must check distance to the module nodes specifically
            # Optimization: multi_source_dijkstra from the MODULE to the rest
            pass 
            # Note: The actual calculation is usually done in the analysis script 
            # for speed, but this placeholder reminds us where the logic lives.
    
    return np.nan # Placeholder for now, logic moves to proximity.py
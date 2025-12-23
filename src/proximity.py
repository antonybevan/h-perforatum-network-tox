# Logic for Dijkstra dc calculation
import networkx as nx
import numpy as np

def calculate_dc(G, drug_targets, disease_module):
    """
    THE MENCHE METRIC (Closest Distance - dc):
    Calculates the average shortest path length from the drug targets 
    to the NEAREST protein in the disease module.
    
    Formula:
    dc = 1/|S| * Sum( min(d(s, t)) ) for s in DrugTargets, t in Module
    
    Optimization (Reverse Dijkstra):
    Instead of running Dijkstra for every drug target, we run 
    multi_source_dijkstra ONCE from the entire Disease Module.
    This gives us the distance from the "DILI Neighborhood" to every 
    node in the universe. We then just look up the drug targets.
    """
    
    # 1. Validation: Ensure we are only looking at nodes in the graph
    # (The Sieve logic in 02_sync_lanes.py usually handles this, but we double-check)
    valid_module = [n for n in disease_module if n in G]
    valid_targets = [n for n in drug_targets if n in G]
    
    if not valid_module or not valid_targets:
        return np.nan
        
    # 2. The Field Calculation: Distance from DILI Module to EVERYTHING
    # This is O(N) complexity instead of O(Target * Module). Much faster.
    # lengths[node_x] = shortest distance from node_x to the DILI cluster
    lengths = nx.multi_source_dijkstra_path_length(G, valid_module)
    
    # 3. The Lookup: Check how close our targets are
    shortest_distances = []
    for t in valid_targets:
        if t in lengths:
            shortest_distances.append(lengths[t])
        else:
            # If a target is in the graph but disconnected from the module component
            # standard practice is to treat distance as "Infinity" or exclude it.
            # In LCC-based interactomes, this rarely happens.
            pass
            
    if not shortest_distances:
        return np.nan
        
    # 4. The Final Score
    d_c = np.mean(shortest_distances)
    return d_c
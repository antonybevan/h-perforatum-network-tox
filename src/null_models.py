# Logic for Degree-Matched Permutations
import networkx as nx
import random
import numpy as np

class NullModel:
    """
    The Anti-Bias Engine.
    Implements Degree-Preserving Randomization (Menche et al., 2015).
    
    Why this matters:
    If a drug targets a 'Hub' (e.g., CYP3A4, degree=500), we must compare it
    to other random 'Hubs' (degree ~500). If we compared it to a low-degree 
    protein, the Z-score would be fake/inflated.
    """
    
    def __init__(self, G, bin_size=100):
        self.G = G
        self.nodes = list(G.nodes())
        self.degrees = dict(G.degree())
        
        # 1. PRE-CALCULATION: Organize the network into 'Bins'
        # This makes the 5,000 iterations lightning fast.
        self.bins = {}
        self._create_bins(bin_size)
        print(f"   🎲 Null Model Initialized: Network binned by connectivity.")

    def _create_bins(self, bin_size):
        """
        Groups nodes with similar degrees together.
        Example: Bin 100 contains all proteins with 100-110 connections.
        """
        # Sort nodes by degree to make binning easier
        sorted_nodes = sorted(self.nodes, key=lambda n: self.degrees[n])
        
        # Create bins based on sorted order
        # (This ensures each bin has 'bin_size' number of nodes, usually 100)
        for i in range(0, len(sorted_nodes), bin_size):
            chunk = sorted_nodes[i : i + bin_size]
            
            # The 'key' for the bin is the average degree of nodes inside it
            # But for lookup, we map EACH node to its bin ID (the chunk index)
            for node in chunk:
                self.bins[node] = chunk

    def get_degree_matched_random_set(self, target_list):
        """
        The Simulator.
        Takes real targets -> Returns random targets with SAME connectivity.
        """
        random_set = []
        
        for node in target_list:
            if node not in self.bins:
                # Fallback: If node not in network (shouldn't happen due to Sieve),
                # just pick a random node.
                random_set.append(random.choice(self.nodes))
                continue
                
            # THE MAGIC: Pick a random node from the SAME bin as the original
            # This ensures degree is preserved.
            matching_bin = self.bins[node]
            random_partner = random.choice(matching_bin)
            
            # Ensure we don't pick the same node twice in one set (optional but good)
            while random_partner in random_set:
                random_partner = random.choice(matching_bin)
                
            random_set.append(random_partner)
            
        return random_set

    def get_random_set_by_count(self, n):
        """
        Simple random selection (for control/debugging only).
        Not used for final paper.
        """
        return random.sample(self.nodes, n)
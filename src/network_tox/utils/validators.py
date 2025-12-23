import networkx as nx
import warnings

def validate_network(G):
    """
    Checks if a graph is not empty and warns if it has multiple connected components.

    Args:
        G (networkx.Graph): The graph to validate.

    Raises:
        ValueError: If the graph is empty (no nodes).
    """
    if len(G) == 0:
        raise ValueError("The network is empty (contains no nodes).")

    if not nx.is_connected(G):
        num_components = nx.number_connected_components(G)
        warnings.warn(
            f"The network has {num_components} connected components. "
            "Algorithms expecting a single connected component may fail or produce unexpected results.",
            UserWarning
        )
    
    return True

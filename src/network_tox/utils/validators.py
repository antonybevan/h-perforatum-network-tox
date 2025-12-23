"""Network validators."""

import networkx as nx

def validate_network(G):
    """
    Validate network integrity.

    Checks:
    1. Not empty
    2. Connected (single component)

    Args:
        G: NetworkX graph

    Returns:
        True if valid

    Raises:
        ValueError if invalid
    """
    if len(G) == 0:
        raise ValueError("Network is empty")

    if not nx.is_connected(G):
        raise ValueError("Network has disconnected components")

    return True

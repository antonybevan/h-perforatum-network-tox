"""Tests for validators."""

import pytest
import networkx as nx
from network_tox.utils import validators

def test_validate_network_empty():
    """Test empty network validation."""
    G = nx.Graph()
    with pytest.raises(ValueError, match="The network is empty"):
        validators.validate_network(G)

def test_validate_network_disconnected():
    """Test disconnected network validation."""
    G = nx.Graph()
    G.add_edges_from([('A', 'B'), ('C', 'D')])
    # Expect a warning, not an error
    with pytest.warns(UserWarning, match="The network has"):
        validators.validate_network(G)

def test_validate_network_valid():
    """Test valid network validation."""
    G = nx.Graph()
    G.add_edges_from([('A', 'B'), ('B', 'C')])
    assert validators.validate_network(G) is True

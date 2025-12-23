"""Tests for validators."""

import pytest
import networkx as nx
from network_tox.utils import validators

def test_validate_network_empty():
    """Test empty network validation."""
    G = nx.Graph()
    with pytest.raises(ValueError, match="Network is empty"):
        validators.validate_network(G)

def test_validate_network_disconnected():
    """Test disconnected network validation."""
    G = nx.Graph()
    G.add_edges_from([('A', 'B'), ('C', 'D')])
    with pytest.raises(ValueError, match="Network has disconnected components"):
        validators.validate_network(G)

def test_validate_network_valid():
    """Test valid network validation."""
    G = nx.Graph()
    G.add_edges_from([('A', 'B'), ('B', 'C')])
    assert validators.validate_network(G) is True

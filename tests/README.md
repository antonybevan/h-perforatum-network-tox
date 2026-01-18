# Unit Tests

Unit tests for the network pharmacology package.

## Running Tests

```bash
# Run all tests
pytest tests/

# Run with coverage
pytest tests/ --cov=network_tox --cov-report=html
```

## Test Structure

- `test_network.py`: Network operation tests
- `test_proximity.py`: Proximity metric tests
- `test_permutation.py`: Permutation testing tests
- `test_rwr_extended.py`: RWR algorithm tests

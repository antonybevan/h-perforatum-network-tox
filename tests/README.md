# Unit Tests

This directory contains unit tests for the network pharmacology package.

## Running Tests

```bash
# Install pytest
pip install pytest pytest-cov

# Run all tests
pytest tests/

# Run with coverage
pytest tests/ --cov=network_tox --cov-report=html

# Run specific test file
pytest tests/test_network.py -v
```

## Test Structure

```
tests/
├── conftest.py           # Pytest configuration
├── test_network.py       # Network operations tests
├── test_proximity.py     # Proximity metrics tests
└── test_permutation.py   # Permutation testing tests
```

## Writing Tests

Follow pytest conventions:
- Test files: `test_*.py`
- Test functions: `def test_*():`
- Use fixtures for common setup
- Assert expected behavior

Example:
```python
def test_my_function():
    result = my_function(input)
    assert result == expected
```

## Current Coverage

Run `pytest --cov` to see test coverage statistics.

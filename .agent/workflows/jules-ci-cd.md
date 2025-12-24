# Jules Workflow: CI/CD and Testing

This document guides Jules (or any AI agent) on handling CI/CD issues.

## Trigger Conditions

Jules should be assigned when:
1. GitHub Actions CI fails
2. pytest tests fail
3. Dependency issues occur

## Common Fixes

### 1. CI Dependency Failures
```bash
# Check requirements.txt for invalid versions
# Use version ranges, not exact pins
# Example: pandas>=2.0.0,<3.0.0 (not pandas==2.3.3)
```

### 2. Git LFS Issues
```yaml
# Ensure checkout has lfs: true
- uses: actions/checkout@v4
  with:
    lfs: true
```

### 3. Unicode Encoding Errors (Windows)
- Replace `≥` with `>=`
- Replace `✓` with `[OK]`
- Replace `✗` with `[ERR]`
- Replace `⚠️` with `[WARN]`

### 4. Test File Not Found
- Check if data files are tracked with Git LFS
- Verify paths use `Path()` not hardcoded strings
- Use `pytest.mark.skipif` for optional data files

## Test Commands

```bash
# Run all tests
pytest tests/ -v

# Run with coverage
pytest tests/ -v --cov=src

# Run specific test
pytest tests/test_rwr.py -v
```

## Post-Fix Checklist

- [ ] All tests pass locally
- [ ] Commit follows conventional format
- [ ] Push triggers green CI

# Contributing to H. Perforatum Network Tox

## Environment Setup

1. **Python 3.10+** required
2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   # Or use conda:
   conda env create -f environment.yml
   ```

## Running Tests

```bash
# Run all tests
pytest tests/ -v

# Run with coverage
pytest tests/ -v --cov=src

# Run specific test file
pytest tests/test_rwr.py -v
```

## Running the Pipeline

```bash
# Full pipeline (data regeneration + validation)
python scripts/run_complete_pipeline.py

# Skip data regeneration (faster)
python scripts/run_complete_pipeline.py --skip-data

# Quick mode (essential steps only)
python scripts/run_complete_pipeline.py --quick
```

## For AI Agents (Jules/Antigravity)

If you are an autonomous agent working on this repository:

1. **Setup**: Install packages from `requirements.txt`
2. **Verify**: Run `pytest tests/ -v` - exit code 0 = success
3. **Style**: Use `ruff check` for linting
4. **Workflow**: See `.agent/workflows/jules-ci-cd.md` for CI/CD fixes

## Scientific Integrity Rules (DO NOT IGNORE)

The following files are **READ-ONLY** unless explicitly instructed:

1. **`docs/METHODOLOGY.md`**: Defines the scientific standard. Code must fit methodology, not vice versa.

2. **`results/`**: Ground truth from verified experiments. Never overwrite without confirmation.

3. **`data/processed/`**: Validated outputs. Use regeneration scripts, don't edit directly.

**Violation of these rules will result in reverted changes.**

## Current Results (Verified 2025-12-24)

| Compound | RWR Z-Score | Significant | Per-Target |
|----------|-------------|-------------|------------|
| Hyperforin | +9.50 | Yes (p<0.0001) | 0.0287 |
| Quercetin | +1.04 | No (p=0.15) | 0.00036 |

**Key Finding:** Hyperforin is ~80x more influential per target than Quercetin.

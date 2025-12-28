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
```

## Running the Pipeline

```bash
# Step 1: Data preprocessing
python scripts/create_lcc_filtered_data.py

# Step 2: Tier 2 analysis (Standard RWI)
python scripts/run_standard_rwr_lcc_permutations.py

# Step 3: Tier 3 analysis (Expression-Weighted EWI)
python scripts/run_expression_weighted_rwr_permutations.py

# Step 4: Negative control (optional)
python scripts/run_chemical_similarity_control.py
```

## For AI Agents

If you are an autonomous agent working on this repository:

1. **Setup**: Install packages from `requirements.txt`
2. **Verify**: Run `pytest tests/ -v` - exit code 0 = success
3. **Style**: Use `ruff check` for linting

## Scientific Integrity Rules (DO NOT IGNORE)

The following files are **READ-ONLY** unless explicitly instructed:

1. **`docs/RESEARCH_SUMMARY.md`**: Defines the scientific standard
2. **`results/tables/`**: Ground truth from verified experiments
3. **`data/processed/`**: Validated outputs

**Violation of these rules will result in reverted changes.**

## Current Results (Verified 2025-12-28)

| Metric | Hyperforin (9 targets) | Quercetin (62 targets) | Ratio |
|--------|------------------------|------------------------|-------|
| **RWI Z-score** | +8.83 | +4.42 | — |
| **EWI Z-score** | +7.99 | +5.56 | — |
| **PTNI (RWI)** | 0.01135 | 0.00052 | **21.9×** |
| **PTNI (EWI)** | 0.0134 | 0.00080 | **16.9×** |

**Key Finding:** Hyperforin is 17–22× more influential per target than Quercetin.

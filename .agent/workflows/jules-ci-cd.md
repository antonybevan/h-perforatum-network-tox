---
description: Jules workflow for handling CI/CD failures and test issues
---

# Jules CI/CD Troubleshooting Workflow

## When to Use
- GitHub Actions CI fails
- pytest tests fail
- Dependency issues occur

## Step 1: Check Test Output
// turbo
```bash
pytest tests/ -v --tb=short
```

## Step 2: Common Fix Categories

### A. Dependency Failures
```bash
# Check requirements.txt for invalid versions
# Use version ranges: pandas>=2.0.0 (not ==2.3.3)
```

Fix location: `requirements.txt`

### B. Git LFS Issues
Ensure checkout has LFS enabled in `.github/workflows/tests.yml`:
```yaml
- uses: actions/checkout@v4
  with:
    lfs: true
```

### C. Unicode Encoding Errors (Windows)
Replace special characters in output strings:
- `≥` → `>=`
- `✓` → `PASS`
- `✗` → `FAIL`
- `×` → `x`

### D. Missing Data Files
Data files should be in:
- `data/processed/` for LCC networks and targets
- `results/tables/` for result CSVs

Key files needed for tests:
- `data/processed/targets_lcc.csv`
- `data/processed/network_900_liver_lcc.parquet`

### E. Deprecated Script References
After project cleanup, these scripts are archived:
- `run_complete_pipeline.py` → archived
- `final_validation_check.py` → archived

Current pipeline scripts (in `scripts/`):
1. `create_lcc_filtered_data.py`
2. `run_standard_rwr_lcc_permutations.py`
3. `run_expression_weighted_rwr_permutations.py`
4. `run_chemical_similarity_control.py`

## Step 3: Run Tests
// turbo
```bash
pytest tests/ -v
```

## Step 4: Verify CI Locally
// turbo
```bash
pip install -r requirements.txt
pytest tests/ -v --cov=src
```

## Step 5: Commit and Push
```bash
git add .
git commit -m "fix(ci): [description of fix]"
git push
```

## Post-Fix Checklist
- [ ] All tests pass locally
- [ ] No deprecated script references
- [ ] requirements.txt has valid version ranges
- [ ] CI should turn green

## Project Structure Reference
```
scripts/                    # 8 essential scripts
results/tables/             # 4 primary result files
data/processed/             # LCC-filtered networks and targets
archive/                    # Deprecated files (don't use in tests)
```

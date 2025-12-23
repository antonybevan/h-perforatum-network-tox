# Contributing to H. Perforatum Net Tox

## Environment Setup

1.  Input: Python 3.10+
2.  Install dependencies:
    ```bash
    pip install -r requirements.txt
    ```

## Running Tests

To run the test suite, simply execute:

```bash
python scripts/cibuild.py
```

Or manually:

```bash
pytest
```

## For AI Agents (Jules)

If you are an autonomous agent working on this repository, please follow these steps to verify your changes:

1.  **Setup**: Ensure all packages in `requirements.txt` are installed.
2.  **Verify**: Run `python scripts/cibuild.py`. This script will return exit code 0 on success and non-zero on failure.
3.  **Style**: Improve code quality where possible, but prioritize correctness verified by the build script.

## ⛔ Scientific Integrity Rules (DO NOT IGNORE)

To preserve the scientific validity of this research, the following files are **STRICTLY READ-ONLY** for autonomous agents unless explicitly instructed otherwise:

1.  **`METHODOLOGY.md`**: This defines the scientific standard. Do not edit it to "fit" your code. Your code must fit the methodology.
2.  **`results/`**: Existing result files represent ground truth from previous experiments. **NEVER overwrite existing files** in this directory without meaningful confirmation. Always save new experiments with new filenames (e.g., `results/experiment_v2.csv`).

**Violation of these rules will result in your changes being reverted.**

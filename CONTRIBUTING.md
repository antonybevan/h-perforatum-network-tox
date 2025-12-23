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

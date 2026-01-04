"""
Pipeline Step 01: Data Ingestion
================================
Downloads and prepares raw datasets for the ECNP DILI prediction pipeline.

Outputs:
    - results/dilirank_full_smiles.csv
    - data/tox21.csv.gz (if not present)

Usage:
    python pipeline/01_data_ingestion.py
"""
import sys
sys.path.insert(0, str(__file__).replace('pipeline\\01_data_ingestion.py', ''))

from src.config import DILIRANK_FULL, TOX21_DATA, DATA_DIR
import pandas as pd
import requests
import gzip
import shutil
from pathlib import Path

print("="*60)
print("STEP 01: DATA INGESTION")
print("="*60)

# --- 1. Download Tox21 if not present ---
if not TOX21_DATA.exists():
    print(f"Downloading Tox21 dataset...")
    urls = [
        "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/tox21.csv.gz",
        "http://bioinf.jku.at/research/DeepTox/tox21_data.csv.gz"
    ]
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    for url in urls:
        try:
            r = requests.get(url, timeout=30)
            r.raise_for_status()
            with open(TOX21_DATA, 'wb') as f:
                f.write(r.content)
            print(f"  Downloaded from {url}")
            break
        except Exception as e:
            print(f"  Failed from {url}: {e}")
else:
    print(f"Tox21 already present: {TOX21_DATA}")

# --- 2. Validate DILIrank ---
if DILIRANK_FULL.exists():
    df = pd.read_csv(DILIRANK_FULL)
    print(f"DILIrank loaded: {len(df)} drugs")
    print(f"  DILI+: {df['is_dili'].sum()}, DILI-: {(df['is_dili']==0).sum()}")
else:
    print(f"ERROR: DILIrank not found: {DILIRANK_FULL}")
    print("  Please run the upstream curation scripts first.")

print("\n[STEP 01 COMPLETE]")

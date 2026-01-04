
import pandas as pd
import gzip
import os

DATA_PATH = r"v:\new\h-perforatum-network-tox\research\ecnp-generalization\data\external\tox21\tox21.csv.gz"

def inspect_tox21():
    print(f"Inspecting: {DATA_PATH}")
    if not os.path.exists(DATA_PATH):
        print("File not found.")
        return

    try:
        # Load first 5 rows
        df = pd.read_csv(DATA_PATH, compression='gzip', nrows=5)
        print("\n--- Columns ---")
        for col in df.columns:
            print(f"- {col}")
        
        print("\n--- First 5 Rows ---")
        print(df.head())
        
        print("\n--- Shape Check ---")
        # Estimate full size? No, just check if it loads.
        print("Load successful.")
        
    except Exception as e:
        print(f"Error reading file: {e}")

if __name__ == "__main__":
    inspect_tox21()


import os
import requests
import gzip
import shutil

# Define URLs for Tox21 Challenge Data (Source: JKU / DeepChem / NCATS)
# Using direct links to the challenge dataset. 
# Updated URL strategy: DeepChem AWS S3 bucket is very reliable.

BASE_URLS = [
    "https://bioinf.jku.at/research/DeepTox/tox21",
    "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets"
]

FILES = [
    "tox21.csv.gz", # Main consolidated file often found on DeepChem
    "tox21_dense_train.csv.gz",
    "tox21_labels_train.csv.gz",
    "tox21_dense_test.csv.gz",
    "tox21_labels_test.csv.gz"
]

OUTPUT_DIR = r"v:\new\h-perforatum-network-tox\research\ecnp-generalization\data\external\tox21"

def download_file(url, dest_path):
    print(f"Attempting download from {url} to {dest_path}...")
    try:
        response = requests.get(url, stream=True, timeout=30)
        if response.status_code == 200:
            with open(dest_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            print(f"SUCCESS: Downloaded {os.path.basename(dest_path)}")
            return True
        else:
            print(f"FAILED: HTTP {response.status_code}")
            return False
    except Exception as e:
        print(f"ERROR: {e}")
        return False

def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
        print(f"Created directory: {OUTPUT_DIR}")

    print(f"Target Directory: {os.path.abspath(OUTPUT_DIR)}")

    for filename in FILES:
        dest_path = os.path.join(OUTPUT_DIR, filename)
        
        # Try primary URL (JKU)
        url1 = f"{BASE_URLS[0]}/{filename}"
        if download_file(url1, dest_path):
            continue
            
        # Try secondary URL (DeepChem) - Filename might vary slightly so we handle the main csv
        if filename == "tox21.csv.gz":
             url2 = f"{BASE_URLS[1]}/{filename}"
             download_file(url2, dest_path)

    # Verify downloads
    print("\n--- Final Directory Listing ---")
    if os.path.exists(OUTPUT_DIR):
        files = os.listdir(OUTPUT_DIR)
        if not files:
            print("Directory is EMPTY.")
        for f in files:
            size_kb = os.path.getsize(os.path.join(OUTPUT_DIR, f)) / 1024
            print(f"- {f}: {size_kb:.2f} KB")
    else:
        print("Directory does not exist!")

if __name__ == "__main__":
    main()

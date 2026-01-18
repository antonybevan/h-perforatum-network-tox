# External Data Sources

This directory contains external data that must be downloaded manually due to size/licensing.

## Required Files

### STRING Database (v12.0)

Download from: https://string-db.org/cgi/download

| File | Download Link |
|------|---------------|
| `string_links.txt.gz` | [9606.protein.links.v12.0.txt.gz](https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz) |
| `string_info.txt.gz` | [9606.protein.info.v12.0.txt.gz](https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz) |

```bash
# Download commands (Linux/macOS)
curl -O https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz
curl -O https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz

# Rename to expected names (Linux/macOS)
mv 9606.protein.links.v12.0.txt.gz string_links.txt.gz
mv 9606.protein.info.v12.0.txt.gz string_info.txt.gz
```

```powershell
# Download commands (Windows PowerShell)
Invoke-WebRequest -Uri "https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz" -OutFile "9606.protein.links.v12.0.txt.gz"
Invoke-WebRequest -Uri "https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz" -OutFile "9606.protein.info.v12.0.txt.gz"

# Rename to expected names (Windows)
Move-Item 9606.protein.links.v12.0.txt.gz string_links.txt.gz
Move-Item 9606.protein.info.v12.0.txt.gz string_info.txt.gz
```

### DILIrank 2.0

Already included: `DILIrank_2.0.xlsx`

Source: [FDA LTKB](https://www.fda.gov/science-research/liver-toxicity-knowledge-base-ltkb)

## After Download

Run the network extraction:
```bash
python scripts/extract_string_network.py
python scripts/create_lcc_filtered_data.py
```

## Note

The processed `.parquet` network files are tracked with Git LFS, so **you don't need to download STRING data** unless you want to regenerate from scratch.

# Corrupted Files Investigation Summary

## Issue Description
Multiple CSV and TXT files in the repository were found to be unreadable using standard tools (`cat`, `more`, `nano`). Investigation revealed these files contain only null bytes (zeros) despite showing correct file sizes.

## Investigation Findings

### Scope of Corruption
- **Total corrupted files**: 135
- **File types affected**: `.csv` and `.txt` files only
- **Unaffected**: `.rdata` files (gzip compressed, intact)

### Corruption Distribution by Directory

| Directory | Corrupted Files |
|-----------|-----------------|
| KETO_PNS/Tables | 80 |
| KETO_CNS/Tables | 33 |
| KETO_CNS/Hippo_GSEA | 11 |
| KETO_PNS/STRING_cluster | 3 |
| KETO_CNS/ProcessedData | 3 |
| KETO_PNS_2/ProcessedData | 2 |
| KETO_PNS/ProcessedData | 2 |

### Key Observations
1. Files show **correct file sizes** but contain only null bytes (`0x00`)
2. Corruption exists in **git blob storage** - files were committed already corrupted
3. The `new_results_all.rdata` file (378MB) is **intact and functional**
4. Some CSV files (e.g., `*_anno.csv`) are valid and readable
5. All newly generated files in `KETO_PNS_2/Tables/` are valid

## Root Cause Analysis

### Most Likely Cause: Cloud Storage Placeholder Files

The symptoms strongly suggest files were copied from a **cloud-synced folder** (OneDrive, Dropbox, Google Drive) where files were set to **"online-only"** or **"files on-demand"**:

- Cloud storage creates "placeholder" files that report correct size but aren't actually downloaded
- When placeholders are copied, the resulting files have correct size but contain null bytes
- Large binary files (`.rdata`) were likely marked "always keep on device" and remained intact
- Smaller CSV/TXT files were "online-only" and became corrupted placeholders

### Alternative Possibilities
1. **Sparse file copy issue** - files were sparse and copied incorrectly via WSL
2. **Interrupted file transfer** - partial write with preallocated file sizes
3. **Storage quota exceeded** during transfer operation

## Impact Assessment

### What Still Works ✅
- `new_results_all.rdata` - main data source containing all DESeq2 objects
- All analysis scripts in `KETO_PNS_2/Code/`
- All newly generated output files in `KETO_PNS_2/Tables/`
- Complete PNS analysis pipeline

### What Is Affected ❌
- Original intermediate CSV/TXT files in KETO_CNS and KETO_PNS
- `KETO2_updated_samples.txt` - sample metadata
- `count_data.csv` - raw count matrix
- Various KEGG, correlation, and VennDetail output tables

## Resolution

### Recommended Actions
1. **Re-copy original files** from source location
2. **Before copying from cloud storage**: Right-click files → "Always keep on this device" (OneDrive) or equivalent
3. **Verify after copying**:
   ```bash
   head -c 100 filename.csv  # Should show actual text content
   file filename.csv         # Should show "ASCII text" or "UTF-8 text"
   ```

### Why Analysis Still Worked
The analysis pipeline loaded data from `new_results_all.rdata`, which contains pre-processed DESeq2 objects:
- `scndds`, `gasdds`, `hipdds` - DESeq2 dataset objects
- `scns`, `gastrocs`, `hippos` - DEG result tables
- `groups` - sample metadata
- `counts` - count matrix

These objects were sufficient for all downstream analyses (DEG, pathway enrichment, Venn diagrams, cross-tissue comparisons).

## Verification Commands

```bash
# Check if a file is corrupted (shows null bytes)
hexdump -C filename.csv | head -5

# Find all corrupted files
find . -type f \( -name "*.csv" -o -name "*.txt" \) \
  -exec sh -c 'head -c 16 "$1" | xxd | grep -q "0000 0000 0000 0000" && echo "$1"' _ {} \;

# Verify .rdata file is intact
file new_results_all.rdata  # Should show "gzip compressed data"
```

## Date
Investigation completed: December 10, 2025
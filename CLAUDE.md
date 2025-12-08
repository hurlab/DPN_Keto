# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a research project analyzing the effects of a ketogenic diet on the nervous system. The project contains RNA-seq data analysis for both central nervous system (CNS - hippocampus) and peripheral nervous system (PNS - sciatic nerves and gastrocnemius muscle) tissues.

## Repository Structure

```
DPN_KETO/
├── KETO_CNS/                # CNS analysis (hippocampus data)
│   ├── Code/               # R analysis scripts for CNS
│   ├── ProcessedData/      # Processed CNS data files
│   ├── Tables/             # CNS result tables
│   └── Figures/            # Generated figures
├── KETO_PNS/               # Original PNS data (incomplete scripts)
│   ├── Code/               # Limited R scripts
│   ├── ProcessedData/      # Raw PNS data
│   ├── new_results_all.rdata  # Main PNS RData file
│   └── Tables/             # Some existing tables
├── KETO_PNS_2/             # Reconstructed PNS analysis (COMPLETE)
│   ├── Code/               # Complete R analysis scripts
│   ├── ProcessedData/      # Input data files
│   ├── Tables/             # Output tables by analysis type
│   └── Figures/            # Generated figures
└── KETO_MANU/              # Manuscript files
```

## Common Commands

### Running CNS Analysis
```bash
cd KETO_CNS/Code
Rscript run_Pong_final_new_paper2.R
```

### Running Complete PNS Analysis (RECONSTRUCTED)
```bash
cd KETO_PNS_2/Code
Rscript 00_Run_All_Analyses.R
```

### Running Individual PNS Analyses
```bash
cd KETO_PNS_2/Code
Rscript 01_DEG_analysis.R          # Differential expression
Rscript 02_Enrichment_analysis.R   # KEGG/GSEA
Rscript 03_Mfuzz_clustering.R      # Soft clustering
Rscript 04_Correlation_analysis.R  # Correlations
Rscript 05_Visualization.R         # Figures
```

## Analysis Pipeline Architecture

### Data Flow
1. **Input**: RNA-seq count data in `new_results_all.rdata` (378MB)
2. **Preprocessing**: Gene symbol annotation, duplicate handling
3. **DESeq2 Analysis**: Differential expression by tissue (SCN, Gastroc, Hippo)
4. **Downstream Analysis**:
   - PCA and sample clustering
   - KEGG pathway enrichment
   - Gene Set Enrichment Analysis (GSEA)
   - Mfuzz soft clustering
   - Inter-tissue correlations

### Key Comparisons
- Diet effects: HFD vs SD, KD vs SD, KD vs HFD
- Interventions: DR vs HFD, KDI vs HFD, EX vs HFD, KDI_EX vs HFD

### Tissue Types
- **SCN**: Sciatic nerve
- **Gastroc**: Gastrocnemius muscle
- **Hippo**: Hippocampus (CNS reference)

## Important Implementation Details

### Color Scheme
The analysis uses a consistent color scheme across all visualizations:
- SD: #DAD9D9 (light gray)
- HFD: #774099 (purple)
- KD: #A5C3DE (light blue)
- DR: #7d7d7d (gray)
- KDI: #3E7DB1 (blue)
- EX: #F9E500 (yellow)
- KDI_EX: #93C954 (green)

### Statistical Thresholds
- DEG significance: p-value < 0.01, log2FC > 1
- Mfuzz clustering: Top 25% most variable genes
- KEGG enrichment: Built-in mouse annotation (mmko)

### Dependencies
- **CRAN packages**: tidyverse, pheatmap, VennDetail, corrplot, Hmisc, ggrepel
- **Bioconductor**: DESeq2, richR, Mfuzz, Biobase, org.Mm.eg.db

## Working with This Repository

### For Analyses
1. Use `KETO_PNS_2` for complete PNS analysis - it contains fully reconstructed scripts
2. Reference `KETO_CNS/Code` for established analysis patterns
3. All scripts are designed to be run from their respective Code directories

### Data Access
- Main data object: `load("KETO_PNS/ProcessedData/new_results_all.rdata")`
- Contains DESeq2 objects: `scndds`, `gasdds`, `hipdds`
- Sample metadata: `groups` dataframe
- Count matrix: `counts`

### Common Tasks
- **Add new tissues**: Extend `create_dds()` function in DEG analysis
- **Custom comparisons**: Modify `make.comparison()` function parameters
- **New visualizations**: Add to `05_Visualization.R` following existing patterns

## Notes

- CIBERSORTx analysis mentioned in manuscript is not included (requires external tool)
- PNS original scripts in `KETO_PNS/Code` are incomplete - use `KETO_PNS_2` instead
- Manuscript files are in `KETO_MANU/` for reference
- All output files are automatically saved with descriptive names
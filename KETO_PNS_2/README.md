# KETO_PNS_2: Reconstructed PNS Analysis Scripts

This folder contains the reconstructed R analysis scripts for the peripheral nervous system (PNS) data from the KETO study. These scripts reproduce the analyses performed for the CNS data, adapted for the PNS tissues (sciatic nerve and gastrocnemius muscle).

## Folder Structure

```
KETO_PNS_2/
├── Code/                    # Analysis scripts
│   ├── 00_Run_All_Analyses.R     # Master script to run all analyses
│   ├── 01_DEG_analysis.R         # Differential expression analysis
│   ├── 02_Enrichment_analysis.R  # KEGG/GO enrichment and GSEA
│   ├── 03_Mfuzz_clustering.R     # Soft clustering of expression patterns
│   ├── 04_Correlation_analysis.R # Inter-tissue and sample correlations
│   ├── 05_Visualization.R        # Additional plots and figures
│   └── examine_data.R            # Data exploration script
├── ProcessedData/           # Input data files
│   ├── new_results_all.rdata     # Main RData file with all objects
│   └── *.txt, *.csv              # Sample information and count data
├── Tables/                  # Output tables
│   ├── DEGs/                   # Differential expression results
│   ├── KEGG/                   # KEGG enrichment results
│   ├── Mfuzz/                  # Clustering results
│   └── Correlation/            # Correlation analysis results
└── Figures/                 # Generated figures
```

## Prerequisites

### Required R Packages

Install the following packages before running the analysis:

```r
# CRAN packages
install.packages(c("tidyverse", "pheatmap", "VennDetail", "corrplot",
                   "Hmisc", "ggrepel", "gridExtra", "cowplot",
                   "circlize", "jsonlite"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("DESeq2", "richR", "org.Mm.eg.db", "Mfuzz",
                       "Biobase", "complexheatmap", "enrichplot",
                       "DOSE", "clusterProfiler", "enrichplot"))
```

## Running the Analysis

### Option 1: Run All Analyses at Once

```r
# Navigate to the KETO_PNS_2 directory
setwd("/path/to/KETO_PNS_2")

# Run the master script
source("Code/00_Run_All_Analyses.R")
```

### Option 2: Run Individual Analyses

```r
# Navigate to the KETO_PNS_2 directory
setwd("/path/to/KETO_PNS_2")

# 1. Explore data structure
source("Code/examine_data.R")

# 2. Differential expression analysis
source("Code/01_DEG_analysis.R")

# 3. Functional enrichment analysis
source("Code/02_Enrichment_analysis.R")

# 4. Mfuzz clustering
source("Code/03_Mfuzz_clustering.R")

# 5. Correlation analysis
source("Code/04_Correlation_analysis.R")

# 6. Visualizations
source("Code/05_Visualization.R")
```

## Analysis Workflow

1. **Data Exploration**: Examine the structure of the RData file to understand available objects
2. **Differential Expression**: DESeq2 analysis for each tissue (SCN, Gastroc, Hippo)
3. **Functional Enrichment**: KEGG pathway analysis and Gene Set Enrichment Analysis (GSEA)
4. **Mfuzz Clustering**: Soft clustering to identify co-expressed gene patterns
5. **Correlation Analysis**: Inter-tissue and inter-sample correlations
6. **Visualization**: Comprehensive figure generation including PCA, heatmaps, volcano plots, etc.

## Key Comparisons

The analysis includes the following key comparisons:

### Diet Comparisons
- HFD vs SD (High-fat diet vs Standard diet)
- KD vs SD (Ketogenic diet vs Standard diet)
- KD vs HFD (Ketogenic diet vs High-fat diet)

### Intervention Comparisons
- DR vs HFD (Dietary restriction vs High-fat diet)
- KDI vs HFD (Ketogenic diet intervention vs High-fat diet)
- EX vs HFD (Exercise vs High-fat diet)
- KDI_EX vs HFD (Combined intervention vs High-fat diet)

## Output Files

### Tables
- **DEGs**: Lists of differentially expressed genes for each comparison
- **KEGG**: KEGG pathway enrichment results and GSEA results
- **Mfuzz**: Gene cluster assignments and pathway enrichment for clusters
- **Correlation**: Inter-tissue and sample correlation matrices

### Figures
- PCA plots with and without confidence ellipses
- Sample distance heatmaps
- DEG heatmaps and volcano plots
- KEGG enrichment dot plots
- Mfuzz clustering visualizations
- Correlation plots and heatmaps

## Notes

1. The analysis assumes the RData file contains DESeq2 objects for each tissue (scndds, gasdds, hipdds)
2. Some analyses (like CIBERSORTx mentioned in the original request) are not included as they require external tools
3. The color scheme matches the CNS analysis for consistency
4. All scripts include error handling and will skip analyses if required data is not available

## Analysis Time

The complete analysis pipeline typically takes:
- DEG analysis: 5-10 minutes
- Enrichment analysis: 10-20 minutes
- Mfuzz clustering: 5-10 minutes per tissue
- Correlation analysis: 2-5 minutes
- Visualization: 5-10 minutes

Total runtime: ~30-60 minutes depending on computational resources

## Reproducibility

- Session information is saved to `Tables/R_session_info.csv`
- A complete analysis summary is generated as `KETO_PNS_Analysis_Summary.json` and `.txt`
- All figures and tables are saved with clear, descriptive names

## Questions or Issues

If you encounter any issues with the analysis:
1. Check that all required packages are installed
2. Verify the input data files are present in `ProcessedData/`
3. Ensure you have sufficient disk space for outputs (~1-2 GB)
4. Check R version compatibility (tested with R 4.0+)
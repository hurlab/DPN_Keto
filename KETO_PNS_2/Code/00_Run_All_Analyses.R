# ========================================
# Master Script for PNS Data Analysis
# ========================================
# This script runs all PNS analyses in sequence
# Based on the CNS analysis pipeline

# Clear workspace
rm(list = ls())

# Load required libraries (check availability first)
required_packages <- c("tidyverse", "DESeq2", "pheatmap", "VennDetail",
                      "richR", "org.Mm.eg.db", "Mfuzz", "Biobase",
                      "corrplot", "Hmisc", "ggrepel", "gridExtra",
                      "cowplot", "complexheatmap", "circlize")

missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if(length(missing_packages) > 0) {
  cat("The following packages are missing:\n")
  cat(paste(missing_packages, collapse = "\n"), "\n")
  cat("\nPlease install them using:\n")
  cat(paste("install.packages('", missing_packages, "')", collapse = "\n"), "\n")
  cat("\nSome packages may need Bioconductor:\n")
  cat("if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')\n")
  cat(paste("BiocManager::install('", missing_packages, "')", collapse = "\n"), "\n")
  stop("Missing required packages")
}

cat("All required packages are available!\n")
cat("Starting PNS data analysis pipeline...\n\n")

# Record start time
start_time <- Sys.time()

# Set working directory
setwd("..")

# ========================================
# 1. Data Exploration
# ========================================

cat("1. Exploring data structure...\n")
source("Code/examine_data.R")

# ========================================
# 2. Differential Expression Analysis
# ========================================

cat("\n2. Performing differential expression analysis...\n")
cat("   This may take several minutes...\n")
source("Code/01_DEG_analysis.R")

# ========================================
# 3. Functional Enrichment Analysis
# ========================================

cat("\n3. Performing functional enrichment analysis...\n")
cat("   KEGG and GSEA analyses...\n")
source("Code/02_Enrichment_analysis.R")

# ========================================
# 4. Mfuzz Clustering Analysis
# ========================================

cat("\n4. Performing Mfuzz clustering analysis...\n")
cat("   This may take several minutes...\n")
source("Code/03_Mfuzz_clustering.R")

# ========================================
# 5. Correlation Analysis
# ========================================

cat("\n5. Performing correlation analysis...\n")
source("Code/04_Correlation_analysis.R")

# ========================================
# 6. Additional Visualizations
# ========================================

cat("\n6. Creating additional visualizations...\n")
source("Code/05_Visualization.R")

# ========================================
# Summary and Cleanup
# ========================================

# Calculate elapsed time
end_time <- Sys.time()
elapsed_time <- difftime(end_time, start_time, units = "hours")

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("PNS DATA ANALYSIS PIPELINE COMPLETED!\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("Total elapsed time:", round(as.numeric(elapsed_time), 2), "hours\n")
cat("\nGenerated outputs:\n")
cat("- DEG tables: Tables/DEGs/\n")
cat("- Enrichment results: Tables/KEGG/\n")
cat("- Clustering results: Tables/Mfuzz/\n")
cat("- Correlation results: Tables/Correlation/\n")
cat("- Figures: Figures/\n")

# Create a summary report
summary_report <- list(
  Analysis_Date = as.character(start_time),
  Completion_Time = as.character(end_time),
  Elapsed_Time_Hours = as.numeric(elapsed_time),
  R_Version = R.version.string,
  Packages_Used = required_packages,
  Output_Directories = c("Tables/DEGs", "Tables/KEGG", "Tables/Mfuzz",
                         "Tables/Correlation", "Figures")
)

# Save summary report
summary_report$json <- jsonlite::toJSON(summary_report, auto_unbox = TRUE, pretty = TRUE)
writeLines(summary_report$json, "KETO_PNS_Analysis_Summary.json")

# Also create a text summary
sink("KETO_PNS_Analysis_Summary.txt")
cat("KETO PNS Data Analysis Summary\n")
cat(paste(rep("=", 50), collapse = ""), "\n")
cat("Start time:", as.character(start_time), "\n")
cat("Completion time:", as.character(end_time), "\n")
cat("Total elapsed time:", round(as.numeric(elapsed_time), 2), "hours\n")
cat("\nAnalyses performed:\n")
cat("1. Data exploration and structure analysis\n")
cat("2. Differential expression analysis (DESeq2)\n")
cat("3. KEGG and GO enrichment analysis\n")
cat("4. Gene Set Enrichment Analysis (GSEA)\n")
cat("5. Mfuzz soft clustering\n")
cat("6. Inter-tissue and inter-sample correlations\n")
cat("7. Comprehensive visualization generation\n")
cat("\nOutput locations:\n")
cat("- DEG results: Tables/DEGs/\n")
cat("- Enrichment results: Tables/KEGG/\n")
cat("- Mfuzz clusters: Tables/Mfuzz/\n")
cat("- Correlation analysis: Tables/Correlation/\n")
cat("- All figures: Figures/\n")
cat("\nR version:", R.version.string, "\n")
sink()

cat("\nAnalysis summary saved to:\n")
cat("- KETO_PNS_Analysis_Summary.json\n")
cat("- KETO_PNS_Analysis_Summary.txt\n")

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("ANALYSIS PIPELINE FINISHED SUCCESSFULLY!\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
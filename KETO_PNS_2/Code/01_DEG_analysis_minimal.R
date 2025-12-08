# ========================================
# PNS Differential Expression Gene Analysis (Minimal)
# ========================================
# This script reproduces the main DEG analysis for PNS data
# Using only base R and available packages

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(pheatmap)
})

# Set working directory - we're already in KETO_PNS_2
# No need to change directory
if(!dir.exists("Figures")) dir.create("Figures")
if(!dir.exists("Tables/DEGs")) dir.create("Tables/DEGs")

# Load the data
load("ProcessedData/new_results_all.rdata")

# Clean up groups data if needed
groups$Intervention <- sub('\\-', '_', groups$Intervention)

# Define comparison groups
dietgroup <- c("HFDvsSD", "KDvsSD", "KDvsHFD")
integroup <- c("DRvsHFD", "KDIvsHFD", "EXvsHFD", "KDI_EXvsHFD")

# Define color scheme
samplecolor <- c("#DAD9D9", "#774099", "#A5C3DE", "#7d7d7d", "#3E7DB1", "#F9E500", "#93C954")
names(samplecolor) <- c("SD", "HFD", "KD", "DR", "KDI", "EX", "KDI_EX")

# ========================================
# PCA Analysis
# ========================================

cat("Generating PCA plots...\n")

# Function to generate PCA plots
generate_pca <- function(dds, tissue, samplecolor) {
  vst <- varianceStabilizingTransformation(dds)
  pcaData <- plotPCA(vst, returnData = TRUE)
  pcaData$group <- sub(paste0(tissue, '_'), '', pcaData$group)

  # Calculate variance percentages
  percentVar <- round(100 * attr(pcaData, "percentVar"))

  # Create basic plot using base R
  pdf(paste0("Figures/PCA_", tissue, "_basic.pdf"))
  plot(pcaData$PC1, pcaData$PC2,
       col = samplecolor[pcaData$group],
       pch = 19, cex = 2,
       xlab = paste0("PC1: ", percentVar[1], "% variance"),
       ylab = paste0("PC2: ", percentVar[2], "% variance"),
       main = paste0("PCA: ", tissue))
  legend("topright", legend = names(samplecolor),
         col = samplecolor, pch = 19, cex = 0.8)
  dev.off()

  cat("PCA plot for", tissue, "saved\n")
}

# Generate PCA for each tissue
if(exists("scndds")) generate_pca(scndds, "SCN", samplecolor)
if(exists("gasdds")) generate_pca(gasdds, "Gastroc", samplecolor)
if(exists("hipdds")) generate_pca(hipdds, "Hippo", samplecolor)

# ========================================
# Extract DEGs (Already computed)
# ========================================

cat("Extracting DEG results...\n")

# The RData file already contains DEG results
# Just save them in organized format

if(exists("scns")) {
  for(name in names(scns)) {
    write.csv(scns[[name]],
              file = paste0("Tables/DEGs/", name, "_DEG_pval.csv"),
              row.names = TRUE)
  }
  cat("SCN DEGs saved\n")
}

if(exists("gastrocs")) {
  for(name in names(gastrocs)) {
    write.csv(gastrocs[[name]],
              file = paste0("Tables/DEGs/", name, "_DEG_pval.csv"),
              row.names = TRUE)
  }
  cat("Gastroc DEGs saved\n")
}

if(exists("hippos")) {
  for(name in names(hippos)) {
    write.csv(hippos[[name]],
              file = paste0("Tables/DEGs/", name, "_DEG_pval.csv"),
              row.names = TRUE)
  }
  cat("Hippo DEGs saved\n")
}

# ========================================
# Summary Statistics
# ========================================

cat("Generating DEG summary statistics...\n")

# Create summary for each tissue
create_deg_summary <- function(deg_list, tissue) {
  summary_df <- data.frame(
    Comparison = names(deg_list),
    Total_DEGs = sapply(deg_list, nrow),
    Upregulated = sapply(deg_list, function(x) sum(x$log2FoldChange > 0)),
    Downregulated = sapply(deg_list, function(x) sum(x$log2FoldChange < 0))
  )
  write.csv(summary_df, file = paste0("Tables/DEGs/", tissue, "_DEG_summary.csv"),
            row.names = FALSE)
  return(summary_df)
}

if(exists("scns")) {
  scn_summary <- create_deg_summary(scns, "SCN")
  print(scn_summary)
}

if(exists("gastrocs")) {
  gastroc_summary <- create_deg_summary(gastrocs, "Gastroc")
  print(gastroc_summary)
}

if(exists("hippos")) {
  hippo_summary <- create_deg_summary(hippos, "Hippo")
  print(hippo_summary)
}

# ========================================
# Sample Distance Heatmap
# ========================================

cat("Creating sample distance heatmap...\n")

if(exists("scndds")) {
  vst_data <- assay(vst(scndds))
  sample_dist <- dist(t(vst_data))
  sample_dist_matrix <- as.matrix(sample_dist)

  # Create heatmap
  pdf("Figures/sample_distance_heatmap.pdf", width = 10, height = 10)
  pheatmap(sample_dist_matrix,
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation",
           main = "Sample Distance Heatmap (SCN)")
  dev.off()

  cat("Sample distance heatmap saved\n")
}

cat("\nDEG analysis completed successfully!\n")
cat("Results saved in Tables/DEGs/ and Figures/\n")
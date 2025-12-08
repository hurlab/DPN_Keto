# Test basic data access without complex packages
cat("Testing PNS data access...\n")

# Load the data
load("ProcessedData/new_results_all.rdata")

# Basic data exploration
cat("\n=== Basic Data Summary ===\n")
cat("Count matrix dimensions:", dim(counts), "\n")
cat("Number of samples:", ncol(counts), "\n")
cat("Number of genes:", nrow(counts), "\n")

cat("\n=== Sample Groups ===\n")
print(table(groups$Group))

cat("\n=== Interventions ===\n")
print(table(groups$Intervention))

cat("\n=== Conditions ===\n")
print(table(groups$condition))

# Check if DESeq2 objects exist
cat("\n=== DESeq2 Objects ===\n")
if(exists("scndds")) {
  cat("SCN DESeq2 object: ", class(scndds), "\n")
  cat("Dimensions: ", dim(scndds), "\n")
}
if(exists("gasdds")) {
  cat("Gastroc DESeq2 object: ", class(gasdds), "\n")
  cat("Dimensions: ", dim(gasdds), "\n")
}
if(exists("hipdds")) {
  cat("Hippo DESeq2 object: ", class(hipdds), "\n")
  cat("Dimensions: ", dim(hipdds), "\n")
}

# Check if results exist
cat("\n=== DEG Results ===\n")
if(exists("scns")) {
  cat("SCN DEGs available for", length(scns), "comparisons\n")
}
if(exists("gastrocs")) {
  cat("Gastroc DEGs available for", length(gastrocs), "comparisons\n")
}
if(exists("hippos")) {
  cat("Hippo DEGs available for", length(hippos), "comparisons\n")
}

# Save a simple summary
summary_data <- list(
  total_genes = nrow(counts),
  total_samples = ncol(counts),
  tissue_counts = table(groups$Group),
  intervention_counts = table(groups$Intervention),
  has_scn = exists("scndds"),
  has_gastroc = exists("gasdds"),
  has_hippo = exists("hipdds")
)

saveRDS(summary_data, file = "ProcessedData/data_summary.rds")
cat("\nData summary saved to ProcessedData/data_summary.rds\n")
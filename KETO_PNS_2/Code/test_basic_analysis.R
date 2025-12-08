# Test basic analysis with available packages
cat("Testing basic PNS analysis with available packages...\n")

# Load what we can
library(DESeq2)
library(pheatmap)

# Load the data
load("ProcessedData/new_results_all.rdata")

# Basic summary
cat("\n=== PNS Data Summary ===\n")
cat("Total genes:", nrow(counts), "\n")
cat("Total samples:", ncol(counts), "\n")
cat("\nSample distribution by tissue:\n")
print(table(groups$Group))
cat("\nSample distribution by intervention:\n")
print(table(groups$Intervention))

# Test working with DESeq2 objects
if(exists("scndds")) {
  cat("\n=== Sciatic Nerve (SCN) Analysis ===\n")
  cat("DESeq2 object dimensions:", dim(scndds), "\n")

  # Get results names
  results_names <- resultsNames(scndds)
  cat("Available comparisons:", length(results_names), "\n")
  print(head(results_names))
}

# Create a simple sample distance heatmap if we have multiple tissues
if(exists("scndds") && exists("gasdds")) {
  cat("\n=== Creating Sample Distance Plot ===\n")

  # Get variance stabilized data for SCN
  vst_scn <- assay(vst(scndds))
  sample_dist_scn <- dist(t(vst_scn))

  # Simple distance matrix
  dist_matrix <- as.matrix(sample_dist_scn)[1:10, 1:10]  # First 10 samples

  # Create a simple heatmap
  pdf("Figures/sample_distance_test.pdf", width = 8, height = 6)
  pheatmap(dist_matrix,
           main = "Sample Distance Matrix (First 10 Samples)",
           display_numbers = TRUE,
           number_format = "%.2f")
  dev.off()

  cat("Sample distance heatmap saved to Figures/sample_distance_test.pdf\n")
}

cat("\nBasic test completed successfully!\n")
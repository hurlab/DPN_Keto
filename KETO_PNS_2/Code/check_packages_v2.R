# Check packages needed for additional analyses
cat("Checking package availability...\n")

# Required packages for additional analyses
packages_needed <- c(
  "VennDetail",      # For Venn diagrams
  "org.Mm.eg.db",    # For mouse gene annotation
  "richR",           # For KEGG enrichment
  "clusterProfiler", # For GO enrichment
  "enrichplot",      # For visualization
  "circlize",       # For circular plots
  "VennDiagram",    # Alternative Venn diagram package
  "ComplexHeatmap"  # For advanced heatmaps
)

installed <- installed.packages()[, "Package"]
status <- data.frame(
  Package = packages_needed,
  Installed = packages_needed %in% installed,
  stringsAsFactors = FALSE
)

print(status)
cat("\n", sum(status$Installed), "out of", length(packages_needed), "packages installed\n")

# Try to load key packages
test_load <- c("VennDetail", "org.Mm.eg.db", "richR", "clusterProfiler")
for(pkg in test_load) {
  if(pkg %in% installed) {
    tryCatch({
      library(pkg, character.only = TRUE)
      cat("✓", pkg, "loaded successfully\n")
    }, error = function(e) {
      cat("✗", pkg, "failed to load:", e$message, "\n")
    })
  } else {
    cat("✗", pkg, "not installed\n")
  }
}
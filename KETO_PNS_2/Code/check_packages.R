# Check which packages are available
cat("Checking installed packages for PNS analysis...\n")

# List all packages we need
required_packages <- c(
  "DESeq2", "pheatmap", "tidyverse", "ggplot2", "dplyr",
  "VennDetail", "corrplot", "Hmisc", "ggrepel",
  "org.Mm.eg.db", "Mfuzz", "Biobase", "richR",
  "clusterProfiler", "enrichplot", "DOSE"
)

cat("\nPackage Status:\n")
cat("================\n")

installed <- installed.packages()[, "Package"]
success <- 0
total <- length(required_packages)

for(pkg in required_packages) {
  if(pkg %in% installed) {
    cat(paste("✓", pkg, "is installed\n"))
    success <- success + 1
  } else {
    cat(paste("✗", pkg, "is NOT installed\n"))
  }
}

cat(paste("\n", success, "out of", total, "packages installed\n"))

# Test loading key packages
cat("\nTesting key packages:\n")
cat("====================\n")

test_load <- c("DESeq2", "pheatmap", "tidyverse")
for(pkg in test_load) {
  if(pkg %in% installed) {
    tryCatch({
      library(pkg, character.only = TRUE)
      cat(paste("✓", pkg, "loaded successfully\n"))
    }, error = function(e) {
      cat(paste("✗", pkg, "failed to load:", e$message, "\n"))
    })
  }
}
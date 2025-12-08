# Install missing packages for PNS analysis
cat("Installing missing packages...\n")

# Check and install CRAN packages
cran_packages <- c("pheatmap", "ggplot2", "tidyverse", "VennDetail",
                   "corrplot", "Hmisc", "ggrepel", "gridExtra",
                   "cowplot", "circlize", "jsonlite")

for(pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste("Installing", pkg, "from CRAN...\n"))
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  } else {
    cat(paste(pkg, "is already installed\n"))
  }
}

# Check and install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager...\n")
  install.packages("BiocManager", repos = "https://cloud.r-project.org/")
}

bioc_packages <- c("org.Mm.eg.db", "Mfuzz", "Biobase", "richR",
                   "clusterProfiler", "enrichplot", "DOSE")

for(pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste("Installing", pkg, "from Bioconductor...\n"))
    BiocManager::install(pkg, update = FALSE)
  } else {
    cat(paste(pkg, "is already installed\n"))
  }
}

# Test loading all packages
cat("\nTesting package loading...\n")
test_packages <- c(cran_packages, bioc_packages)
success <- 0
for(pkg in test_packages) {
  tryCatch({
    library(pkg, character.only = TRUE)
    cat(paste("✓", pkg, "loaded successfully\n"))
    success <- success + 1
  }, error = function(e) {
    cat(paste("✗", pkg, "failed to load:", e$message, "\n"))
  })
}

cat(paste("\nSuccessfully loaded", success, "out of", length(test_packages), "packages\n"))
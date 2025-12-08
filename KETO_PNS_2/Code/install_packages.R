# Install required packages for PNS analysis
cat("Installing required packages...\n")

# Install BiocManager if not present
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cran.r-project.org")
}

# CRAN packages
cran_packages <- c("tidyverse", "pheatmap", "VennDetail", "corrplot",
                   "Hmisc", "ggrepel", "gridExtra", "cowplot",
                   "circlize", "jsonlite")

for(pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste("Installing", pkg, "...\n"))
    install.packages(pkg, repos = "https://cran.r-project.org")
  }
}

# Bioconductor packages
bioc_packages <- c("DESeq2", "richR", "org.Mm.eg.db", "Mfuzz",
                   "Biobase", "clusterProfiler", "enrichplot",
                   "DOSE")

for(pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste("Installing", pkg, "from Bioconductor...\n"))
    BiocManager::install(pkg, update = FALSE)
  }
}

cat("Package installation completed!\n")
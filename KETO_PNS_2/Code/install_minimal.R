# Install only essential packages for PNS analysis
cat("Installing essential packages for PNS analysis...\n")

# Install BiocManager if not present
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org/")
}

# Minimal essential packages - without tidyverse for now
essential_packages <- c("DESeq2", "pheatmap", "org.Mm.eg.db", "Mfuzz", "Biobase")

cat("Installing essential Bioconductor packages...\n")
for(pkg in essential_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste("Installing", pkg, "...\n"))
    BiocManager::install(pkg, update = FALSE, force = TRUE)
  }
}

# Try installing base R packages separately
base_packages <- c("RColorBrewer", "gplots")
for(pkg in base_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste("Installing", pkg, "...\n"))
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  }
}

cat("Checking installed packages...\n")
installed <- rownames(installed.packages())
for(pkg in c(essential_packages, base_packages)) {
  if(pkg %in% installed) {
    cat(paste("✓", pkg, "is installed\n"))
  } else {
    cat(paste("✗", pkg, "is NOT installed\n"))
  }
}
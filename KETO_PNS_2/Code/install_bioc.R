# Install Bioconductor packages
cat("Installing Bioconductor packages...\n")

# Install BiocManager if not present
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org/")
}

# Install packages
bioc_packages <- c("org.Mm.eg.db", "Mfuzz", "richR",
                   "clusterProfiler", "enrichplot", "DOSE",
                   "VennDetail", "Hmisc")

for(pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste("Installing", pkg, "...\n"))
    BiocManager::install(pkg, update = FALSE)
  } else {
    cat(paste(pkg, "is already installed\n"))
  }
}

cat("\nBioconductor package installation completed!\n")
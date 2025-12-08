# Install Bioconductor 3.20 packages for R 4.4
cat("Installing Bioconductor 3.20 packages...\n")

# Install BiocManager if not present
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org/")
}

# Force Bioconductor 3.20
BiocManager::install(version = "3.20", update = FALSE, ask = FALSE)

# Install key packages
bioc_packages <- c("org.Mm.eg.db", "Mfuzz", "richR",
                   "clusterProfiler", "enrichplot", "DOSE")

for(pkg in bioc_packages) {
  cat(paste("Installing", pkg, "...\n"))
  BiocManager::install(pkg, version = "3.20", update = FALSE, ask = FALSE)
}

# Install CRAN packages
cran_packages <- c("VennDetail", "Hmisc")
for(pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste("Installing", pkg, "from CRAN...\n"))
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  }
}

cat("\nPackage installation completed!\n")
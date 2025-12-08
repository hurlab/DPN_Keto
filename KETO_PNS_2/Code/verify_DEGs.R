# Verify DEG numbers and tissue assignments
cat("Verifying DEG numbers...\n")

# Load data
load("ProcessedData/new_results_all.rdata")

# Check what's actually in scns, gastrocs, hippos
cat("\nChecking DEG objects:\n")
if(exists("scns")) {
  cat("SCN DEGs available for comparisons:", paste(names(scns), collapse=", "), "\n")
  for(name in names(scns)) {
    cat("  ", name, ":", nrow(scns[[name]]), "DEGs\n")
  }
}

if(exists("gastrocs")) {
  cat("\nGastroc DEGs available for comparisons:", paste(names(gastrocs), collapse=", "), "\n")
  for(name in names(gastrocs)) {
    cat("  ", name, ":", nrow(gastrocs[[name]]), "DEGs\n")
  }
}

if(exists("hippos")) {
  cat("\nHippo DEGs available for comparisons:", paste(names(hippos), collapse=", "), "\n")
  for(name in names(hippos)) {
    cat("  ", name, ":", nrow(hippos[[name]]), "DEGs\n")
  }
}

# Summary
cat("\n=== SUMMARY ===\n")
if(exists("scns")) cat("Total SCN DEG objects:", length(scns), "\n")
if(exists("gastrocs")) cat("Total Gastroc DEG objects:", length(gastrocs), "\n")
if(exists("hippos")) cat("Total Hippo DEG objects:", length(hippos), "\n")
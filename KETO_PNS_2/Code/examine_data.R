# Script to examine PNS data structure
# This script will load the RData file and explore what's available

# Load the RData file
load("../ProcessedData/new_results_all.rdata")

# List all objects in the workspace
cat("Objects in RData file:\n")
print(ls())

# Check if key objects exist and examine their structure
if (exists("counts")) {
  cat("\n=== Counts data ===\n")
  print(dim(counts))
  print(head(colnames(counts)))
  print(head(rownames(counts)))
}

if (exists("groups")) {
  cat("\n=== Groups data ===\n")
  print(str(groups))
  print(unique(groups$Group))
  print(unique(groups$Intervention))
}

if (exists("scndds")) {
  cat("\n=== SCN DESeq2 object ===\n")
  print(class(scndds))
  print(dim(scndds))
}

if (exists("gasdds")) {
  cat("\n=== Gastroc DESeq2 object ===\n")
  print(class(gasdds))
  print(dim(gasdds))
}

# Check for other common objects
common_objects <- c("scns", "gastrocs", "scnk", "gastrock", "hippo", "hippos", "hippok")
for (obj in common_objects) {
  if (exists(obj)) {
    cat(sprintf("\n=== %s ===\n", obj))
    print(str(get(obj)))
  }
}
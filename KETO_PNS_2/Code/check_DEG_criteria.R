# Check what criteria were used for DEGs in the original analysis
cat("Checking DEG criteria...\n")

# Load data
load("ProcessedData/new_results_all.rdata")

# Check scns - this should be the significant DEGs
cat("\n=== Checking scns object ===\n")
if(exists("scns")) {
  for(name in names(scns)) {
    deg <- scns[[name]]
    cat("\n", name, ":\n")
    cat("  Total genes in object:", nrow(deg), "\n")
    cat("  Min p-value:", min(deg$pvalue, na.rm=TRUE), "\n")
    cat("  Max p-value:", max(deg$pvalue, na.rm=TRUE), "\n")
    cat("  Genes with p < 0.01:", sum(deg$pvalue < 0.01, na.rm=TRUE), "\n")
    cat("  Genes with padj < 0.01:", sum(deg$padj < 0.01, na.rm=TRUE), "\n")
    cat("  Min |log2FC|:", min(abs(deg$log2FoldChange), na.rm=TRUE), "\n")
    cat("  Genes with |log2FC| > 1:", sum(abs(deg$log2FoldChange) > 1, na.rm=TRUE), "\n")
  }
}

# Also check the full results (scn object)
cat("\n=== Checking scn object (full results) ===\n")
if(exists("scn")) {
  # Check one comparison
  comp_name <- names(scn)[1]
  res <- scn[[comp_name]]
  cat("\n", comp_name, " - Full results:\n")
  cat("  Total genes:", nrow(res), "\n")

  # Apply different filters
  cat("  Genes with p < 0.01:", sum(res$pvalue < 0.01, na.rm=TRUE), "\n")
  cat("  Genes with p < 0.01 AND |log2FC| > 1:", sum(res$pvalue < 0.01 & abs(res$log2FoldChange) > 1, na.rm=TRUE), "\n")
  cat("  Genes with padj < 0.01:", sum(res$padj < 0.01, na.rm=TRUE), "\n")
  cat("  Genes with padj < 0.01 AND |log2FC| > 1:", sum(res$padj < 0.01 & abs(res$log2FoldChange) > 1, na.rm=TRUE), "\n")
}

# Summary
cat("\n=== CRITICAL FINDING ===\n")
cat("The manuscript uses: p < 0.01 (unadjusted) with NO log2FC cutoff\n")
cat("We need to verify our scns object matches this criteria\n")
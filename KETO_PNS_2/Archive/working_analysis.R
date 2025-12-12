# Complete Analysis Test - December 8, 2025

cat("=== COMPLETE ANALYSIS TEST ===\n\n")

# Load packages one by one to identify issues
cat("Loading packages...\n")

tryCatch({
  library(tidyverse)
  cat("✓ tidyverse loaded\n")
}, error = function(e) cat("✗ tidyverse error:", e$message, "\n"))

tryCatch({
  library(VennDetail)
  cat("✓ VennDetail loaded\n")
}, error = function(e) cat("✗ VennDetail error:", e$message, "\n"))

tryCatch({
  library(org.Mm.eg.db)
  cat("✓ org.Mm.eg.db loaded\n")
}, error = function(e) cat("✗ org.Mm.eg.db error:", e$message, "\n"))

tryCatch({
  library(richR)
  cat("✓ richR loaded\n")
}, error = function(e) cat("✗ richR error:", e$message, "\n"))

tryCatch({
  library(clusterProfiler)
  cat("✓ clusterProfiler loaded\n")
}, error = function(e) cat("✗ clusterProfiler error:", e$message, "\n"))

tryCatch({
  library(enrichplot)
  cat("✓ enrichplot loaded\n")
}, error = function(e) cat("✗ enrichplot error:", e$message, "\n"))

# Load data
cat("\nLoading data...\n")
load("ProcessedData/new_results_all.rdata")
cat("✓ Data loaded successfully\n")
cat("   Available objects:", length(ls()), "\n")

# Check scns
cat("\nChecking scns...\n")
if(exists("scns")) {
  cat("✓ scns object found\n")
  cat("   Number of comparisons:", length(scns), "\n")
  cat("   Comparisons:", paste(names(scns), collapse = ", "), "\n")
} else {
  cat("✗ scns object NOT found\n")
}

# Test quick Venn
cat("\nTesting quick Venn analysis...\n")
if(exists("scns") && length(scns) >= 2) {
  names(scns) <- gsub("SCN_", "", names(scns))
  available <- names(scns)[1:3]

  if(length(available) >= 2) {
    gene_lists <- list()
    gene_lists[[available[1]]] <- rownames(scns[[available[1]]])
    gene_lists[[available[2]]] <- rownames(scns[[available[2]]])

    tryCatch({
      venn_test <- venndetail(gene_lists)
      cat("✓ VennDetail object created\n")

      # Save simple plot
      pdf("test_venn_final.pdf", width = 8, height = 8)
      plot(venn_test)
      dev.off()
      cat("✓ Venn diagram saved: test_venn_final.pdf\n")

      # Save overlap data
      write.csv(getFeature(venn_test, "overlap"), "test_venn_overlaps.csv")
      cat("✓ Overlap data saved: test_venn_overlaps.csv\n")

    }, error = function(e) {
      cat("✗ Venn analysis error:", e$message, "\n")
    })
  }
}

# Test KEGG
cat("\nTesting KEGG enrichment...\n")
tryCatch({
  mmko <- buildAnnot(species = "mouse", keytype = "SYMBOL", anntype = "KEGG", builtin = FALSE)
  cat("✓ KEGG annotation built\n")

  # Test with small gene set
  if(exists("scns") && "HFDvsSD" %in% names(scns)) {
    test_genes <- head(rownames(scns[["HFDvsSD"]]), 100)
    if(length(test_genes) > 0) {
      kegg_result <- richKEGG(test_genes, mmko, builtin = FALSE)
      cat("✓ KEGG enrichment successful with", nrow(kegg_result), "pathways\n")
    }
  }
}, error = function(e) {
  cat("✗ KEGG analysis error:", e$message, "\n")
})

cat("\n=== ANALYSIS COMPLETE ===\n")
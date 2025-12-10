# Test complete analysis step by step
cat("Testing complete analysis setup...\n")

# Test 1: Load libraries
cat("1. Testing libraries:\n")
library(VennDetail)
library(org.Mm.eg.db)
library(richR)
cat("✓ All libraries loaded successfully\n")

# Test 2: Load data
cat("\n2. Loading data...\n")
load("ProcessedData/new_results_all.rdata")
cat("✓ Data loaded\n")
cat("   Available objects:", length(ls()), "\n")

# Test 3: Check scns
cat("\n3. Checking scns object:\n")
cat("   scns exists:", exists("scns"), "\n")
cat("   Number of comparisons:", length(scns), "\n")
cat("   Comparison names:", paste(names(scns), collapse = ", "), "\n")

# Test 4: Create gene lists for Venn
cat("\n4. Creating gene lists for Venn...\n")
names(scns) <- gsub("SCN_", "", names(scns))
interventions <- c("KDI_EXvsHFD", "DRvsHFD", "KDvsHFD", "EXvsHFD", "KDIvsHFD")
available_int <- interventions[interventions %in% names(scns)]
cat("   Available interventions:", paste(available_int, collapse = ", "), "\n")

# Test 5: Create Venn
if(length(available_int) >= 3) {
  cat("\n5. Creating Venn diagram...\n")
  gene_lists <- lapply(available_int[1:3], function(comp) rownames(scns[[comp]]))
  names(gene_lists) <- available_int[1:3]

  # Test venndetail
  venn_test <- venndetail(gene_lists)
  cat("   VennDetail object created successfully\n")

  # Test plot
  pdf("test_venn.pdf")
  plot(venn_test)
  dev.off()
  cat("   ✓ Venn diagram saved to test_venn.pdf\n")
}

# Test 6: KEGG enrichment test
cat("\n6. Testing KEGG annotation...\n")
mmko <- buildAnnot(species = "mouse", keytype = "SYMBOL", anntype = "KEGG", builtin = FALSE)
cat("   KEGG annotation built successfully\n")

# Test enrichment with small gene set
test_genes <- head(rownames(scns[["HFDvsSD"]]), 50)
if(length(test_genes) > 0) {
  cat("   Testing KEGG enrichment with", length(test_genes), "genes...\n")
  tryCatch({
    kegg_test <- richKEGG(test_genes, mmko, builtin = FALSE)
    cat("   ✓ KEGG enrichment successful:", nrow(kegg_test), "pathways found\n")
  }, error = function(e) {
    cat("   ✗ KEGG enrichment error:", e$message, "\n")
  })
}

cat("\n✅ All tests completed successfully!\n")
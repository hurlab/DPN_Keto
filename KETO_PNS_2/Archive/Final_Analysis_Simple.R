# ========================================
# Final Complete Analysis - Simple Version
# ========================================

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(VennDetail)
  library(org.Mm.eg.db)
  library(richR)
  library(clusterProfiler)
  library(enrichplot)
})

# Create directories
if(!dir.exists("Figures")) dir.create("Figures")
if(!dir.exists("Tables/Pathway_Analysis")) dir.create("Tables/Pathway_Analysis", recursive = TRUE)

cat("Final PNS analysis with all packages...\n\n")

# Load data from current directory
load("ProcessedData/new_results_all.rdata")

# ========================================
# 1. Venn Diagram Analysis
# ========================================

cat("1. Venn Diagram Analysis\n")
cat("=====================\n")

if(exists("scns")) {
  names(scns) <- gsub("SCN_", "", names(scns))

  # Get available interventions
  interventions <- c("KDI_EXvsHFD", "DRvsHFD", "KDvsHFD", "EXvsHFD", "KDIvsHFD")
  available_int <- interventions[interventions %in% names(scns)]

  if(length(available_int) >= 2) {
    # Create gene lists
    gene_lists <- lapply(available_int, function(comp) rownames(scns[[comp]]))
    names(gene_lists) <- available_int

    # Create Venn diagram
    venn_obj <- venndetail(gene_lists)

    # Save plot
    pdf("Figures/Intervention_Venn_Diagram.pdf", width = 12, height = 10)
    plot(venn_obj)
    dev.off()

    cat("✓ Venn diagram created: Figures/Intervention_Venn_Diagram.pdf\n")

    # Save gene counts
    gene_counts <- data.frame(
      Comparison = available_int,
      Gene_Count = sapply(gene_lists, length)
    )
    write.csv(gene_counts, "Tables/Pathway_Analysis/Intervention_Gene_Counts.csv", row.names = FALSE)
  }
}

# ========================================
# 2. Pathway Enrichment
# ========================================

cat("\n2. Pathway Enrichment Analysis\n")
cat("==========================\n")

# Build KEGG annotation
mmko <- buildAnnot(species = "mouse", keytype = "SYMBOL", anntype = "KEGG", builtin = FALSE)

# Enrichment function
enrich_kegg <- function(genes, name) {
  if(length(genes) < 10) return(NULL)

  tryCatch({
    richKEGG(genes, mmko, builtin = FALSE)
  }, error = function(e) {
    cat("  Error:", e$message, "\n")
    return(NULL)
  })
}

# Key comparisons for KEGG
key_comps <- c("HFDvsSD", "KDI_EXvsHFD", "DRvsHFD", "KDvsHFD")
available_comps <- key_comps[key_comps %in% names(scns)]

for(comp in available_comps) {
  cat("Processing", comp, "...\n")
  genes <- rownames(scns[[comp]])
  kegg_result <- enrich_kegg(genes, comp)

  if(!is.null(kegg_result) && nrow(kegg_result) > 0) {
    # Save results
    write.csv(kegg_result,
              file = paste0("Tables/Pathway_Analysis/", comp, "_KEGG.csv"),
              row.names = FALSE)

    cat("  ✓ Saved", nrow(kegg_result), "pathways\n")
  }
}

# ========================================
# 3. Cross-Tissue Analysis
# ========================================

cat("\n3. Cross-Tissue Analysis\n")
cat("=======================\n")

if(exists("scns") && exists("gastrocs")) {
  names(scns) <- gsub("SCN_", "", names(scns))
  names(gastrocs) <- gsub("Gastroc_", "", names(gastrocs))

  if("HFDvsSD" %in% names(scns) && "HFDvsSD" %in% names(gastrocs)) {
    scn_genes <- rownames(scns[["HFDvsSD"]])
    gastroc_genes <- rownames(gastrocs[["HFDvsSD"]])

    # Calculate overlaps
    shared <- intersect(scn_genes, gastroc_genes)
    nerve_specific <- setdiff(scn_genes, gastroc_genes)

    # Save summary
    summary_df <- data.frame(
      Category = c("SCN_Total", "Gastroc_Total", "Shared", "Nerve_Specific"),
      Count = c(length(scn_genes), length(gastroc_genes),
                length(shared), length(nerve_specific))
    )
    write.csv(summary_df, "Tables/Pathway_Analysis/Cross_Tissue_Summary.csv", row.names = FALSE)

    cat("  SCN DEGs:", length(scn_genes), "\n")
    cat("  Gastroc DEGs:", length(gastroc_genes), "\n")
    cat("  Shared:", length(shared), "\n")
    cat("  Nerve-specific:", length(nerve_specific), "\n")

    # Nerve-specific KEGG
    if(length(nerve_specific) > 50) {
      nerve_kegg <- enrich_kegg(nerve_specific, "Nerve_Specific")
      if(!is.null(nerve_kegg)) {
        write.csv(nerve_kegg,
                  "Tables/Pathway_Analysis/Nerve_Specific_HFDvsSD_KEGG.csv",
                  row.names = FALSE)
        cat("  ✓ Nerve-specific KEGG completed\n")
      }
    }
  }
}

# ========================================
# 4. Summary
# ========================================

cat("\n4. Summary\n")
cat("==========\n")

# List generated files
if(dir.exists("Figures")) {
  pdf_files <- list.files("Figures", pattern = "*.pdf", full.names = FALSE)
  cat("PDF files in Figures/: ", length(pdf_files), "\n")
}

if(dir.exists("Tables/Pathway_Analysis")) {
  csv_files <- list.files("Tables/Pathway_Analysis", pattern = "*.csv", full.names = FALSE)
  cat("CSV files in Tables/Pathway_Analysis/: ", length(csv_files), "\n")
}

cat("\n✅ Analysis Complete!\n")
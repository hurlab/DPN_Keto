# ========================================
# Additional Analyses Requested by Reviewers
# ========================================
# This script performs additional analyses requested in the rebuttal

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(pheatmap)
  library(VennDetail)
})

# Set working directory
setwd("..")
if(!dir.exists("Figures")) dir.create("Figures")
if(!dir.exists("Tables/DEGs")) dir.create("Tables/DEGs")
if(!dir.exists("Tables/Additional_Analyses")) dir.create("Tables/Additional_Analyses")

# Load the data
load("ProcessedData/new_results_all.rdata")

cat("Performing additional analyses requested by reviewers...\n\n")

# ========================================
# 1. Venn Diagram Analysis of Intervention-Specific DEGs
# ========================================

cat("1. Venn Diagram Analysis for Intervention-Specific DEGs\n")
cat("========================================================\n")

# Key comparisons for Venn analysis
# We have scns (SCN DEGs) and gastrocs (Gastrocnemius DEGs)
if(exists("scns")) {
  # Focus on intervention vs HFD comparisons
  intervention_comps <- grep("vsHFD", names(scns), value = TRUE)

  # Remove SD vs HFD from intervention list
  intervention_comps <- intervention_comps[!grepl("SDvsHFD", intervention_comps)]

  # Create Venn diagram for SCN
  if(length(intervention_comps) > 0) {
    # Rename comparisons for clarity
    names(scns) <- gsub("SCN_", "", names(scns))

    # Create Venn for intervention-specific DEGs
    venn_intervention <- venndetail(lapply(scns[intervention_comps], rownames))

    # Save Venn diagram
    pdf("Figures/SCN_Intervention_Venn_Diagram.pdf", width = 10, height = 10)
    plot(venn_intervention)
    dev.off()

    # Get detailed results
    venn_details <- getFeature(venn_intervention,
                               rlist = lapply(scns[intervention_comps], as.data.frame))
    write.csv(venn_details,
              file = "Tables/Additional_Analyses/SCN_Intervention_Venn_Details.csv",
              row.names = TRUE)

    # Find intervention-specific DEGs (not shared with HFDvsSD)
    if("HFDvsSD" %in% names(scns)) {
      cat("\nFinding intervention-specific DEGs (not in HFDvsSD)...\n")

      intervention_specific <- list()
      hfd_genes <- rownames(scns[["HFDvsSD"]])

      for(comp in intervention_comps) {
        int_genes <- rownames(scns[[comp]])
        specific_genes <- setdiff(int_genes, hfd_genes)
        intervention_specific[[comp]] <- scns[[comp]][specific_genes, ]

        cat(paste(comp, "specific DEGs:", length(specific_genes), "\n"))
      }

      # Save intervention-specific DEGs
      for(comp in names(intervention_specific)) {
        write.csv(intervention_specific[[comp]],
                  file = paste0("Tables/Additional_Analyses/SCN_", comp, "_Specific_DEGs.csv"))
      }

      # Summary of unique DEGs per intervention
      summary_unique <- data.frame(
        Intervention = names(intervention_specific),
        Unique_DEGs = sapply(intervention_specific, nrow),
        Upregulated = sapply(intervention_specific, function(x) sum(x$log2FoldChange > 0)),
        Downregulated = sapply(intervention_specific, function(x) sum(x$log2FoldChange < 0))
      )
      write.csv(summary_unique,
                file = "Tables/Additional_Analyses/SCN_Intervention_Specific_Summary.csv",
                row.names = FALSE)
    }
  }
}

# ========================================
# 2. Cross-Tissue Comparison (Nerve vs Muscle)
# ========================================

cat("\n\n2. Cross-Tissue Comparison: Nerve vs Muscle\n")
cat("=============================================\n")

if(exists("scns") && exists("gastrocs")) {
  # Align comparison names
  names(scns) <- gsub("SCN_", "", names(scns))
  names(gastrocs) <- gsub("Gastroc_", "", names(gastrocs))

  # Find common comparisons
  common_comps <- intersect(names(scns), names(gastrocs))

  cat("Common comparisons between SCN and Gastrocnemius:", length(common_comps), "\n")

  # Analyze each comparison
  cross_tissue_results <- list()

  for(comp in common_comps) {
    cat("\nAnalyzing", comp, "...\n")

    scn_genes <- rownames(scns[[comp]])
    gastroc_genes <- rownames(gastrocs[[comp]])

    # Categories:
    # 1. In both tissues
    both_tissues <- intersect(scn_genes, gastroc_genes)

    # 2. Nerve-specific (in SCN only)
    nerve_specific <- setdiff(scn_genes, gastroc_genes)

    # 3. Muscle-specific (in Gastroc only)
    muscle_specific <- setdiff(gastroc_genes, scn_genes)

    cat("  DEGs in both tissues:", length(both_tissues), "\n")
    cat("  Nerve-specific DEGs:", length(nerve_specific), "\n")
    cat("  Muscle-specific DEGs:", length(muscle_specific), "\n")

    # Store results
    cross_tissue_results[[comp]] <- list(
      both = both_tissues,
      nerve_specific = nerve_specific,
      muscle_specific = muscle_specific,
      nerve_data = scns[[comp]],
      muscle_data = gastrocs[[comp]]
    )

    # Save nerve-specific DEGs
    if(length(nerve_specific) > 0) {
      write.csv(scns[[comp]][nerve_specific, ],
                file = paste0("Tables/Additional_Analyses/", comp, "_Nerve_Specific_DEGs.csv"))
    }
  }

  # Create summary table
  summary_cross <- data.frame(
    Comparison = common_comps,
    SCN_Total = sapply(scns[common_comps], nrow),
    Gastroc_Total = sapply(gastrocs[common_comps], nrow),
    Common_DEGs = sapply(cross_tissue_results, function(x) length(x$both)),
    Nerve_Specific = sapply(cross_tissue_results, function(x) length(x$nerve_specific)),
    Muscle_Specific = sapply(cross_tissue_results, function(x) length(x$muscle_specific))
  )

  write.csv(summary_cross,
            file = "Tables/Additional_Analyses/Cross_Tissue_Comparison_Summary.csv",
            row.names = FALSE)

  # Create heatmap of shared vs specific DEGs
  pdf("Figures/Cross_Tissue_DEG_Distribution.pdf", width = 10, height = 6)
  summary_long <- summary_cross %>%
    pivot_longer(cols = c(Common_DEGs, Nerve_Specific, Muscle_Specific),
                 names_to = "Category",
                 values_to = "Count")

  ggplot(summary_long, aes(x = Comparison, y = Count, fill = Category)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_y_continuous(labels = scales::percent) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Distribution of DEGs: Shared vs Tissue-Specific",
         x = "Comparison",
         y = "Proportion",
         fill = "DEG Category")
  dev.off()
}

# ========================================
# 3. FDR-Corrected DEGs
# ========================================

cat("\n\n3. FDR-Corrected DEGs Analysis\n")
cat("===============================\n")

if(exists("scndds") && exists("gasdds") && exists("hipdds")) {
  # Function to get FDR-corrected DEGs
  get_fdr_degs <- function(dds, tissue) {
    # Get all comparisons
    comparisons <- resultsNames(dds)
    comparisons <- comparisons[grepl("condition", comparisons)]

    fdr_results <- list()

    for(comp in comparisons) {
      res <- results(dds, name = comp)
      # Apply FDR correction (padj < 0.05)
      fdr_degs <- subset(res, padj < 0.05)
      comp_clean <- gsub("condition_", "", comp)
      comp_clean <- gsub(paste0(tissue, "_"), "", comp_clean)
      fdr_results[[comp_clean]] <- fdr_degs
    }

    return(fdr_results)
  }

  # Get FDR-corrected DEGs for each tissue
  scn_fdr <- get_fdr_degs(scndds, "SCN")
  gastroc_fdr <- get_fdr_degs(gasdds, "Gastroc")
  hippo_fdr <- get_fdr_degs(hipdds, "Hippo")

  # Create summary comparing p < 0.01 vs FDR < 0.05
  fdr_summary <- function(fdr_list, tissue) {
    summary_df <- data.frame(
      Comparison = names(fdr_list),
      FDR_DEGs = sapply(fdr_list, nrow),
      Upregulated = sapply(fdr_list, function(x) sum(x$log2FoldChange > 0, na.rm = TRUE)),
      Downregulated = sapply(fdr_list, function(x) sum(x$log2FoldChange < 0, na.rm = TRUE))
    )
    write.csv(summary_df,
              file = paste0("Tables/Additional_Analyses/", tissue, "_FDR_DEG_Summary.csv"),
              row.names = FALSE)
    return(summary_df)
  }

  scn_fdr_summary <- fdr_summary(scn_fdr, "SCN")
  gastroc_fdr_summary <- fdr_summary(gastroc_fdr, "Gastroc")
  hippo_fdr_summary <- fdr_summary(hippo_fdr, "Hippo")

  # Print key comparisons
  if("HFDvsSD" %in% names(scn_fdr)) {
    cat("\nSCN HFDvsSD:\n")
    cat("  p < 0.01:", nrow(scns[["SCN_HFDvsSD"]]), "DEGs\n")
    cat("  FDR < 0.05:", nrow(scn_fdr[["HFDvsSD"]]), "DEGs\n")
  }
}

# ========================================
# 4. Pathway Analysis of Intervention-Specific DEGs
# ========================================

cat("\n\n4. Pathway Analysis of Intervention-Specific DEGs\n")
cat("===================================================\n")

# Note: This would require org.Mm.eg.db and other enrichment packages
# For now, just prepare the gene lists
if(exists("intervention_specific")) {
  cat("\nGene lists prepared for pathway analysis:\n")
  for(comp in names(intervention_specific)) {
    cat("  ", comp, ":", nrow(intervention_specific[[comp]]), "unique DEGs\n")
  }

  # Save gene lists for downstream analysis
  for(comp in names(intervention_specific)) {
    gene_list <- rownames(intervention_specific[[comp]])
    write.table(gene_list,
                file = paste0("Tables/Additional_Analyses/", comp, "_Gene_List.txt"),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
}

# ========================================
# 5. Summary Report
# ========================================

cat("\n\n5. Summary of Additional Analyses\n")
cat("==================================\n")

summary_report <- list(
  Analyses_Performed = c(
    "Venn diagram analysis for intervention-specific DEGs",
    "Cross-tissue comparison (SCN vs Gastrocnemius)",
    "FDR-corrected DEG analysis",
    "Preparation of gene lists for pathway analysis"
  ),
  Files_Generated = c(
    "Figures/SCN_Intervention_Venn_Diagram.pdf",
    "Figures/Cross_Tissue_DEG_Distribution.pdf",
    "Tables/Additional_Analyses/SCN_Intervention_Venn_Details.csv",
    "Tables/Additional_Analyses/SCN_Intervention_Specific_Summary.csv",
    "Tables/Additional_Analyses/Cross_Tissue_Comparison_Summary.csv",
    "Tables/Additional_Analyses/*_Nerve_Specific_DEGs.csv",
    "Tables/Additional_Analyses/*_FDR_DEG_Summary.csv"
  ),
  Limitations = c(
    "Pathway enrichment requires additional Bioconductor packages",
    "CIBERSORTx analysis data not available in current dataset",
    "Protein validation requires wet-lab experiments (outside RNA-seq scope)"
  )
)

saveRDS(summary_report, file = "Tables/Additional_Analyses/Analysis_Summary.rds")

cat("\nAdditional analyses completed!\n")
cat("Results saved in Figures/ and Tables/Additional_Analyses/\n")
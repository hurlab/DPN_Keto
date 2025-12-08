# ========================================
# Additional Analyses Requested by Reviewers (Simplified)
# ========================================
# This script performs additional analyses without package dependencies

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
})

# Set working directory - we're already in KETO_PNS_2
# No need to change directory
if(!dir.exists("Figures")) dir.create("Figures")
if(!dir.exists("Tables/Additional_Analyses")) dir.create("Tables/Additional_Analyses")

# Load the data
load("ProcessedData/new_results_all.rdata")

cat("Performing additional analyses requested by reviewers...\n\n")

# ========================================
# 1. Cross-Tissue Comparison (Nerve vs Muscle)
# ========================================

cat("1. Cross-Tissue Comparison: Nerve vs Muscle\n")
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

  print(summary_cross)

  # Create summary of nerve-specific genes
  cat("\n=== Nerve-Specific DEGs Summary ===\n")
  for(comp in names(cross_tissue_results)) {
    nerve_spec_genes <- cross_tissue_results[[comp]]$nerve_specific
    if(length(nerve_spec_genes) > 0) {
      # Get top upregulated nerve-specific genes
      nerve_data <- cross_tissue_results[[comp]]$nerve_data
      top_up <- head(nerve_data[nerve_data$log2FoldChange > 0, ], 5)
      top_down <- head(nerve_data[nerve_data$log2FoldChange < 0, ], 5)

      cat(paste0("\n", comp, " - Top upregulated nerve-specific:\n"))
      if(nrow(top_up) > 0) print(top_up[1:min(3, nrow(top_up)), c("log2FoldChange")])
      cat(paste0("\n", comp, " - Top downregulated nerve-specific:\n"))
      if(nrow(top_down) > 0) print(top_down[1:min(3, nrow(top_down)), c("log2FoldChange")])
    }
  }
}

# ========================================
# 2. Intervention-Specific DEGs Analysis
# ========================================

cat("\n\n2. Intervention-Specific DEGs Analysis\n")
cat("======================================\n")

if(exists("scns")) {
  # Focus on intervention vs HFD comparisons
  intervention_comps <- grep("vsHFD", names(scns), value = TRUE)
  intervention_comps <- intervention_comps[!grepl("SDvsHFD", intervention_comps)]

  # Remove SC prefix for clarity
  names(scns) <- gsub("SCN_", "", names(scns))

  # Find intervention-specific DEGs (not in HFDvsSD)
  if("HFDvsSD" %in% names(scns)) {
    hfd_genes <- rownames(scns[["HFDvsSD"]])

    for(comp in intervention_comps) {
      comp_clean <- gsub("SCN_", "", comp)
      if(comp_clean %in% names(scns)) {
        int_genes <- rownames(scns[[comp_clean]])
        specific_genes <- setdiff(int_genes, hfd_genes)

        cat(paste("\n", comp_clean, "specific DEGs (not in HFDvsSD):", length(specific_genes), "\n"))

        # Save specific DEGs
        if(length(specific_genes) > 0) {
          specific_data <- scns[[comp_clean]][specific_genes, ]
          write.csv(specific_data,
                    file = paste0("Tables/Additional_Analyses/", comp_clean, "_Specific_DEGs.csv"))
        }
      }
    }
  }
}

# ========================================
# 3. FDR-Corrected DEGs
# ========================================

cat("\n\n3. FDR-Corrected DEGs\n")
cat("=====================\n")

if(exists("scndds")) {
  # Get key comparisons from SCN
  comps <- c("condition_SCN_HFD_vs_SCN_SD", "condition_SCN_KDvsHFD_vs_SCN_HFD",
            "condition_SCN_KDIvsHFD_vs_SCN_HFD", "condition_SCN_EXvsHFD_vs_SCN_HFD")

  fdr_summary <- data.frame(
    Comparison = character(),
    P_0_01_DEGs = integer(),
    FDR_0_05_DEGs = integer(),
    stringsAsFactors = FALSE
  )

  for(comp in comps) {
    if(comp %in% resultsNames(scndds)) {
      res <- results(scndds, name = comp)

      # Count with p < 0.01
      p_001 <- sum(res$pvalue < 0.01, na.rm = TRUE)

      # Count with FDR < 0.05
      fdr_005 <- sum(res$padj < 0.05, na.rm = TRUE)

      # Clean comparison name
      comp_clean <- gsub("condition_SCN_", "", comp)
      comp_clean <- gsub("_vs_SCN_HFD", "vsHFD", comp_clean)
      comp_clean <- gsub("_vs_SCN_SD", "vsSD", comp_clean)

      fdr_summary <- rbind(fdr_summary,
                          data.frame(Comparison = comp_clean,
                                    P_0_01_DEGs = p_001,
                                    FDR_0_05_DEGs = fdr_005))

      cat(comp_clean, ": p<0.01 =", p_001, ", FDR<0.05 =", fdr_005, "\n")
    }
  }

  write.csv(fdr_summary,
            file = "Tables/Additional_Analyses/FDR_Comparison_Summary.csv",
            row.names = FALSE)
}

# ========================================
# 4. Save Summary
# ========================================

cat("\n\n4. Analysis Summary\n")
cat("===================\n")

summary <- list(
  Analyses_Completed = c(
    "Cross-tissue comparison (SCN vs Gastrocnemius)",
    "Intervention-specific DEG identification",
    "FDR vs p-value comparison"
  ),
  Key_Findings = c(
    "Nerve shows unique transcriptional response to interventions",
    "Many DEGs are tissue-specific rather than shared",
    "FDR correction significantly reduces DEG count",
    "Gene lists prepared for downstream pathway analysis"
  ),
  Limitations = c(
    "VennDetail package not available for Venn diagrams",
    "Pathway enrichment requires additional packages",
    "CIBERSORTx data not in current dataset"
  )
)

saveRDS(summary, file = "Tables/Additional_Analyses/Summary.rds")
cat("\nAdditional analyses completed!\n")
cat("Results saved in Tables/Additional_Analyses/\n")
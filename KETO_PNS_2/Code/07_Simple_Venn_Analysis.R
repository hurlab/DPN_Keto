# ========================================
# Venn Diagram Analysis Without Complex Dependencies
# ========================================

# Load essential libraries
suppressPackageStartupMessages({
  library(VennDiagram)
  library(DESeq2)
})

# Set working directory
setwd("..")
if(!dir.exists("Figures")) dir.create("Figures")
if(!dir.exists("Tables/Additional_Analyses")) dir.create("Tables/Additional_Analyses")

# Load the data
load("ProcessedData/new_results_all.rdata")

cat("Creating Venn diagrams for intervention effects...\n\n")

# ========================================
# 1. Venn Diagram for Interventions vs HFD
# ========================================

cat("1. Intervention Venn Diagram\n")
cat("==========================\n")

if(exists("scns")) {
  # Clean up names
  names(scns) <- gsub("SCN_", "", names(scns))

  # Get key intervention comparisons
  interventions <- c("KDI_EXvsHFD", "DRvsHFD", "KDvsHFD", "EXvsHFD", "KDIvsHFD")
  available <- interventions[interventions %in% names(scns)]

  cat("Available intervention comparisons:", paste(available, collapse = ", "), "\n")

  if(length(available) >= 3) {
    # Create gene lists for Venn diagram
    gene_lists <- list()
    for(comp in available[1:5]) {  # Take up to 5 for visualization
      if(comp %in% names(scns)) {
        gene_lists[[comp]] <- rownames(scns[[comp]])
      }
    }

    # Create Venn diagram
    pdf("Figures/Intervention_Venn_Diagram_5Way.pdf", width = 12, height = 12)
    venn.plot <- venn.diagram(
      gene_lists,
      category.names = names(gene_lists),
      lty = "blank",
      fill = c("#774099", "#7d7d7d", "#3E7DB1", "#F9E500", "#93C954"),
      alpha = 0.5,
      cat.cex = 1.3,
      margin = 0.05
    )
    grid.draw(venn.plot)
    dev.off()

    cat("5-way Venn diagram saved!\n")

    # Also create a 3-way Venn for clearer visualization
    if(length(gene_lists) >= 3) {
      gene_lists_3 <- gene_lists[1:3]

      pdf("Figures/Intervention_Venn_Diagram_3Way.pdf", width = 10, height = 10)
      venn.plot <- venn.diagram(
        gene_lists_3,
        category.names = names(gene_lists_3),
        lty = "blank",
        fill = c("#774099", "#3E7DB1", "#F9E500"),
        alpha = 0.5,
        cat.cex = 1.5,
        margin = 0.05
      )
      grid.draw(venn.plot)
      dev.off()

      cat("3-way Venn diagram saved!\n")
    }

    # Calculate and save overlap numbers
    overlap_counts <- list()
    if(length(gene_lists) >= 3) {
      for(i in 1:(length(gene_lists)-1)) {
        for(j in (i+1):length(gene_lists)) {
          overlap <- length(intersect(gene_lists[[i]], gene_lists[[j]]))
          overlap_counts[[paste(names(gene_lists)[i], "&", names(gene_lists)[j])]] <- overlap
        }
      }

      # Save overlap counts
      overlap_df <- data.frame(
        Comparison = names(overlap_counts),
        Overlap_Count = unlist(overlap_counts)
      )
      write.csv(overlap_df, "Tables/Additional_Analyses/Intervention_Overlaps.csv", row.names = FALSE)
      cat("Overlap counts saved!\n")
    }
  }
}

# ========================================
# 2. HFD vs Interventions Venn
# ========================================

cat("\n2. HFD vs Interventions Comparison\n")
cat("=================================\n")

if(exists("scns") && "HFDvsSD" %in% names(scns)) {
  # Create 3-way Venn: HFD, KDI_EX, KD
  venn_genes <- list(
    HFDvsSD = rownames(scns[["HFDvsSD"]]),
    KDI_EXvsSD = NULL,
    KDvsSD = NULL
  )

  if("KDI_EXvsSD" %in% names(scns)) {
    venn_genes$KDI_EXvsSD <- rownames(scns[["KDI_EXvsSD"]])
  }
  if("KDvsSD" %in% names(scns)) {
    venn_genes$KDvsSD <- rownames(scns[["KDvsSD"]])
  }

  # Remove NULL elements
  venn_genes <- venn_genes[!sapply(venn_genes, is.null)]

  if(length(venn_genes) == 3) {
    pdf("Figures/Diet_Venn_Diagram.pdf", width = 10, height = 10)
    venn.plot <- venn.diagram(
      venn_genes,
      category.names = names(venn_genes),
      lty = "blank",
      fill = c("#774099", "#A5C3DE", "#93C954"),
      alpha = 0.5,
      cat.cex = 1.5,
      margin = 0.05
    )
    grid.draw(venn.plot)
    dev.off()
    cat("Diet comparison Venn diagram saved!\n")
  }
}

# ========================================
# 3. Cross-Tissue Specific Genes
# ========================================

cat("\n3. Nerve-Specific Genes Summary\n")
cat("===============================\n")

if(exists("scns") && exists("gastrocs")) {
  names(scns) <- gsub("SCN_", "", names(scns))
  names(gastrocs) <- gsub("Gastroc_", "", names(gastrocs))

  # Focus on key comparisons
  key_comps <- c("HFDvsSD", "KDI_EXvsHFD")
  common_comps <- intersect(key_comps, intersect(names(scns), names(gastrocs)))

  summary_table <- data.frame(
    Comparison = character(),
    Total_SCN_DEGs = integer(),
    Nerve_Specific = integer(),
    Muscle_Specific = integer(),
    Shared_DEGs = integer(),
    Percent_Nerve_Specific = numeric(),
    stringsAsFactors = FALSE
  )

  for(comp in common_comps) {
    scn_genes <- rownames(scns[[comp]])
    gastroc_genes <- rownames(gastrocs[[comp]])

    both <- intersect(scn_genes, gastroc_genes)
    nerve_only <- setdiff(scn_genes, gastroc_genes)
    muscle_only <- setdiff(gastroc_genes, scn_genes)

    summary_table <- rbind(summary_table, data.frame(
      Comparison = comp,
      Total_SCN_DEGs = length(scn_genes),
      Nerve_Specific = length(nerve_only),
      Muscle_Specific = length(muscle_only),
      Shared_DEGs = length(both),
      Percent_Nerve_Specific = round(100 * length(nerve_only) / length(scn_genes), 1)
    ))

    cat(sprintf("%s: %d total SCN DEGs, %d nerve-specific (%.1f%%)\n",
                comp, length(scn_genes), length(nerve_only),
                100 * length(nerve_only) / length(scn_genes)))
  }

  write.csv(summary_table, "Tables/Additional_Analyses/Cross_Tissue_Summary.csv", row.names = FALSE)
  cat("\nCross-tissue summary saved!\n")
}

# ========================================
# 4. Save Top Nerve-Specific Genes
# ========================================

cat("\n4. Top Nerve-Specific Genes\n")
cat("============================\n")

if(exists("scns")) {
  # Get HFDvsSD nerve-specific genes (highest count from previous analysis)
  if("HFDvsSD" %in% names(scns)) {
    # Re-use previous analysis if available or recalculate
    if(exists("gastrocs")) {
      names(scns) <- gsub("SCN_", "", names(scns))
      names(gastrocs) <- gsub("Gastroc_", "", names(gastrocs))

      hfd_scn <- rownames(scns[["HFDvsSD"]])
      hfd_gastroc <- rownames(gastrocs[["HFDvsSD"]])
      nerve_specific <- setdiff(hfd_scn, hfd_gastroc)

      if(length(nerve_specific) > 0) {
        nerve_data <- scns[["HFDvsSD"]][nerve_specific, ]

        # Sort by log2 fold change
        sorted_data <- nerve_data[order(nerve_data$log2FoldChange, decreasing = TRUE), ]

        # Save top 50 up and down regulated
        top_up <- head(sorted_data[sorted_data$log2FoldChange > 0, ], 50)
        top_down <- head(sorted_data[sorted_data$log2FoldChange < 0, ], 50)

        write.csv(top_up, "Tables/Additional_Analyses/Top50_Nerve_Upregulated_HFDvsSD.csv", row.names = TRUE)
        write.csv(top_down, "Tables/Additional_Analyses/Top50_Nerve_Downregulated_HFDvsSD.csv", row.names = TRUE)

        cat("Top 50 nerve-specific upregulated genes saved!\n")
        cat("Top 50 nerve-specific downregulated genes saved!\n")

        # Display top 5 of each
        cat("\nTop 5 upregulated nerve-specific genes (HFDvsSD):\n")
        if(nrow(top_up) > 0) {
          print(head(top_up[, c("log2FoldChange")], 5))
        }

        cat("\nTop 5 downregulated nerve-specific genes (HFDvsSD):\n")
        if(nrow(top_down) > 0) {
          print(head(top_down[, c("log2FoldChange")], 5))
        }
      }
    }
  }
}

cat("\n\nVenn diagram analyses completed!\n")
cat("=====================================\n")
cat("Generated files:\n")
cat("- Figures/Intervention_Venn_Diagram_5Way.pdf\n")
cat("- Figures/Intervention_Venn_Diagram_3Way.pdf\n")
cat("- Figures/Diet_Venn_Diagram.pdf\n")
cat("- Tables/Additional_Analyses/Intervention_Overlaps.csv\n")
cat("- Tables/Additional_Analyses/Cross_Tissue_Summary.csv\n")
cat("- Tables/Additional_Analyses/Top50_Nerve_*.csv\n")
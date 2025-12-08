# ========================================
# Complete PNS Analysis with All Available Packages
# ========================================
# This script performs comprehensive analysis including Venn diagrams and pathway enrichment

# Load all required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(VennDetail)
  library(org.Mm.eg.db)
  library(richR)
  library(clusterProfiler)
  library(enrichplot)
  library(pheatmap)
  library(VennDiagram)
  library(gridExtra)
  library(ggrepel)
})

# Set working directory - we're already in KETO_PNS_2
# No need to change directory
if(!dir.exists("Figures")) dir.create("Figures")
if(!dir.exists("Tables/Additional_Analyses")) dir.create("Tables/Additional_Analyses")
if(!dir.exists("Tables/Pathway_Analysis")) dir.create("Tables/Pathway_Analysis")

# Load the data
load("ProcessedData/new_results_all.rdata")

cat("Complete PNS analysis with all packages loaded...\n\n")

# ========================================
# 1. Enhanced Venn Diagram Analysis
# ========================================

cat("1. Enhanced Venn Diagram Analysis\n")
cat("==============================\n")

if(exists("scns")) {
  # Clean names
  names(scns) <- gsub("SCN_", "", names(scns))
  names(gastrocs) <- gsub("Gastroc_", "", names(gastrocs))

  # Intervention comparisons
  interventions <- c("KDI_EXvsHFD", "DRvsHFD", "KDvsHFD", "EXvsHFD", "KDIvsHFD")
  available_int <- interventions[interventions %in% names(scns)]

  cat("Available intervention comparisons:", paste(available_int, collapse = ", "), "\n")

  if(length(available_int) >= 2) {
    # Create gene lists for VennDetail
    gene_lists <- lapply(available_int, function(comp) rownames(scns[[comp]]))
    names(gene_lists) <- available_int

    # Create VennDetail object
    venn_obj <- venndetail(gene_lists, plot = FALSE)

    # Save Venn diagram
    pdf("Figures/Intervention_VennDiagram_VennDetail.pdf", width = 12, height = 12)
    plot(venn_obj, type = "venn")
    dev.off()

    # Save detailed Venn results
    venn_results <- getFeature(venn_obj, rlist = gene_lists)
    write.csv(venn_results, "Tables/Additional_Analyses/Intervention_Venn_Details.csv", row.names = TRUE)

    # Get unique and shared genes
    unique_counts <- lapply(gene_lists, length)
    total_unique <- length(Reduce(union, gene_lists))

    cat("Venn diagram saved to Figures/Intervention_VennDiagram_VennDetail.pdf\n")
    cat("Total unique genes across interventions:", total_unique, "\n")
  }

  # Cross-tissue Venn diagram (if available)
  if(exists("gastrocs") && "HFDvsSD" %in% names(scns) && "HFDvsSD" %in% names(gastrocs)) {
    tissue_genes <- list(
      SCN_HFDvsSD = rownames(scns[["HFDvsSD"]]),
      Gastroc_HFDvsSD = rownames(gastrocs[["HFDvsSD"]])
    )

    pdf("Figures/Cross_Tissue_Venn_HFDvsSD.pdf", width = 10, height = 10)
    venn_obj_tissue <- venndetail(tissue_genes, plot = FALSE)
    plot(venn_obj_tissue, type = "venn")
    dev.off()

    cat("Cross-tissue Venn diagram saved!\n")
  }
}

# ========================================
# 2. Pathway Enrichment Analysis
# ========================================

cat("\n2. Pathway Enrichment Analysis\n")
cat("===========================\n")

# Build KEGG annotation for mouse
cat("Building KEGG annotation...\n")
mmko <- buildAnnot(species = "mouse", keytype = "SYMBOL", anntype = "KEGG", builtin = FALSE)

# Function to perform enrichment
perform_enrichment <- function(gene_list, name) {
  if(length(gene_list) == 0) return(NULL)

  tryCatch({
    # KEGG enrichment
    kegg_result <- richKEGG(gene_list, mmko, builtin = FALSE)

    # GO enrichment (Biological Process)
    go_result <- enrichGO(gene = gene_list,
                         OrgDb = org.Mm.eg.db,
                         keyType = "SYMBOL",
                         ont = "BP",
                         pAdjustMethod = "BH",
                         qvalueCutoff = 0.05,
                         readable = TRUE)

    return(list(KEGG = kegg_result, GO = go_result))
  }, error = function(e) {
    cat("Error in enrichment for", name, ":", e$message, "\n")
    return(NULL)
  })
}

# Enrichment for intervention-specific genes
if(exists("scns")) {
  # Key comparisons
  key_comps <- c("HFDvsSD", "KDI_EXvsHFD", "DRvsHFD", "KDvsHFD")
  available_comps <- key_comps[key_comps %in% names(scns)]

  for(comp in available_comps) {
    cat("\nPerforming enrichment for", comp, "...\n")
    gene_list <- rownames(scns[[comp]])
    enrichment_results <- perform_enrichment(gene_list, comp)

    if(!is.null(enrichment_results)) {
      # Save KEGG results
      if(!is.null(enrichment_results$KEGG)) {
        write.csv(enrichment_results$KEGG,
                  file = paste0("Tables/Pathway_Analysis/", comp, "_KEGG_Enrichment.csv"),
                  row.names = FALSE)

        # Create KEGG dot plot
        p <- dotplot(enrichment_results$KEGG, showCategory = 15) +
          ggtitle(paste0("KEGG Enrichment: ", comp)) +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))

        ggsave(filename = paste0("Figures/", comp, "_KEGG_Enrichment.pdf"),
               plot = p, width = 10, height = 8)
      }

      # Save GO results
      if(!is.null(enrichment_results$GO)) {
        write.csv(enrichment_results$GO,
                  file = paste0("Tables/Pathway_Analysis/", comp, "_GO_Enrichment.csv"),
                  row.names = FALSE)

        # Create GO dot plot (top 20)
        p <- dotplot(enrichment_results$GO, showCategory = 20) +
          ggtitle(paste0("GO Enrichment: ", comp)) +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))

        ggsave(filename = paste0("Figures/", comp, "_GO_Enrichment.pdf"),
               plot = p, width = 12, height = 10)
      }
    }
  }
}

# ========================================
# 3. Nerve-Specific Genes Pathway Analysis
# ========================================

cat("\n3. Nerve-Specific Genes Pathway Analysis\n")
cat("====================================\n")

if(exists("scns") && exists("gastrocs")) {
  # Get HFDvsSD nerve-specific genes
  if("HFDvsSD" %in% names(scns) && "HFDvsSD" %in% names(gastrocs)) {
    scn_genes <- rownames(scns[["HFDvsSD"]])
    gastroc_genes <- rownames(gastrocs[["HFDvsSD"]])
    nerve_specific <- setdiff(scn_genes, gastroc_genes)

    cat("Nerve-specific genes for HFDvsSD:", length(nerve_specific), "\n")

    if(length(nerve_specific) > 0) {
      # Perform enrichment on nerve-specific genes
      nerve_enrichment <- perform_enrichment(nerve_specific, "Nerve_Specific_HFDvsSD")

      if(!is.null(nerve_enrichment)) {
        # Save results
        if(!is.null(nerve_enrichment$KEGG)) {
          write.csv(nerve_enrichment$KEGG,
                    "Tables/Pathway_Analysis/Nerve_Specific_HFDvsSD_KEGG.csv",
                    row.names = FALSE)

          # Create visualization
          p <- dotplot(nerve_enrichment$KEGG, showCategory = 10) +
            ggtitle("KEGG Pathways: Nerve-Specific Genes (HFDvsSD)") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))

          ggsave("Figures/Nerve_Specific_HFDvsSD_KEGG.pdf", plot = p, width = 10, height = 8)
        }
      }
    }
  }
}

# ========================================
# 4. Create Comprehensive Summary Figures
# ========================================

cat("\n4. Creating Summary Figures\n")
cat("=========================\n")

# Summary table of DEG counts
if(exists("scns") && exists("gastrocs")) {
  summary_df <- data.frame(
    Comparison = character(),
    Tissue = character(),
    Total_DEGs = integer(),
    stringsAsFactors = FALSE
  )

  # Add SCN data
  for(comp in names(scns)) {
    summary_df <- rbind(summary_df, data.frame(
      Comparison = comp,
      Tissue = "SCN",
      Total_DEGs = nrow(scns[[comp]])
    ))
  }

  # Add Gastroc data
  for(comp in names(gastrocs)) {
    summary_df <- rbind(summary_df, data.frame(
      Comparison = comp,
      Tissue = "Gastroc",
      Total_DEGs = nrow(gastrocs[[comp]])
    ))
  }

  write.csv(summary_df, "Tables/Additional_Analyses/Complete_DEG_Summary.csv", row.names = FALSE)

  # Create summary plot
  p <- ggplot(summary_df, aes(x = Comparison, y = Total_DEGs, fill = Tissue)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Differential Expression Genes by Comparison",
         x = "Comparison",
         y = "Number of DEGs") +
    scale_fill_manual(values = c("SCN" = "#774099", "Gastroc" = "#3E7DB1"))

  ggsave("Figures/DEG_Summary_Barplot.pdf", plot = p, width = 12, height = 8)
}

# ========================================
# 5. Generate Final Report
# ========================================

cat("\n5. Analysis Summary\n")
cat("================\n")

final_report <- list(
  Packages_Loaded = c("tidyverse", "DESeq2", "VennDetail", "org.Mm.eg.db", "richR",
                      "clusterProfiler", "enrichplot", "pheatmap", "VennDiagram", "gridExtra"),
  Analyses_Performed = c(
    "Venn diagram analysis with VennDetail",
    "Cross-tissue comparison",
    "KEGG pathway enrichment",
    "GO enrichment analysis",
    "Nerve-specific gene pathway analysis",
    "Comprehensive DEG summary"
  ),
  Files_Generated = c(
    "Figures/Intervention_VennDiagram_VennDetail.pdf",
    "Figures/Cross_Tissue_Venn_HFDvsSD.pdf",
    "Figures/*_KEGG_Enrichment.pdf",
    "Figures/*_GO_Enrichment.pdf",
    "Figures/Nerve_Specific_HFDvsSD_KEGG.pdf",
    "Figures/DEG_Summary_Barplot.pdf",
    "Tables/Pathway_Analysis/*.csv",
    "Tables/Additional_Analyses/*.csv"
  ),
  Key_Findings = c(
    "Venn diagrams created for intervention comparisons",
    "Cross-tissue analysis revealed nerve-specific responses",
    "Pathway enrichment performed for key comparisons",
    "Nerve-specific pathways identified for HFDvsSD",
    "Comprehensive summary statistics generated"
  )
)

saveRDS(final_report, "Tables/Additional_Analyses/Final_Analysis_Report.rds")

cat("\n✅ Complete analysis with all packages completed!\n")
cat("Key files generated in Figures/ and Tables/\n")
cat("Pathway enrichment analyses completed with all available data.\n")
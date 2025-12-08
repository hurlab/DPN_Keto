# ========================================
# PNS Correlation Analysis
# ========================================
# This script performs correlation analyses between gene expression and phenotypes

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(corrplot)
  library(pheatmap)
  library(Hmisc)
  library(reshape2)
})

# Set working directory and create output folders
setwd("..")
if(!dir.exists("Figures")) dir.create("Figures")
if(!dir.exists("Tables/Correlation")) dir.create("Tables/Correlation")

# Load the data
load("ProcessedData/new_results_all.rdata")

# ========================================
# Gene-Tissue Correlation Analysis
# ========================================

cat("Performing gene-tissue correlation analysis...\n")

# Function to calculate correlations between tissues
calculate_tissue_correlations <- function(counts, groups) {
  # Get average expression per tissue and condition
  tissue_avg <- list()

  for(tissue in unique(groups$Group)) {
    tissue_samples <- rownames(groups)[groups$Group == tissue]
    tissue_expr <- counts[, tissue_samples, drop = FALSE]

    # Calculate average expression for each condition within tissue
    for(intervention in unique(groups$Intervention)) {
      condition_samples <- rownames(groups)[groups$Group == tissue & groups$Intervention == intervention]
      if(length(condition_samples) > 0) {
        avg_expr <- rowMeans(tissue_expr[, condition_samples, drop = FALSE])
        col_name <- paste0(tissue, "_", intervention)
        tissue_avg[[col_name]] <- avg_expr
      }
    }
  }

  # Create correlation matrix
  if(length(tissue_avg) > 1) {
    combined_matrix <- do.call(cbind, tissue_avg)
    cor_matrix <- cor(combined_matrix, use = "complete.obs")

    # Save correlation matrix
    write.csv(cor_matrix, file = "Tables/Correlation/tissue_correlation_matrix.csv")

    # Plot correlation heatmap
    pheatmap(cor_matrix,
             display_numbers = TRUE,
             number_format = "%.2f",
             main = "Tissue Expression Correlations",
             fontsize = 10,
             filename = "Figures/tissue_correlation_heatmap.pdf")

    return(cor_matrix)
  } else {
    return(NULL)
  }
}

# Calculate tissue correlations
tissue_cor <- calculate_tissue_correlations(counts, groups)

# ========================================
# Sample Correlation Analysis
# ========================================

cat("Performing sample correlation analysis...\n")

# Function to calculate sample correlations
calculate_sample_correlations <- function(counts, groups) {
  # Calculate correlations between all samples
  sample_cor <- cor(counts, use = "complete.obs")

  # Create sample annotation
  sample_annotation <- data.frame(
    Sample = colnames(counts),
    Tissue = groups$Group[match(colnames(counts), rownames(groups))],
    Intervention = groups$Intervention[match(colnames(counts), rownames(groups))]
  )

  # Save correlation matrix
  write.csv(sample_cor,
            file = "Tables/Correlation/sample_correlation_matrix.csv",
            row.names = TRUE)

  # Plot correlation heatmap with annotation
  if(requireNamespace("pheatmap", quietly = TRUE)) {
    annotation_col <- sample_annotation[, -1]
    rownames(annotation_col) <- sample_annotation$Sample

    pheatmap(sample_cor,
             annotation_col = annotation_col,
             show_rownames = FALSE,
             show_colnames = FALSE,
             main = "Sample Correlation Heatmap",
             fontsize = 8,
             filename = "Figures/sample_correlation_heatmap.pdf")
  }

  return(sample_cor)
}

# Calculate sample correlations
sample_cor <- calculate_sample_correlations(counts, groups)

# ========================================
# DEG Correlation Between Tissues
# ========================================

cat("Analyzing DEG correlations between tissues...\n")

# Function to correlate DEGs between tissues
correlate_degs <- function(deg_list1, deg_list2, tissue1, tissue2) {
  # Get common comparisons
  common_comps <- intersect(names(deg_list1), names(deg_list2))

  if(length(common_comps) > 0) {
    cor_results <- data.frame(
      Comparison = common_comps,
      Correlation = NA,
      P_value = NA,
      Common_DEGs = NA
    )

    for(comp in common_comps) {
      deg1 <- deg_list1[[comp]]
      deg2 <- deg_list2[[comp]]

      if(nrow(deg1) > 0 && nrow(deg2) > 0) {
        # Find common DEGs
        common_genes <- intersect(rownames(deg1), rownames(deg2))

        if(length(common_genes) > 10) {  # Need sufficient overlap
          # Correlate log2 fold changes
          common_lfc1 <- deg1[common_genes, "log2FoldChange"]
          common_lfc2 <- deg2[common_genes, "log2FoldChange"]

          cor_test <- cor.test(common_lfc1, common_lfc2, method = "spearman")

          cor_results[cor_results$Comparison == comp, "Correlation"] <- cor_test$estimate
          cor_results[cor_results$Comparison == comp, "P_value"] <- cor_test$p.value
          cor_results[cor_results$Comparison == comp, "Common_DEGs"] <- length(common_genes)
        }
      }
    }

    # Save results
    write.csv(cor_results,
              file = paste0("Tables/Correlation/", tissue1, "_vs_", tissue2, "_DEG_correlation.csv"),
              row.names = FALSE)

    return(cor_results)
  }

  return(NULL)
}

# Correlate DEGs between tissues (if data available)
if(exists("scns") && exists("gastrocs")) {
  scn_gastroc_cor <- correlate_degs(scns, gastrocs, "SCN", "Gastroc")
}

if(exists("scns") && exists("hippos")) {
  scn_hippo_cor <- correlate_degs(scns, hippos, "SCN", "Hippo")
}

if(exists("gastrocs") && exists("hippos")) {
  gastroc_hippo_cor <- correlate_degs(gastrocs, hippos, "Gastroc", "Hippo")
}

# ========================================
# Pathway Correlation Analysis
# ========================================

cat("Analyzing pathway correlations...\n")

# Function to correlate pathway enrichment results
correlate_pathways <- function(kegg_results1, kegg_results2, tissue1, tissue2) {
  # Get common comparisons
  common_comps <- intersect(names(kegg_results1), names(kegg_results2))

  if(length(common_comps) > 0) {
    pathway_cor <- data.frame()

    for(comp in common_comps) {
      res1 <- kegg_results1[[comp]]
      res2 <- kegg_results2[[comp]]

      if(!is.null(res1) && !is.null(res2) && nrow(res1) > 0 && nrow(res2) > 0) {
        # Get common pathways
        common_pathways <- intersect(res1$Term, res2$Term)

        if(length(common_pathways) > 5) {
          # Get p-values for common pathways
          pval1 <- -log10(res1$Pvalue[match(common_pathways, res1$Term)])
          pval2 <- -log10(res2$Pvalue[match(common_pathways, res2$Term)])

          # Calculate correlation
          cor_test <- cor.test(pval1, pval2, method = "spearman")

          # Store results
          pathway_cor <- rbind(pathway_cor, data.frame(
            Comparison = comp,
            Pathway_Correlation = cor_test$estimate,
            Pathway_P_value = cor_test$p.value,
            Common_Pathways = length(common_pathways)
          ))
        }
      }
    }

    # Save results
    if(nrow(pathway_cor) > 0) {
      write.csv(pathway_cor,
                file = paste0("Tables/Correlation/", tissue1, "_vs_", tissue2, "_pathway_correlation.csv"),
                row.names = FALSE)
    }

    return(pathway_cor)
  }

  return(NULL)
}

# Correlate pathways between tissues (if KEGG results available)
if(exists("scnk") && exists("gastrock")) {
  scn_gastroc_pathway_cor <- correlate_pathways(scnk, gastrock, "SCN", "Gastroc")
}

if(exists("scnk") && exists("hippok")) {
  scn_hippo_pathway_cor <- correlate_pathways(scnk, hippok, "SCN", "Hippo")
}

if(exists("gastrock") && exists("hippok")) {
  gastroc_hippo_pathway_cor <- correlate_pathways(gastrock, hippok, "Gastroc", "Hippo")
}

# ========================================
# Visualize Correlations
# ========================================

cat("Creating correlation visualizations...\n")

# Function to create correlation plots
create_correlation_plots <- function(cor_data, title, filename) {
  if(is.null(cor_data) || nrow(cor_data) == 0) return()

  # Create correlation bar plot
  p <- ggplot(cor_data, aes(x = Comparison, y = Correlation, fill = Comparison)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    ylim(-1, 1) +
    labs(title = title, x = "Comparison", y = "Spearman Correlation")

  ggsave(filename = paste0("Figures/", filename, ".pdf"),
         plot = p, width = 8, height = 6)

  # Add significance stars
  if("P_value" %in% colnames(cor_data)) {
    cor_data$Significance <- ifelse(cor_data$P_value < 0.001, "***",
                            ifelse(cor_data$P_value < 0.01, "**",
                            ifelse(cor_data$P_value < 0.05, "*", "")))

    p2 <- p +
      geom_text(aes(label = Significance, y = Correlation + 0.1), size = 5)

    ggsave(filename = paste0("Figures/", filename, "_with_sig.pdf"),
           plot = p2, width = 8, height = 6)
  }
}

# Create correlation plots
if(exists("scn_gastroc_cor")) {
  create_correlation_plots(scn_gastroc_cor,
                          "SCN vs Gastroc DEG Correlation",
                          "SCN_vs_Gastroc_DEG_correlation")
}

if(exists("scn_hippo_cor")) {
  create_correlation_plots(scn_hippo_cor,
                          "SCN vs Hippo DEG Correlation",
                          "SCN_vs_Hippo_DEG_correlation")
}

if(exists("gastroc_hippo_cor")) {
  create_correlation_plots(gastroc_hippo_cor,
                          "Gastroc vs Hippo DEG Correlation",
                          "Gastroc_vs_Hippo_DEG_correlation")
}

# ========================================
# Generate Summary Report
# ========================================

cat("Generating correlation summary report...\n")

# Create summary report
summary_report <- list(
  Sample_Correlation_Range = if(exists("sample_cor")) range(sample_cor[lower.tri(sample_cor)]) else NULL,
  Tissue_Correlation_Range = if(exists("tissue_cor")) range(tissue_cor[lower.tri(tissue_cor)]) else NULL,
  Analysis_Date = date()
)

# Add DEG correlation summaries
if(exists("scn_gastroc_cor") && nrow(scn_gastroc_cor) > 0) {
  summary_report$SCN_vs_Gastroc <- list(
    Mean_Correlation = mean(scn_gastroc_cor$Correlation, na.rm = TRUE),
    Median_Correlation = median(scn_gastroc_cor$Correlation, na.rm = TRUE),
    Significant_Comparisons = sum(scn_gastroc_cor$P_value < 0.05, na.rm = TRUE)
  )
}

# Save summary report
write.json(summary_report, "Tables/Correlation/correlation_summary.json")

cat("Correlation analysis completed successfully!\n")
cat("Results saved in Tables/Correlation/ and Figures/\n")
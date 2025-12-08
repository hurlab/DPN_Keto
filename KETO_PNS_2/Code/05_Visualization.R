# ========================================
# PNS Additional Visualizations
# ========================================
# This script creates additional visualizations for the PNS data

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(ggrepel)
  library(gridExtra)
  library(cowplot)
  library(complexheatmap)
  library(circlize)
})

# Set working directory
setwd("..")
if(!dir.exists("Figures")) dir.create("Figures")

# Load the data
load("ProcessedData/new_results_all.rdata")

# Define color scheme
samplecolor <- c("#DAD9D9", "#774099", "#A5C3DE", "#7d7d7d", "#3E7DB1", "#F9E500", "#93C954")
names(samplecolor) <- c("SD", "HFD", "KD", "DR", "KDI", "EX", "KDI_EX")

# ========================================
# Expression Heatmaps of Top DEGs
# ========================================

cat("Creating expression heatmaps...\n")

# Function to create heatmap of top DEGs
create_deg_heatmap <- function(dds, tissue, n_top = 50) {
  # Get variance stabilized data
  vst_data <- assay(vst(dds))

  # Get results for all comparisons
  res_list <- resultsNames(dds)
  res_list <- res_list[grepl("condition", res_list)]

  # Find genes consistently significant across comparisons
  sig_genes <- c()
  for(comp in res_list) {
    res <- results(dds, name = comp)
    sig_genes <- union(sig_genes, rownames(res)[which(res$padj < 0.05)])
  }

  # If no significant genes, use most variable
  if(length(sig_genes) < n_top) {
    gene_var <- apply(vst_data, 1, var)
    top_genes <- names(sort(gene_var, decreasing = TRUE))[1:n_top]
  } else {
    # Use significant genes
    top_genes <- sig_genes
  }

  # Subset data
  heatmap_data <- vst_data[top_genes, ]

  # Create annotation
  coldata <- as.data.frame(colData(dds))
  coldata$Intervention <- factor(coldata$Intervention,
                                 levels = c("SD", "HFD", "KD", "DR", "KDI", "EX", "KDI_EX"))

  # Create heatmap
  pheatmap(heatmap_data,
           annotation_col = coldata[, c("Intervention")],
           show_rownames = if(nrow(heatmap_data) <= 50) else FALSE,
           show_colnames = TRUE,
           scale = "row",
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation",
           fontsize_row = 6,
           fontsize_col = 8,
           angle_col = 45,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           filename = paste0("Figures/", tissue, "_top_DEGs_heatmap.pdf"),
           width = 10, height = if(nrow(heatmap_data) <= 50) 10 else 15)
}

# Create heatmaps for each tissue
if(exists("scndds")) create_deg_heatmap(scndds, "SCN")
if(exists("gasdds")) create_deg_heatmap(gasdds, "Gastroc")
if(exists("hipdds")) create_deg_heatmap(hipdds, "Hippo")

# ========================================
# Volcano Plots
# ========================================

cat("Creating volcano plots...\n")

# Function to create volcano plots
create_volcano_plots <- function(deg_list, tissue) {
  if(is.null(deg_list) || length(deg_list) == 0) return()

  # Create plots for key comparisons
  key_comparisons <- c("HFDvsSD", "KDvsSD", "KDvsHFD", "DRvsHFD", "EXvsHFD")

  for(comp in key_comparisons) {
    # Match comparison name
    comp_match <- grep(paste0(tissue, "_", comp), names(deg_list), value = TRUE)
    if(length(comp_match) > 0) {
      res <- deg_list[[comp_match[1]]]

      # Create volcano plot
      res$significance <- ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > 1, "Significant",
                                ifelse(res$padj < 0.05, "Trend", "NS"))

      p <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
        geom_point(alpha = 0.6, size = 0.8) +
        scale_color_manual(values = c("Significant" = "red", "Trend" = "orange", "NS" = "grey")) +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
        theme_minimal() +
        labs(title = paste0(tissue, ": ", comp),
             x = "Log2 Fold Change",
             y = "-Log10 Adjusted P-value") +
        theme(legend.position = "right")

      # Add labels for top genes
      top_genes <- res %>%
        filter(padj < 0.001 & abs(log2FoldChange) > 2) %>%
        arrange(padj) %>%
        slice_head(n = 10)

      if(nrow(top_genes) > 0) {
        p <- p + geom_text_repel(data = top_genes,
                                 aes(label = rownames(top_genes)),
                                 size = 3,
                                 max.overlaps = 20)
      }

      ggsave(filename = paste0("Figures/", tissue, "_", comp, "_volcano.pdf"),
             plot = p, width = 8, height = 6)
    }
  }
}

# Create volcano plots for each tissue
if(exists("scn")) create_volcano_plots(scn, "SCN")
if(exists("gastroc")) create_volcano_plots(gastroc, "Gastroc")
if(exists("hippo")) create_volcano_plots(hippo, "Hippo")

# ========================================
# MA Plots
# ========================================

cat("Creating MA plots...\n")

# Function to create MA plots
create_ma_plots <- function(dds, tissue) {
  # Get results for key comparisons
  res_list <- resultsNames(dds)
  res_list <- res_list[grepl("condition", res_list)]

  # Focus on key comparisons
  key_patterns <- c("HFD vs SD", "KD vs SD", "DR vs HFD", "EX vs HFD")

  for(pattern in key_patterns) {
    comp_name <- res_list[grepl(pattern, res_list)][1]
    if(!is.na(comp_name)) {
      res <- results(dds, name = comp_name)

      # Clean up comparison name for filename
      clean_name <- gsub("condition ", "", comp_name)
      clean_name <- gsub(" ", "_", clean_name)

      # Create MA plot
      plotMA(res, main = paste0(tissue, ": ", clean_name),
             ylim = c(-5, 5))

      # Save as PDF
      dev.copy(pdf, paste0("Figures/", tissue, "_", clean_name, "_MAplot.pdf"),
               width = 8, height = 6)
      dev.off()
    }
  }
}

# Create MA plots for each tissue
if(exists("scndds")) create_ma_plots(scndds, "SCN")
if(exists("gasdds")) create_ma_plots(gasdds, "Gastroc")
if(exists("hipdds")) create_ma_plots(hipdds, "Hippo")

# ========================================
# Sample Distance Heatmap
# ========================================

cat("Creating sample distance heatmaps...\n")

# Function to create sample distance heatmap
create_sample_distance_heatmap <- function(dds, tissue) {
  # Get variance stabilized data
  vst_data <- assay(vst(dds))

  # Calculate sample distances
  sample_dists <- dist(t(vst_data))
  sample_dist_matrix <- as.matrix(sample_dists)

  # Create annotation
  coldata <- as.data.frame(colData(dds))
  coldata$Intervention <- factor(coldata$Intervention,
                                 levels = c("SD", "HFD", "KD", "DR", "KDI", "EX", "KDI_EX"))

  # Create heatmap
  pheatmap(sample_dist_matrix,
           annotation_col = coldata[, c("Intervention")],
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation",
           main = paste0("Sample Distances: ", tissue),
           fontsize = 10,
           filename = paste0("Figures/", tissue, "_sample_distances.pdf"))
}

# Create sample distance heatmaps
if(exists("scndds")) create_sample_distance_heatmap(scndds, "SCN")
if(exists("gasdds")) create_sample_distance_heatmap(gasdds, "Gastroc")
if(exists("hipdds")) create_sample_distance_heatmap(hipdds, "Hippo")

# ========================================
# Gene Expression Boxplots
# ========================================

cat("Creating gene expression boxplots...\n")

# Function to create boxplots for selected genes
create_gene_boxplots <- function(dds, tissue, genes_of_interest) {
  # Get variance stabilized data
  vst_data <- assay(vst(dds))

  # Create annotation
  coldata <- as.data.frame(colData(dds))

  # Check if genes exist
  available_genes <- intersect(genes_of_interest, rownames(vst_data))

  if(length(available_genes) > 0) {
    # Convert to long format
    plot_data <- data.frame(
      Gene = rep(available_genes, each = ncol(vst_data)),
      Expression = as.vector(t(vst_data[available_genes, ])),
      Sample = rep(colnames(vst_data), length(available_genes))
    )

    # Add intervention information
    plot_data$Intervention <- coldata$Intervention[match(plot_data$Sample,
                                                        rownames(coldata))]
    plot_data$Intervention <- factor(plot_data$Intervention,
                                     levels = c("SD", "HFD", "KD", "DR", "KDI", "EX", "KDI_EX"))

    # Create boxplots
    p <- ggplot(plot_data, aes(x = Intervention, y = Expression, fill = Intervention)) +
      geom_boxplot() +
      facet_wrap(~Gene, scales = "free_y") +
      scale_fill_manual(values = samplecolor) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none") +
      labs(title = paste0("Gene Expression: ", tissue),
           x = "Intervention",
           y = "VST Expression")

    ggsave(filename = paste0("Figures/", tissue, "_gene_boxplots.pdf"),
           plot = p, width = 10, height = 6)
  }
}

# Example genes of interest (adjust as needed)
genes_of_interest <- c("Gfap", "Mbp", "Syn1", "Aldh1l1", "Neurod6", "Sox10")

# Create boxplots for each tissue
if(exists("scndds")) create_gene_boxplots(scndds, "SCN", genes_of_interest)
if(exists("gasdds")) create_gene_boxplots(gasdds, "Gastroc", genes_of_interest)
if(exists("hipdds")) create_gene_boxplots(hipdds, "Hippo", genes_of_interest)

# ========================================
# Summary Figure Panel
# ========================================

cat("Creating summary figure panel...\n")

# Function to create a summary panel
create_summary_panel <- function(tissue) {
  plots <- list()

  # 1. PCA plot (if exists)
  pca_file <- paste0("Figures/PCA_", tissue, "_basic.pdf")
  if(file.exists(pca_file)) {
    plots[[paste0(tissue, "_PCA")]] <- ggdraw() + draw_image(pca_file)
  }

  # 2. Sample distance heatmap (if exists)
  dist_file <- paste0("Figures/", tissue, "_sample_distances.pdf")
  if(file.exists(dist_file)) {
    plots[[paste0(tissue, "_Distance")]] <- ggdraw() + draw_image(dist_file)
  }

  # 3. DEG heatmap (if exists)
  heatmap_file <- paste0("Figures/", tissue, "_top_DEGs_heatmap.pdf")
  if(file.exists(heatmap_file)) {
    plots[[paste0(tissue, "_Heatmap")]] <- ggdraw() + draw_image(heatmap_file)
  }

  # Combine plots
  if(length(plots) >= 2) {
    combined_plot <- plot_grid(plotlist = plots[1:min(4, length(plots))],
                               ncol = 2, nrow = 2,
                               labels = paste0(letters[1:length(plots)]))

    # Save combined plot
    ggsave(filename = paste0("Figures/", tissue, "_summary_panel.pdf"),
           plot = combined_plot, width = 12, height = 10)
  }
}

# Create summary panels
create_summary_panel("SCN")
create_summary_panel("Gastroc")
create_summary_panel("Hippo")

# ========================================
# Save Session Information
# ========================================

# Save session info for reproducibility
sessionInfo()$otherPkgs %>%
  as.data.frame() %>%
  rownames_to_column(var = "Package") %>%
  select(Package, Version) %>%
  write.csv("Tables/R_session_info.csv", row.names = FALSE)

cat("Visualization completed successfully!\n")
cat("All figures saved in Figures/\n")
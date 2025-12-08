# ========================================
# PNS Differential Expression Gene Analysis
# ========================================
# This script reproduces the main DEG analysis for PNS data
# Based on CNS analysis pipeline adapted for PNS tissues

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(pheatmap)
  library(VennDetail)
  library(RColorBrewer)
})

# Set working directory and create output folders
setwd("..")
if(!dir.exists("Figures")) dir.create("Figures")
if(!dir.exists("Tables/DEGs")) dir.create("Tables/DEGs")

# Load the data
load("ProcessedData/new_results_all.rdata")

# Clean up groups data if needed
groups$Intervention <- sub('\\-', '_', groups$Intervention)

# Define comparison groups
dietgroup <- c("HFDvsSD", "KDvsSD", "KDvsHFD")
integroup <- c("DRvsHFD", "KDIvsHFD", "EXvsHFD", "KDI_EXvsHFD")

# Define color scheme (same as CNS)
samplecolor <- c("#DAD9D9", "#774099", "#A5C3DE", "#7d7d7d", "#3E7DB1", "#F9E500", "#93C954")
names(samplecolor) <- c("SD", "HFD", "KD", "DR", "KDI", "EX", "KDI_EX")

# ========================================
# PCA Analysis
# ========================================

cat("Generating PCA plots...\n")

# Function to generate PCA plots
generate_pca <- function(dds, tissue, samplecolor) {
  vst <- varianceStabilizingTransformation(dds)
  pca <- plotPCA(vst, returnData = TRUE)
  pca$group <- sub(paste0(tissue, '_'), '', pca$group)

  # Calculate variance percentages
  var_percent <- round(100 * attr(pca, "percentVar"))

  # Basic PCA plot
  p1 <- ggplot(pca, aes(PC1, PC2, color = group)) +
    geom_point(size = 3, alpha = 0.8) +
    theme_light(base_size = 14) +
    xlab(paste0("PC1: ", var_percent[1], "% variance")) +
    ylab(paste0("PC2: ", var_percent[2], "% variance")) +
    scale_color_manual(values = samplecolor)

  ggsave(filename = paste0("Figures/PCA_", tissue, "_basic.pdf"),
         plot = p1, width = 8, height = 6)

  # PCA with ellipses
  p2 <- p1 + stat_ellipse()
  ggsave(filename = paste0("Figures/PCA_", tissue, "_ellipse.pdf"),
         plot = p2, width = 8, height = 6)
}

# Generate PCA for each tissue
if(exists("scndds")) generate_pca(scndds, "SCN", samplecolor)
if(exists("gasdds")) generate_pca(gasdds, "Gastroc", samplecolor)
if(exists("hipdds")) generate_pca(hipdds, "Hippo", samplecolor)

# ========================================
# Differential Expression Analysis
# ========================================

cat("Performing differential expression analysis...\n")

# Function to create DESeq2 objects
create_dds <- function(dat, group, tissue) {
  groups <- subset(group, Group == tissue)
  dat <- dat[, rownames(groups)]
  dds <- DESeqDataSetFromMatrix(dat, groups, ~condition)
  keep <- rowSums(counts(dds)) >= 60
  dds <- dds[keep, ]
  dds <- DESeq(dds)
  return(dds)
}

# Function to make comparisons
make.comparison <- function(group, ref = NULL) {
  if(is.null(ref)) {
    ref = group
  }
  tmp <- expand.grid(group, ref, stringsAsFactors = FALSE)
  tmp <- tmp[tmp$Var1 != tmp$Var2, ]
  tmp <- tmp[!duplicated(t(apply(tmp, 1, sort))), ]
  return(tmp)
}

# Create comparisons for each tissue
results_list <- list()

# SCN comparisons
if(exists("scndds")) {
  com <- make.comparison(grep('SCN', unique(condition), value = T),
                        ref = c("SCN_SD", "SCN_HFD"))
  scn <- apply(com, 1, function(x) results(scndds, contrast = c("condition", x[1], x[2])))
  names(scn) <- paste0(com$Var1, "vs", com$Var2)

  # Filter by p-value
  scns <- lapply(scn, function(x) subset(x, pvalue < 0.01))
  scns1 <- lapply(scn, function(x) subset(x, pvalue < 0.001))

  results_list$SCN <- list(full = scn, pval01 = scns, pval001 = scns1)

  # Save DEG tables
  sapply(names(scn), function(x) write.csv(scn[[x]],
                                          file = paste0("Tables/DEGs/", x, "_DEG.csv")))
  sapply(names(scns), function(x) write.csv(scns[[x]],
                                           file = paste0("Tables/DEGs/", x, "_DEG_pval.csv")))
  sapply(names(scns1), function(x) write.csv(scns[[x]],
                                             file = paste0("Tables/DEGs/", x, "_DEG_pval0001.csv")))
}

# Gastroc comparisons
if(exists("gasdds")) {
  com <- make.comparison(grep('Gastroc', unique(condition), value = T),
                        ref = c("Gastroc_SD", "Gastroc_HFD"))
  gastroc <- apply(com, 1, function(x) results(gasdds, contrast = c("condition", x[1], x[2])))
  names(gastroc) <- paste0(com$Var1, "vs", com$Var2)

  # Filter by p-value
  gastrocs <- lapply(gastroc, function(x) subset(x, pvalue < 0.01))
  gastrocs1 <- lapply(gastroc, function(x) subset(x, pvalue < 0.001))

  results_list$Gastroc <- list(full = gastroc, pval01 = gastrocs, pval001 = gastrocs1)

  # Save DEG tables
  sapply(names(gastroc), function(x) write.csv(gastroc[[x]],
                                               file = paste0("Tables/DEGs/", x, "_DEG.csv")))
  sapply(names(gastrocs), function(x) write.csv(gastrocs[[x]],
                                                file = paste0("Tables/DEGs/", x, "_DEG_pval.csv")))
  sapply(names(gastrocs1), function(x) write.csv(gastrocs[[x]],
                                                  file = paste0("Tables/DEGs/", x, "_DEG_pval0001.csv")))
}

# Hippo comparisons
if(exists("hipdds")) {
  com <- make.comparison(grep('Hippo', unique(condition), value = T),
                        ref = c("Hippo_SD", "Hippo_HFD"))
  hippo <- apply(com, 1, function(x) results(hipdds, contrast = c("condition", x[1], x[2])))
  names(hippo) <- paste0(com$Var1, "vs", com$Var2)

  # Filter by p-value
  hippos <- lapply(hippo, function(x) subset(x, pvalue < 0.01))
  hippos1 <- lapply(hippo, function(x) subset(x, pvalue < 0.001))

  results_list$Hippo <- list(full = hippo, pval01 = hippos, pval001 = hippos1)

  # Save DEG tables
  sapply(names(hippo), function(x) write.csv(hippo[[x]],
                                             file = paste0("Tables/DEGs/", x, "_DEG.csv")))
  sapply(names(hippos), function(x) write.csv(hippos[[x]],
                                              file = paste0("Tables/DEGs/", x, "_DEG_pval.csv")))
  sapply(names(hippos1), function(x) write.csv(hippos[[x]],
                                               file = paste0("Tables/DEGs/", x, "_DEG_pval0001.csv")))
}

# ========================================
# DEG Summary Statistics
# ========================================

cat("Generating DEG summary statistics...\n")

# Function to create DEG summary table
create_deg_summary <- function(deg_list, tissue) {
  summary_df <- data.frame(
    Comparison = names(deg_list$full),
    Total_DEGs_pval01 = sapply(deg_list$pval01, nrow),
    Total_DEGs_pval001 = sapply(deg_list$pval001, nrow),
    Upregulated_pval01 = sapply(deg_list$pval01, function(x) sum(x$log2FoldChange > 0)),
    Downregulated_pval01 = sapply(deg_list$pval01, function(x) sum(x$log2FoldChange < 0)),
    Upregulated_pval001 = sapply(deg_list$pval001, function(x) sum(x$log2FoldChange > 0)),
    Downregulated_pval001 = sapply(deg_list$pval001, function(x) sum(x$log2FoldChange < 0))
  )
  write.csv(summary_df, file = paste0("Tables/DEGs/", tissue, "_DEG_summary.csv"), row.names = FALSE)
  return(summary_df)
}

# Create summary tables
for(tissue in names(results_list)) {
  create_deg_summary(results_list[[tissue]], tissue)
}

# ========================================
# Venn Diagram Analysis
# ========================================

cat("Generating Venn diagrams...\n")

# Function to create Venn diagrams
create_venn <- function(deg_list, tissue, dietgroup, intergroup) {
  # Rename for plotting
  names(deg_list$pval01) <- gsub(paste0(tissue, '_'), '', names(deg_list$pval01))

  # Diet group Venn
  if(all(dietgroup %in% names(deg_list$pval01))) {
    venn_diet <- venndetail(lapply(deg_list$pval01[dietgroup], rownames))
    pdf(paste0("Figures/", tissue, "_Venn_diet_DEG.pdf"), width = 8, height = 8)
    plot(venn_diet)
    dev.off()

    # Save Venn details
    write.csv(getFeature(venn_diet, rlist = lapply(deg_list$pval01[dietgroup], as.data.frame)),
              file = paste0("Tables/DEGs/", tissue, "_DEGs_diet_venndetail.csv"))
  }

  # Intervention group Venn
  if(all(integroup %in% names(deg_list$pval01))) {
    venn_inte <- venndetail(lapply(deg_list$pval01[integroup], rownames))
    pdf(paste0("Figures/", tissue, "_Venn_inte_DEG.pdf"), width = 10, height = 10)
    plot(venn_inte)
    dev.off()

    # Save Venn details
    write.csv(getFeature(venn_inte, rlist = lapply(deg_list$pval01[integroup], as.data.frame)),
              file = paste0("Tables/DEGs/", tissue, "_DEGs_inte_venndetail.csv"))
  }
}

# Create Venn diagrams for each tissue
for(tissue in names(results_list)) {
  create_venn(results_list[[tissue]], tissue, dietgroup, intergroup)
}

# ========================================
# Heatmap of Shared DEGs
# ========================================

cat("Generating heatmaps of shared DEGs...\n")

# Function to create shared DEG heatmap
create_shared_heatmap <- function(deg_list, tissue, dietgroup, intergroup) {
  names(deg_list$pval01) <- gsub(paste0(tissue, '_'), '', names(deg_list$pval01))

  # Diet group heatmap
  if(all(dietgroup %in% names(deg_list$pval01))) {
    venn_diet <- venndetail(lapply(deg_list$pval01[dietgroup], rownames))
    shared_diet <- getFeature(venn_diet, rlist = lapply(deg_list$pval01[dietgroup], as.data.frame), wide = TRUE)

    mat_breaks <- seq(-4, 4, length.out = 128)

    # Heatmap for genes shared in >=2 comparisons
    if(nrow(shared_diet) > 0) {
      heatmap_data <- shared_diet %>%
        filter(SharedSets >= 2) %>%
        dplyr::select(contains("log2FoldChange"), Detail) %>%
        rename_with(~gsub('_log2FoldChange', '', .x)) %>%
        column_to_rownames(var = "Detail")

      if(nrow(heatmap_data) > 0) {
        pheatmap(heatmap_data,
                 na_col = "grey",
                 breaks = mat_breaks,
                 color = colorRampPalette(c("darkblue", "white", "red"))(127),
                 cluster_cols = FALSE,
                 fontsize_row = 3,
                 fontsize_col = 4,
                 filename = paste0("Figures/", tissue, "_diet_shared_heatmap.pdf"))
      }
    }
  }
}

# Create heatmaps for each tissue
for(tissue in names(results_list)) {
  create_shared_heatmap(results_list[[tissue]], tissue, dietgroup, intergroup)
}

cat("DEG analysis completed successfully!\n")
cat("Results saved in Tables/DEGs/ and Figures/\n")
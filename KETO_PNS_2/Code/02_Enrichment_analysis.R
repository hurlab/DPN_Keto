# ========================================
# PNS Functional Enrichment Analysis
# ========================================
# This script performs KEGG and GO enrichment analysis for DEGs

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(richR)
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(enrichplot)
  library(DOSE)
})

# Set working directory and create output folders
setwd("..")
if(!dir.exists("Figures")) dir.create("Figures")
if(!dir.exists("Tables/KEGG")) dir.create("Tables/KEGG")

# Load the data
load("ProcessedData/new_results_all.rdata")

# Define comparison groups
dietgroup <- c("HFDvsSD", "KDvsSD", "KDvsHFD")
integroup <- c("DRvsHFD", "KDIvsHFD", "EXvsHFD", "KDI_EXvsHFD")

# Build KEGG annotation for mouse
cat("Building KEGG annotation...\n")
mmko <- buildAnnot(species = "mouse", keytype = "SYMBOL", anntype = "KEGG", builtin = FALSE)

# ========================================
# KEGG Enrichment Analysis
# ========================================

cat("Performing KEGG enrichment analysis...\n")

# Function to perform KEGG enrichment
perform_kegg_enrichment <- function(deg_list, tissue) {
  # Use significant DEGs (p < 0.01)
  degs <- deg_list$pval01

  # Perform KEGG enrichment for each comparison
  kegg_results <- lapply(degs, function(x) {
    if(nrow(x) > 0) {
      tryCatch({
        richKEGG(rownames(x), mmko, builtin = FALSE)
      }, error = function(e) {
        cat("Error in KEGG enrichment:", e$message, "\n")
        return(NULL)
      })
    } else {
      return(NULL)
    }
  })

  # Filter out NULL results
  kegg_results <- kegg_results[!sapply(kegg_results, is.null)]

  # Save results
  sapply(names(kegg_results), function(x) {
    if(!is.null(kegg_results[[x]])) {
      write.csv(kegg_results[[x]], file = paste0("Tables/KEGG/", x, "_KEGG_pval.csv"))
    }
  })

  return(kegg_results)
}

# Perform KEGG enrichment for each tissue
kegg_results <- list()

if(exists("scns")) {
  cat("Performing KEGG enrichment for SCN...\n")
  # Create a list structure similar to CNS
  scn_list <- list(pval01 = scns)
  kegg_results$SCN <- perform_kegg_enrichment(scn_list, "SCN")
}

if(exists("gastrocs")) {
  cat("Performing KEGG enrichment for Gastroc...\n")
  gastroc_list <- list(pval01 = gastrocs)
  kegg_results$Gastroc <- perform_kegg_enrichment(gastroc_list, "Gastroc")
}

if(exists("hippos")) {
  cat("Performing KEGG enrichment for Hippo...\n")
  hippo_list <- list(pval01 = hippos)
  kegg_results$Hippo <- perform_kegg_enrichment(hippo_list, "Hippo")
}

# ========================================
# GSEA Analysis
# ========================================

cat("Performing GSEA analysis...\n")

# Function to perform GSEA
perform_gsea <- function(deg_full, tissue) {
  gsea_results <- lapply(deg_full, function(x) {
    tryCatch({
      parGSEA(as.data.frame(x), mmko, KEGG = TRUE)
    }, error = function(e) {
      cat("Error in GSEA:", e$message, "\n")
      return(NULL)
    })
  })

  # Filter out NULL results
  gsea_results <- gsea_results[!sapply(gsea_results, is.null)]

  # Save results
  sapply(names(gsea_results), function(x) {
    if(!is.null(gsea_results[[x]])) {
      write.csv(gsea_results[[x]], file = paste0("Tables/KEGG/", x, "_KEGG_GSEA_pval.csv"))
    }
  })

  return(gsea_results)
}

# Perform GSEA for each tissue
gsea_results <- list()

if(exists("scn")) {
  cat("Performing GSEA for SCN...\n")
  gsea_results$SCN <- perform_gsea(scn, "SCN")
}

if(exists("gastroc")) {
  cat("Performing GSEA for Gastroc...\n")
  gsea_results$Gastroc <- perform_gsea(gastroc, "Gastroc")
}

if(exists("hippo")) {
  cat("Performing GSEA for Hippo...\n")
  gsea_results$Hippo <- perform_gsea(hippo, "Hippo")
}

# ========================================
# Generate KEGG Enrichment Plots
# ========================================

cat("Generating KEGG enrichment plots...\n")

# Function to generate KEGG dot plots
generate_kegg_dotplots <- function(kegg_results, tissue, dietgroup, intergroup) {
  if(length(kegg_results) == 0) return()

  # Remove tissue prefix from names
  names(kegg_results) <- gsub(paste0(tissue, '_'), '', names(kegg_results))

  # Diet group comparisons
  diet_available <- intersect(dietgroup, names(kegg_results))
  if(length(diet_available) > 0) {
    # Basic dot plot
    p1 <- comparedot(compareResult(kegg_results[diet_available]), usePadj = FALSE) +
      theme(axis.text.x = element_text(angle = 90, size = 8),
            axis.text.y = element_text(size = 8)) +
      xlim(diet_available)

    ggsave(filename = paste0("Figures/", tissue, "_KEGG_diet_compare.pdf"),
           plot = p1, width = 10, height = 8)

    # Top 10 pathways
    top_terms <- compareResult(kegg_results[diet_available]) %>%
      group_by(group) %>%
      slice(1:10) %>%
      pull(Annot) %>%
      unique()

    p2 <- compareResult(kegg_results[diet_available]) %>%
      filter(Annot %in% top_terms) %>%
      ggplot(aes(group, Term, color = -log10(Pvalue), size = Significant)) +
      geom_point() +
      scale_colour_gradient(low = "pink", high = "red") +
      theme_light(base_size = 12) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      xlim(diet_available) +
      xlab("") + ylab("") + labs(size = "Gene Number")

    ggsave(filename = paste0("Figures/", tissue, "_KEGG_diet_top10.pdf"),
           plot = p2, width = 10, height = 8)
  }

  # Intervention group comparisons
  inte_available <- intersect(integroup, names(kegg_results))
  if(length(inte_available) > 0) {
    # Basic dot plot
    p3 <- comparedot(compareResult(kegg_results[inte_available]), usePadj = FALSE) +
      theme(axis.text.x = element_text(angle = 90, size = 8),
            axis.text.y = element_text(size = 8)) +
      xlim(inte_available)

    ggsave(filename = paste0("Figures/", tissue, "_KEGG_inte_compare.pdf"),
           plot = p3, width = 12, height = 8)

    # Top 10 pathways
    top_terms <- compareResult(kegg_results[inte_available]) %>%
      group_by(group) %>%
      slice(1:10) %>%
      pull(Annot) %>%
      unique()

    p4 <- compareResult(kegg_results[inte_available]) %>%
      filter(Annot %in% top_terms) %>%
      ggplot(aes(group, Term, color = -log10(Pvalue), size = Significant)) +
      geom_point() +
      scale_colour_gradient(low = "pink", high = "red") +
      theme_light(base_size = 12) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      xlim(inte_available) +
      xlab("") + ylab("") + labs(size = "Gene Number")

    ggsave(filename = paste0("Figures/", tissue, "_KEGG_inte_top10.pdf"),
           plot = p4, width = 12, height = 8)
  }
}

# Generate plots for each tissue
for(tissue in names(kegg_results)) {
  generate_kegg_dotplots(kegg_results[[tissue]], tissue, dietgroup, intergroup)
}

# ========================================
# Generate GSEA Plots
# ========================================

cat("Generating GSEA plots...\n")

# Function to generate GSEA dot plots
generate_gsea_dotplots <- function(gsea_results, tissue) {
  if(length(gsea_results) == 0) return()

  # Combine results
  gsea_combined <- do.call(rbind, lapply(gsea_results, result))
  gsea_combined$Group <- sub('\\..*', '', rownames(gsea_combined))
  gsea_combined$Group <- gsub(paste0(tissue, '_'), '', gsea_combined$Group)

  # Filter for common comparisons
  common_comparisons <- c("HFDvsSD", "DRvsHFD", "EXvsHFD", "EXvsDR")
  available_comparisons <- intersect(common_comparisons, unique(gsea_combined$Group))

  if(length(available_comparisons) > 0) {
    gsea_filtered <- gsea_combined %>%
      filter(Group %in% available_comparisons)

    p <- ggplot(gsea_filtered, aes(Group, pathway, color = NES, size = -log10(pval))) +
      geom_point() +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            axis.text.y = element_text(size = if(tissue == "Hippo") 10 else 6)) +
      scale_color_gradient2(low = "darkblue", high = "red", mid = "white", midpoint = 0) +
      xlab("") + ylab("") +
      xlim(available_comparisons)

    ggsave(filename = paste0("Figures/", tissue, "_GSEA_dot.pdf"),
           plot = p, width = 10, height = if(tissue == "Hippo") 8 else 6)

    # Level 2 category plots
    if("Level2" %in% colnames(gsea_filtered)) {
      for(level2_cat in unique(gsea_filtered$Level2)) {
        cat_data <- gsea_filtered %>% filter(Level2 == level2_cat)
        if(nrow(cat_data) > 0) {
          p_level2 <- ggplot(cat_data, aes(Group, pathway, color = NES, size = -log10(pval))) +
            geom_point() +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                  axis.text.y = element_text(size = 10)) +
            scale_color_gradient2(low = "darkblue", high = "red", mid = "white", midpoint = 0) +
            xlab("") + ylab("") +
            xlim(available_comparisons)

          ggsave(filename = paste0("Figures/", tissue, "_", level2_cat, "_GSEA.pdf"),
                 plot = p_level2, width = 10, height = 6)
        }
      }
    }
  }
}

# Generate GSEA plots for each tissue
for(tissue in names(gsea_results)) {
  generate_gsea_dotplots(gsea_results[[tissue]], tissue)
}

# ========================================
# GO Enrichment Analysis (Optional)
# ========================================

cat("Performing GO enrichment analysis...\n")

# Function to perform GO enrichment
perform_go_enrichment <- function(deg_list, tissue) {
  degs <- deg_list$pval01

  go_results <- lapply(degs, function(x) {
    if(nrow(x) > 0) {
      tryCatch({
        gene_list <- rownames(x)
        ego <- enrichGO(gene = gene_list,
                       OrgDb = org.Mm.eg.db,
                       keyType = "SYMBOL",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = TRUE)
        return(ego)
      }, error = function(e) {
        cat("Error in GO enrichment:", e$message, "\n")
        return(NULL)
      })
    } else {
      return(NULL)
    }
  })

  go_results <- go_results[!sapply(go_results, is.null)]

  # Save results
  sapply(names(go_results), function(x) {
    if(!is.null(go_results[[x]])) {
      write.csv(as.data.frame(go_results[[x]]),
                file = paste0("Tables/KEGG/", x, "_GO_BP_pval.csv"))
    }
  })

  return(go_results)
}

# Perform GO enrichment for each tissue (optional - can be computationally intensive)
go_results <- list()
# Uncomment the following lines if GO enrichment is needed
# if(exists("scns")) {
#   scn_list <- list(pval01 = scns)
#   go_results$SCN <- perform_go_enrichment(scn_list, "SCN")
# }

cat("Functional enrichment analysis completed successfully!\n")
cat("Results saved in Tables/KEGG/ and Figures/\n")
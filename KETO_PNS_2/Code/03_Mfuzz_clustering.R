# ========================================
# PNS Mfuzz Clustering Analysis
# ========================================
# This script performs soft clustering of gene expression patterns using Mfuzz

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(Mfuzz)
  library(Biobase)
  library(richR)
})

# Set working directory and create output folders
setwd("..")
if(!dir.exists("Figures")) dir.create("Figures")
if(!dir.exists("Tables/Mfuzz")) dir.create("Tables/Mfuzz")

# Load the data
load("ProcessedData/new_results_all.rdata")

# Build KEGG annotation for mouse
mmko <- buildAnnot(species = "mouse", keytype = "SYMBOL", anntype = "KEGG", builtin = FALSE)

# ========================================
# Function for Mfuzz Clustering
# ========================================

perform_mfuzz_clustering <- function(dds, tissue_name, n_clusters = 8) {
  cat(paste("Performing Mfuzz clustering for", tissue_name, "...\n"))

  # Get variance stabilized data
  vst_data <- assay(vst(dds))

  # Filter for variable genes (top 25% most variable)
  gene_var <- apply(vst_data, 1, var)
  var_threshold <- quantile(gene_var, 0.75)
  vst_filtered <- vst_data[gene_var > var_threshold, ]

  # Create ExpressionSet
  eset <- new("ExpressionSet", exprs = as.matrix(vst_filtered))

  # Standardize the data
  eset <- standardise(eset)

  # Estimate optimum number of clusters and fuzzification parameter
  # Calculate minimum centroid distance and Dunn index for different c values
  m_test <- seq(1.1, 2.5, by = 0.1)
  c_test <- seq(2, 12, by = 1)

  # Find optimal m (fuzzification parameter)
  min_distances <- list()
  dunn_indices <- list()

  for(c_val in c_test) {
    min_dist_c <- c()
    dunn_index_c <- c()

    for(m_val in m_test) {
      cl <- mfuzz(eset, c = c_val, m = m_val)
      min_dist_c[c(min_dist_c) + 1] <- min(distcenters(cl)[distcenters(cl) > 0])
      dunn_index_c[c(dunn_index_c) + 1] <- Dindex(cl, eset)
    }

    min_distances[[paste0("c", c_val)]] <- min_dist_c
    dunn_indices[[paste0("c", c_val)]] <- dunn_index_c
  }

  # Plot validation metrics
  # Minimum distance plots
  pdf(paste0("Figures/", tissue_name, "_minimum_distance_mfuzz.pdf"), width = 10, height = 8)
  par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))
  for(i in 1:length(c_test)) {
    if(i <= 9) {
      plot(m_test, min_distances[[i]], type = "b",
           xlab = "Fuzzification parameter m", ylab = "Minimum distance",
           main = paste("c =", c_test[i]))
    }
  }
  dev.off()

  # Dunn index plots
  pdf(paste0("Figures/", tissue_name, "_Dunn_distance_mfuzz.pdf"), width = 10, height = 8)
  par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))
  for(i in 1:length(c_test)) {
    if(i <= 9) {
      plot(m_test, dunn_indices[[i]], type = "b",
           xlab = "Fuzzification parameter m", ylab = "Dunn index",
           main = paste("c =", c_test[i]))
    }
  }
  dev.off()

  # Select optimal parameters (default: c = 8, m = 1.75)
  c_optimal <- n_clusters
  m_optimal <- 1.75

  # Perform final clustering
  cl <- mfuzz(eset, c = c_optimal, m = m_optimal)

  # Create cluster assignment matrix
  mat <- sapply(1:c_optimal, function(x) cl$membership[, x] > 0.5)
  rownames(mat) <- rownames(vst_filtered)
  colnames(mat) <- paste0("Cluster", 1:c_optimal)
  mat <- as.data.frame(mat)

  # Save cluster assignments
  write.csv(mat, file = paste0("Tables/Mfuzz/", tissue_name, "_mfuzz_cluster.csv"), row.names = TRUE)

  # Plot clusters
  pdf(paste0("Figures/", tissue_name, "_mfuzz_cluster.pdf"), width = 11.6, height = 5.2)
  mfuzz.plot(eset, cl = cl, mfrow = c(2, 4), time.labels = colnames(vst_filtered), new.window = FALSE)
  dev.off()

  # ========================================
  # KEGG Enrichment for Clusters
  # ========================================

  # Perform KEGG enrichment for each cluster
  cluster_kegg <- list()

  for(i in 1:c_optimal) {
    cluster_genes <- rownames(mat)[mat[, paste0("Cluster", i)] == TRUE]

    if(length(cluster_genes) > 0) {
      tryCatch({
        kegg_result <- richKEGG(cluster_genes, mmko, builtin = FALSE)
        cluster_kegg[[paste0("Cluster", i)]] <- kegg_result
      }, error = function(e) {
        cat("Error in KEGG enrichment for cluster", i, ":", e$message, "\n")
        cluster_kegg[[paste0("Cluster", i)]] <- NULL
      })
    }
  }

  # Combine cluster KEGG results
  if(length(cluster_kegg) > 0) {
    matrr <- do.call(rbind, lapply(cluster_kegg, function(x) {
      if(!is.null(x) && nrow(x) > 0) {
        head(x, 10)
      } else {
        return(data.frame())
      }
    }))

    # Save cluster enrichment results
    if(nrow(matrr) > 0) {
      write.csv(matrr, file = paste0("Tables/Mfuzz/", tissue_name, "_mfuzz_KEGG_cluster.csv"), row.names = FALSE)

      # Plot cluster enrichment
      matrr$Cluster <- sub('\\..*', '', rownames(matrr))

      p <- ggplot(matrr, aes(Cluster, Term, color = -log10(Pvalue), size = Significant)) +
        geom_point() +
        scale_colour_gradient(low = "pink", high = "red") +
        theme_light(base_size = 12) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.text.y = element_text(size = 8)) +
        xlab("") + ylab("") + labs(size = "Gene Number")

      ggsave(filename = paste0("Figures/", tissue_name, "_mfuzz_cluster_KEGG.pdf"),
             plot = p, width = 6.16, height = 8.99)

      # Top 10 pathways per cluster
      if(length(cluster_kegg) > 0) {
        top10_df <- do.call(rbind, lapply(names(cluster_kegg), function(cl) {
          if(!is.null(cluster_kegg[[cl]]) && nrow(cluster_kegg[[cl]]) > 0) {
            top10 <- head(cluster_kegg[[cl]], 10)
            top10$Cluster <- cl
            return(top10)
          }
        }))

        if(nrow(top10_df) > 0) {
          p2 <- ggplot(top10_df, aes(Cluster, Term, color = -log10(Pvalue), size = Significant)) +
            geom_point() +
            scale_colour_gradient(low = "pink", high = "red") +
            theme_light(base_size = 12) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                  axis.text.y = element_text(size = 8)) +
            xlab("") + ylab("") + labs(size = "Gene Number")

          ggsave(filename = paste0("Figures/", tissue_name, "_mfuzz_cluster_KEGG_top10.pdf"),
                 plot = p2, width = 6.16, height = 8.99)
        }
      }
    }
  }

  # Return results
  return(list(
    clusters = mat,
    cluster_kegg = cluster_kegg,
    membership = cl$membership,
    eset = eset,
    cl = cl
  ))
}

# ========================================
# Perform Clustering for Each Tissue
# ========================================

# SCN clustering
if(exists("scndds")) {
  scn_results <- perform_mfuzz_clustering(scndds, "SCN")
}

# Gastroc clustering
if(exists("gasdds")) {
  gas_results <- perform_mfuzz_clustering(gasdds, "Gastroc")
}

# Hippo clustering
if(exists("hipdds")) {
  hippo_results <- perform_mfuzz_clustering(hipdds, "Hippo")
}

# ========================================
# Generate Summary of Clustering Results
# ========================================

cat("Generating clustering summary...\n")

# Function to create clustering summary
create_clustering_summary <- function(results, tissue_name) {
  if(is.null(results)) return(NULL)

  cluster_summary <- data.frame(
    Cluster = colnames(results$clusters),
    Gene_Count = colSums(results$clusters)
  )

  write.csv(cluster_summary,
            file = paste0("Tables/Mfuzz/", tissue_name, "_cluster_summary.csv"),
            row.names = FALSE)

  return(cluster_summary)
}

# Create summaries for each tissue
if(exists("scn_results")) {
  create_clustering_summary(scn_results, "SCN")
}

if(exists("gas_results")) {
  create_clustering_summary(gas_results, "Gastroc")
}

if(exists("hippo_results")) {
  create_clustering_summary(hippo_results, "Hippo")
}

cat("Mfuzz clustering analysis completed successfully!\n")
cat("Results saved in Tables/Mfuzz/ and Figures/\n")
# ========================================
# Revision Analyses for Reviewer Requests
# ========================================
# - Analysis 1: Remove DRvsHFD genes, 3-way Venn (EXvsHFD, KDIvsHFD, KDI_EXvsHFD),
#               unique sets, GO/KEGG enrichment, summary heatmaps.
# - Analysis 2: SDI (SD) vs interventions (EX, KDI, KDI_EX), Venn + enrichment + heatmaps.
# - Analysis 3: Nerve-specific DEGs (sciatic minus gastroc) for intervention and maintenance schemes,
#               Venns, gene lists, GO/KEGG enrichment, heatmaps.
#
# Outputs are written to output_revision/:
#   Figures/
#   Tables/Pathway_Analysis/

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(ggplot2)
library(readr)
library(stringr)
library(VennDetail)
library(richR)
library(pheatmap)
library(ggplot2)
})

# Output directories
out_base <- "output_revision"
fig_dir <- file.path(out_base, "Figures")
path_dir <- file.path(out_base, "Tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(path_dir, recursive = TRUE, showWarnings = FALSE)

# Clear old figure outputs so only the latest variants remain
old_figs <- list.files(fig_dir, full.names = TRUE)
if (length(old_figs) > 0) unlink(old_figs, recursive = TRUE, force = TRUE)

# Common heatmap palette (light to dark red)
heat_colors <- c("#fff5f0", "#fc9272", "#cb181d")
heatmap_tops <- c(10, 20)

# Load data
load("ProcessedData/new_results_all.rdata")
if (!exists("scns")) stop("scns object missing in new_results_all.rdata")
names(scns) <- gsub("SCN_", "", names(scns))
if (exists("gastrocs")) {
  names(gastrocs) <- gsub("Gastroc_", "", names(gastrocs))
}

# Build annotations (using online KEGG/GO; builtin=FALSE as preferred)
mmko <- buildAnnot(species = "mouse", keytype = "SYMBOL", anntype = "KEGG", builtin = FALSE)
mmgo <- buildAnnot(species = "mouse", keytype = "SYMBOL", anntype = "GO", builtin = FALSE)

# Helpers ---------------------------------------------------------------
save_genes <- function(genes, label, subfolder = path_dir) {
  if (length(genes) == 0) return(invisible(NULL))
  write.csv(
    data.frame(Gene = genes),
    file = file.path(subfolder, paste0(label, ".csv")),
    row.names = FALSE
  )
}

run_kegg <- function(genes) {
  if (length(genes) == 0) return(NULL)
  tryCatch({
    res <- richKEGG(genes, mmko, builtin = FALSE)
    if (!is.null(res)) {
      res <- as.data.frame(res)
      if (!"Description" %in% names(res) && "Term" %in% names(res)) {
        res$Description <- res$Term
      }
    }
    res
  }, error = function(e) {
    message("KEGG error: ", e$message)
    NULL
  })
}

run_go <- function(genes) {
  if (length(genes) == 0) return(NULL)
  tryCatch({
    res <- richGO(genes, mmgo, ont = "BP", keytype = "SYMBOL")
    if (is.null(res)) return(NULL)
    df <- as.data.frame(res)
    if (!"Description" %in% names(df) && "Term" %in% names(df)) {
      df$Description <- df$Term
    }
    if (!"Description" %in% names(df)) return(NULL)
    df
  }, error = function(e) {
    message("GO error: ", e$message)
    NULL
  })
}

make_heatmap <- function(res_list, top_n = 15, prefix, value_col = c("Padj", "Qvalue", "Pvalue", "p.adjust", "qvalue", "pvalue")) {
  if (length(res_list) == 0) return(NULL)

  # Collect top terms across all comparisons
  all_terms <- res_list |>
    purrr::imap(function(df, nm) {
      if (is.null(df) || nrow(df) == 0) return(character(0))
      col_use <- intersect(value_col, names(df))[1]
      if (is.na(col_use) || !"Description" %in% names(df)) return(character(0))
      df |>
        dplyr::arrange(.data[[col_use]]) |>
        dplyr::slice_head(n = top_n) |>
        dplyr::pull(Description)
    }) |>
    unlist() |>
    unique()

  if (length(all_terms) == 0) return(NULL)

  mat <- matrix(NA_real_, nrow = length(all_terms), ncol = length(res_list))
  rownames(mat) <- all_terms
  colnames(mat) <- names(res_list)
  col_label <- NULL

  for (nm in names(res_list)) {
    df <- res_list[[nm]]
    if (is.null(df) || nrow(df) == 0) next
    col_use <- intersect(value_col, names(df))[1]
    if (is.na(col_use) || !"Description" %in% names(df)) next
    if (is.null(col_label)) col_label <- col_use
    vals <- df |>
      dplyr::filter(Description %in% all_terms) |>
      dplyr::select(Description, !!sym(col_use))
    for (i in seq_len(nrow(vals))) {
      mat[vals$Description[i], nm] <- -log10(vals[[col_use]][i])
    }
  }

  # Drop rows that are all NA to avoid blank plots
  keep_rows <- which(rowSums(!is.na(mat)) > 0)
  if (length(keep_rows) == 0) return(NULL)
  mat <- mat[keep_rows, , drop = FALSE]

  list(
    mat = mat,
    col_label = if (is.null(col_label)) value_col[1] else col_label,
    prefix = prefix
  )
}

write_heatmap_files <- function(hm_data, prefix, with_values = FALSE, digits = 2) {
  fill_label <- paste0("-log10(", hm_data$col_label, ")")
  df <- as.data.frame(as.table(hm_data$mat))
  colnames(df) <- c("Term", "Comparison", "Value")
  df$Term <- factor(df$Term, levels = rev(rownames(hm_data$mat)))
  df$Comparison <- factor(df$Comparison, levels = colnames(hm_data$mat))

  p <- ggplot(df, aes(x = Comparison, y = Term, fill = Value)) +
    geom_tile(color = "white") +
    scale_fill_gradientn(
      colors = heat_colors,
      na.value = "#f0f0f0",
      name = fill_label
    ) +
    labs(x = NULL, y = NULL, title = if (with_values) paste0(prefix, " (with values)") else prefix) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )

  if (with_values) {
    p <- p + geom_text(aes(label = ifelse(is.na(Value), "", round(Value, digits))), size = 3, color = "black")
  }

  pdf_path <- file.path(fig_dir, paste0(prefix, if (with_values) "_Heatmap_with_values.pdf" else "_Heatmap.pdf"))
  png_path <- file.path(fig_dir, paste0(prefix, if (with_values) "_Heatmap_with_values.png" else "_Heatmap.png"))

  ggsave(pdf_path, p, width = 9, height = 7)
  ggsave(png_path, p, width = 9, height = 7, dpi = 300)
}

save_heatmap <- function(res_list, top_n, prefix, value_col) {
  hm_data <- make_heatmap(res_list, top_n = top_n, prefix = prefix, value_col = value_col)
  if (is.null(hm_data)) {
    message("Skipping heatmap for ", prefix, ": no data to plot")
    return(invisible(NULL))
  }
  write_heatmap_files(hm_data, prefix, with_values = FALSE)
}

save_heatmap_with_values <- function(res_list, top_n, prefix, value_col, digits = 2) {
  hm_data <- make_heatmap(res_list, top_n = top_n, prefix = prefix, value_col = value_col)
  if (is.null(hm_data)) {
    message("Skipping heatmap with values for ", prefix, ": no data to plot")
    return(invisible(NULL))
  }
  write_heatmap_files(hm_data, prefix, with_values = TRUE, digits = digits)
}

save_heatmap_variants <- function(res_list, base_prefix, value_col, top_ns = heatmap_tops) {
  for (tn in top_ns) {
    prefix <- paste0(base_prefix, "_Top", tn)
    save_heatmap(res_list, top_n = tn, prefix = prefix, value_col = value_col)
    save_heatmap_with_values(res_list, top_n = tn, prefix = prefix, value_col = value_col)
  }
}

make_venn <- function(gene_lists, prefix) {
  gene_lists <- gene_lists[lengths(gene_lists) > 0]
  if (length(gene_lists) < 2) {
    message("Skipping Venn for ", prefix, ": need at least two non-empty gene sets")
    return(NULL)
  }

  venn_obj <- venndetail(gene_lists)
  plt <- plot(venn_obj, numberSize = 8, labelSize = 5)

  pdf_path <- file.path(fig_dir, paste0(prefix, "_Venn.pdf"))
  png_path <- file.path(fig_dir, paste0(prefix, "_Venn.png"))
  ggplot2::ggsave(pdf_path, plot = plt, width = 10, height = 9)
  ggplot2::ggsave(png_path, plot = plt, width = 10, height = 9, dpi = 300)

  overlaps <- VennDetail::result(venn_obj)
  write.csv(overlaps, file.path(path_dir, paste0(prefix, "_overlaps.csv")), row.names = FALSE)
  venn_obj
}

unique_only <- function(gene_lists) {
  comps <- names(gene_lists)
  purrr::map(comps, function(comp) {
    others <- setdiff(comps, comp)
    setdiff(gene_lists[[comp]], Reduce(union, gene_lists[others]))
  }) |>
    rlang::set_names(comps)
}

# Analysis 1: DR-subtracted interventions --------------------------------
dr_genes <- if ("DRvsHFD" %in% names(scns)) rownames(scns[["DRvsHFD"]]) else character(0)
interv_comps <- c("EXvsHFD", "KDIvsHFD", "KDI_EXvsHFD")
interv_avail <- interv_comps[interv_comps %in% names(scns)]
interv_filtered <- lapply(interv_avail, function(comp) setdiff(rownames(scns[[comp]]), dr_genes))
names(interv_filtered) <- interv_avail

if (length(interv_filtered) >= 2) {
  venn1 <- make_venn(interv_filtered, "Analysis1_DR_Subtracted_Interventions")

  # Save filtered gene lists
  purrr::iwalk(interv_filtered, ~save_genes(.x, paste0("Analysis1_", .y, "_Filtered_DEGs")))

  # Unique sets per comparison
  interv_unique <- unique_only(interv_filtered)
  purrr::iwalk(interv_unique, ~save_genes(.x, paste0("Analysis1_", .y, "_Unique_DEGs")))

  # Enrichment for unique sets
  interv_kegg <- purrr::imap(interv_unique, ~run_kegg(.x))
  interv_go <- purrr::imap(interv_unique, ~run_go(.x))
  purrr::iwalk(interv_kegg, ~if (!is.null(.x)) write.csv(.x, file.path(path_dir, paste0("Analysis1_", .y, "_KEGG.csv")), row.names = FALSE))
  purrr::iwalk(interv_go, ~if (!is.null(.x)) write.csv(as.data.frame(.x), file.path(path_dir, paste0("Analysis1_", .y, "_GO.csv")), row.names = FALSE))

  # Heatmaps summarizing top terms
  save_heatmap_variants(interv_go, base_prefix = "Analysis1_GO", value_col = c("Padj", "Pvalue", "p.adjust", "qvalue", "pvalue"))
  save_heatmap_variants(interv_kegg, base_prefix = "Analysis1_KEGG", value_col = c("Qvalue", "Pvalue"))
}

# Analysis 2: SD vs interventions ----------------------------------------
sdi_comps <- c("EXvsSD", "KDIvsSD", "KDI_EXvsSD")
sdi_avail <- sdi_comps[sdi_comps %in% names(scns)]
sdi_lists <- lapply(sdi_avail, function(comp) rownames(scns[[comp]]))
names(sdi_lists) <- sdi_avail

if (length(sdi_lists) >= 2) {
  venn2 <- make_venn(sdi_lists, "Analysis2_SDI_vs_Interventions")
  purrr::iwalk(sdi_lists, ~save_genes(.x, paste0("Analysis2_", .y, "_DEGs")))

  sdi_unique <- unique_only(sdi_lists)
  purrr::iwalk(sdi_unique, ~save_genes(.x, paste0("Analysis2_", .y, "_Unique_DEGs")))

  sdi_kegg <- purrr::imap(sdi_unique, ~run_kegg(.x))
  sdi_go <- purrr::imap(sdi_unique, ~run_go(.x))
  purrr::iwalk(sdi_kegg, ~if (!is.null(.x)) write.csv(.x, file.path(path_dir, paste0("Analysis2_", .y, "_KEGG.csv")), row.names = FALSE))
  purrr::iwalk(sdi_go, ~if (!is.null(.x)) write.csv(as.data.frame(.x), file.path(path_dir, paste0("Analysis2_", .y, "_GO.csv")), row.names = FALSE))

  save_heatmap_variants(sdi_go, base_prefix = "Analysis2_GO", value_col = c("Padj", "Pvalue", "p.adjust", "qvalue", "pvalue"))
  save_heatmap_variants(sdi_kegg, base_prefix = "Analysis2_KEGG", value_col = c("Qvalue", "Pvalue"))
}

# Analysis 3: Nerve-specific (SCN minus gastroc) -------------------------
nerve_specific_sets <- function(comps) {
  purrr::map(comps, function(comp) {
    if (!exists("gastrocs")) return(NULL)
    if (!(comp %in% names(scns)) || !(comp %in% names(gastrocs))) return(NULL)
    setdiff(rownames(scns[[comp]]), rownames(gastrocs[[comp]]))
  }) |>
    rlang::set_names(comps) |>
    purrr::discard(is.null)
}

# Intervention scheme nerve-specific
interv_scheme <- c("DRvsHFD", "EXvsHFD", "KDIvsHFD", "KDI_EXvsHFD")
interv_ns <- nerve_specific_sets(interv_scheme)
if (length(interv_ns) >= 2) {
  venn3 <- make_venn(interv_ns, "Analysis3_NerveSpecific_Interventions")
  purrr::iwalk(interv_ns, ~save_genes(.x, paste0("Analysis3_", .y, "_NerveSpecific_DEGs")))

  interv_ns_unique <- unique_only(interv_ns)
  purrr::iwalk(interv_ns_unique, ~save_genes(.x, paste0("Analysis3_", .y, "_NerveSpecific_Unique_DEGs")))

  interv_ns_kegg <- purrr::imap(interv_ns_unique, ~run_kegg(.x))
  interv_ns_go <- purrr::imap(interv_ns_unique, ~run_go(.x))
  purrr::iwalk(interv_ns_kegg, ~if (!is.null(.x)) write.csv(.x, file.path(path_dir, paste0("Analysis3_", .y, "_NerveSpecific_KEGG.csv")), row.names = FALSE))
  purrr::iwalk(interv_ns_go, ~if (!is.null(.x)) write.csv(as.data.frame(.x), file.path(path_dir, paste0("Analysis3_", .y, "_NerveSpecific_GO.csv")), row.names = FALSE))

  save_heatmap_variants(interv_ns_go, base_prefix = "Analysis3_Intervention_GO", value_col = c("Padj", "Pvalue", "p.adjust", "qvalue", "pvalue"))
  save_heatmap_variants(interv_ns_kegg, base_prefix = "Analysis3_Intervention_KEGG", value_col = c("Qvalue", "Pvalue"))
}

# Maintenance scheme nerve-specific
maint_scheme <- c("HFDvsSD", "KDvsSD", "KDvsHFD")
maint_ns <- nerve_specific_sets(maint_scheme)
if (length(maint_ns) >= 2) {
  venn4 <- make_venn(maint_ns, "Analysis3_NerveSpecific_Maintenance")
  purrr::iwalk(maint_ns, ~save_genes(.x, paste0("Analysis3_", .y, "_NerveSpecific_DEGs")))

  maint_ns_unique <- unique_only(maint_ns)
  purrr::iwalk(maint_ns_unique, ~save_genes(.x, paste0("Analysis3_", .y, "_NerveSpecific_Unique_DEGs")))

  maint_ns_kegg <- purrr::imap(maint_ns_unique, ~run_kegg(.x))
  maint_ns_go <- purrr::imap(maint_ns_unique, ~run_go(.x))
  purrr::iwalk(maint_ns_kegg, ~if (!is.null(.x)) write.csv(.x, file.path(path_dir, paste0("Analysis3_", .y, "_NerveSpecific_KEGG.csv")), row.names = FALSE))
  purrr::iwalk(maint_ns_go, ~if (!is.null(.x)) write.csv(as.data.frame(.x), file.path(path_dir, paste0("Analysis3_", .y, "_NerveSpecific_GO.csv")), row.names = FALSE))

  save_heatmap_variants(maint_ns_go, base_prefix = "Analysis3_Maintenance_GO", value_col = c("Padj", "Pvalue", "p.adjust", "qvalue", "pvalue"))
  save_heatmap_variants(maint_ns_kegg, base_prefix = "Analysis3_Maintenance_KEGG", value_col = c("Qvalue", "Pvalue"))
}

cat("Revision analyses complete. Outputs in ", out_base, "\n")

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
})

# Output directories
out_base <- "output_revision"
fig_dir <- file.path(out_base, "Figures")
path_dir <- file.path(out_base, "Tables", "Pathway_Analysis")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(path_dir, recursive = TRUE, showWarnings = FALSE)

# Load data
load("ProcessedData/new_results_all.rdata")
if (!exists("scns")) stop("scns object missing in new_results_all.rdata")
names(scns) <- gsub("SCN_", "", names(scns))
if (exists("gastrocs")) {
  names(gastrocs) <- gsub("Gastroc_", "", names(gastrocs))
}

# Build annotations (use builtin to avoid network dependency)
mmko <- buildAnnot(species = "mouse", keytype = "SYMBOL", anntype = "KEGG", builtin = TRUE)
mmgo <- buildAnnot(species = "mouse", keytype = "SYMBOL", anntype = "GO", builtin = TRUE)

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
    res <- richKEGG(genes, mmko, builtin = TRUE)
    if (!is.null(res)) res <- as.data.frame(res)
    res
  }, error = function(e) {
    message("KEGG error: ", e$message)
    NULL
  })
}

run_go <- function(genes) {
  NULL
}

make_heatmap <- function(res_list, top_n = 15, prefix, value_col = c("Qvalue", "Pvalue")) {
  if (length(res_list) == 0) return(invisible(NULL))

  # Collect top terms across all comparisons
  all_terms <- res_list |>
    purrr::imap(function(df, nm) {
      if (is.null(df) || nrow(df) == 0) return(character(0))
      col_use <- intersect(value_col, names(df))[1]
      if (is.na(col_use)) return(character(0))
      df |>
        dplyr::arrange(.data[[col_use]]) |>
        dplyr::slice_head(n = top_n) |>
        dplyr::pull(Description)
    }) |>
    unlist() |>
    unique()

  if (length(all_terms) == 0) return(invisible(NULL))

  mat <- matrix(NA_real_, nrow = length(all_terms), ncol = length(res_list))
  rownames(mat) <- all_terms
  colnames(mat) <- names(res_list)

  for (nm in names(res_list)) {
    df <- res_list[[nm]]
    if (is.null(df) || nrow(df) == 0) next
    col_use <- intersect(value_col, names(df))[1]
    if (is.na(col_use)) next
    vals <- df |>
      dplyr::filter(Description %in% all_terms) |>
      dplyr::select(Description, !!sym(col_use))
    for (i in seq_len(nrow(vals))) {
      mat[vals$Description[i], nm] <- -log10(vals[[col_use]][i])
    }
  }

  pheatmap(
    mat,
    color = colorRampPalette(c("#f7fbff", "#6baed6", "#08306b"))(50),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    na_col = "#f0f0f0",
    fontsize_row = 8,
    fontsize_col = 10,
    main = prefix
  )
}

save_heatmap <- function(res_list, top_n, prefix, value_col) {
  pdf(file.path(fig_dir, paste0(prefix, "_Heatmap.pdf")), width = 9, height = 7)
  make_heatmap(res_list, top_n = top_n, prefix = prefix, value_col = value_col)
  dev.off()
}

make_venn <- function(gene_lists, prefix) {
  if (length(gene_lists) < 2) return(NULL)
  venn_obj <- venndetail(gene_lists)
  pdf(file.path(fig_dir, paste0(prefix, "_Venn.pdf")), width = 10, height = 9)
  plot(venn_obj)
  dev.off()
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
  save_heatmap(interv_go, top_n = 20, prefix = "Analysis1_GO", value_col = c("p.adjust", "qvalue", "pvalue"))
  save_heatmap(interv_kegg, top_n = 20, prefix = "Analysis1_KEGG", value_col = c("Qvalue", "Pvalue"))
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

  save_heatmap(sdi_go, top_n = 20, prefix = "Analysis2_GO", value_col = c("p.adjust", "qvalue", "pvalue"))
  save_heatmap(sdi_kegg, top_n = 20, prefix = "Analysis2_KEGG", value_col = c("Qvalue", "Pvalue"))
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

  save_heatmap(interv_ns_go, top_n = 15, prefix = "Analysis3_Intervention_GO", value_col = c("p.adjust", "qvalue", "pvalue"))
  save_heatmap(interv_ns_kegg, top_n = 15, prefix = "Analysis3_Intervention_KEGG", value_col = c("Qvalue", "Pvalue"))
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

  save_heatmap(maint_ns_go, top_n = 15, prefix = "Analysis3_Maintenance_GO", value_col = c("p.adjust", "qvalue", "pvalue"))
  save_heatmap(maint_ns_kegg, top_n = 15, prefix = "Analysis3_Maintenance_KEGG", value_col = c("Qvalue", "Pvalue"))
}

cat("Revision analyses complete. Outputs in ", out_base, "\n")

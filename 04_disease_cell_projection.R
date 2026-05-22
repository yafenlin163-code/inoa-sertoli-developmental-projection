#!/usr/bin/env Rscript

# ==============================================================================
# 04_disease_cell_projection.R
#
# Purpose
#   Project disease-derived cells from idiopathic non-obstructive azoospermia
#   and azoospermia factor a deletion samples onto normal cell-type-specific
#   developmental references generated in 03_developmental_reference_construction.
#   The script produces:
#     1. cell-level projected pseudotime and maturation estimates;
#     2. sample-level and disease-group summaries;
#     3. projection quality-control tables;
#     4. cell-style main and supplementary figures.
#
# Interpretation boundary
#   This analysis tests whether disease cells resemble earlier normal
#   developmental states. Disease cells are mapped onto normal references rather
#   than used to construct the reference trajectory. Results support a
#   developmental-reference association model, not causal proof of developmental
#   arrest.
#
# Expected upstream inputs
#   results/03_developmental_reference_construction/<celltype>/rds/
#     seurat_with_pseudotime_PAGA_DPT.rds
#   Legacy Part 3 output directory names are also supported.
#
#   Disease objects, discovered automatically when present:
#     NOA_results/iNOA_merged_final.rds
#     NOA_results/iNOA_merged.rds
#     iNOA_merged_final.rds
#     iNOA_merged.rds
#     NOA_results/AZFa_Del_1_processed.rds
#     NOA_results/AZFa_Del_1_final.rds
#     AZFa_Del_1_processed.rds
#     AZFa_Del_1_final.rds
#
# Main outputs
#   results/04_disease_cell_projection/
#     tables/disease_projected_pseudotime_values.csv
#     tables/disease_projection_summary_by_sample.csv
#     tables/disease_projection_summary_by_group.csv
#     tables/disease_projection_quality_by_celltype.csv
#     tables/part4_figure_manifest.csv
#     rds/combined_final_with_projections.rds
#     figures/main/*.pdf
#     figures/supplementary/*.pdf
#
# Default paths can be overridden with environment variables:
#   PROJECT_DIR, PART4_BASE_DIR, PART4_PART3_DIR, PART4_OUTPUT_DIR,
#   PART4_DISEASE_RDS, PART4_iNOA_RDS, PART4_AZFA_RDS, PART4_CELLTYPE_COL
#
# Recommended figure placement
#   Main:
#     Fig5A_multicelltype_iNOA_developmental_deficit_priority.pdf
#     Fig5B_Sertoli_projection_on_normal_axis.pdf
#     Fig5C_Sertoli_pseudotime_deficit_by_sample.pdf
#     Fig5D_projection_workflow_schematic.pdf
#   Supplementary:
#     FigS5A_all_celltypes_projection_axis.pdf
#     FigS5B_cell_level_projected_pseudotime_density.pdf
#     FigS5C_projection_quality_by_celltype.pdf
#     FigS5D_nearest_reference_age_composition.pdf
#     FigS5E_disease_celltype_composition_for_projection.pdf
#     FigS5F_celltype_annotation_score_qc.pdf
#     FigS5G_Sertoli_projection_quality_qc.pdf
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

HAS_PATCHWORK <- requireNamespace("patchwork", quietly = TRUE)
HAS_GGREPEL <- requireNamespace("ggrepel", quietly = TRUE)
HAS_IRLBA <- requireNamespace("irlba", quietly = TRUE)
HAS_FNN <- requireNamespace("FNN", quietly = TRUE)
HAS_RANN <- requireNamespace("RANN", quietly = TRUE)

# ------------------------------------------------------------------------------
# 0. Configuration
# ------------------------------------------------------------------------------

USER_BASE_DIR <- Sys.getenv("PART4_BASE_DIR", unset = Sys.getenv("PROJECT_DIR", unset = ""))
USER_DISEASE_RDS <- Sys.getenv("PART4_DISEASE_RDS", unset = "")
USER_iNOA_RDS <- Sys.getenv("PART4_iNOA_RDS", unset = "")
USER_AZFA_RDS <- Sys.getenv("PART4_AZFA_RDS", unset = "")
USER_CELLTYPE_COL <- Sys.getenv("PART4_CELLTYPE_COL", unset = "")
FORCE_MARKER_ANNOTATION <- tolower(Sys.getenv("PART4_FORCE_MARKER_ANNOTATION", unset = "FALSE")) %in%
  c("1", "true", "yes", "y")

OUTDIR_NAME <- Sys.getenv("PART4_OUTPUT_DIR", unset = file.path("results", "04_disease_cell_projection"))
PART3_DIR_NAME <- Sys.getenv(
  "PART4_PART3_DIR",
  unset = file.path("results", "03_developmental_reference_construction")
)
LEGACY_PART3_DIR_NAME <- "part3_multicelltype_pseudotime_results"
LEGACY_PART4_DIR_NAME <- "part4_disease_projection_results"

CELLTYPES_TO_PROJECT <- c(
  "Sertoli", "Leydig", "PMC", "Endothelial", "Macrophage",
  "T_cell", "Spermatogonia", "Spermatocyte", "Spermatid"
)

CELLTYPE_DIR_IDS <- c(
  Sertoli       = "sertoli",
  Leydig        = "leydig",
  PMC           = "pmc",
  Endothelial   = "endothelial",
  Macrophage    = "macrophage",
  T_cell        = "t_cell",
  Spermatogonia = "spermatogonia",
  Spermatocyte  = "spermatocyte",
  Spermatid     = "spermatid"
)

ASSAY_USE <- "RNA"
SLOT_USE <- "data"
AGE_COL <- "age_numeric"
SAMPLE_COL_CANDIDATES <- c("sample_id", "sample", "orig.ident", "donor_id", "patient_id")
DISEASE_COL_CANDIDATES <- c("disease_group", "disease_type", "condition", "group", "diagnosis")
CELLTYPE_COL_CANDIDATES <- c("cell_type", "celltype", "cell_type_final", "annotation", "major_celltype")

N_PCS_PROJECTION <- 30
N_FEATURES_PROJECTION <- 2500
K_NEIGHBORS <- 30
MIN_REF_CELLS <- 50
MIN_DISEASE_CELLS <- 20
MAX_REF_CELLS_FOR_PCA <- 9000
MAX_DISEASE_CELLS_FOR_PLOTS <- 5000
WASSERSTEIN_GRID_N <- 200
PROJECTION_DISTANCE_QC_QUANTILE <- 0.95
set.seed(20260506)

# Cell-style palette adapted to the user's visual direction.
PAL <- list(
  pink_dark = "#C36E8C",
  pink_mid = "#D89FB3",
  pink_light = "#F1D9E0",
  pink_pale = "#FBEEF1",
  blue_dark = "#245399",
  blue_mid = "#3772A9",
  blue_light = "#99B6D2",
  blue_lavender = "#6B7AB2",
  grey_dark = "#7F7F7F",
  grey_mid = "#D2D2D2",
  grey_light = "#E7E5EE",
  lavender = "#A497BA",
  black = "#000000",
  white = "#FFFFFF"
)

disease_palette <- c(
  "iNOA" = PAL$pink_dark,
  "AZFa_Del" = PAL$blue_mid,
  "Disease" = PAL$pink_mid,
  "Unknown" = PAL$grey_dark
)

celltype_palette <- c(
  Sertoli       = PAL$pink_dark,
  Leydig        = "#C08A4B",
  PMC           = "#9A9087",
  Endothelial   = PAL$blue_light,
  Macrophage    = "#B0186B",
  T_cell        = "#E6B7C6",
  Spermatogonia = "#BFDDE3",
  Spermatocyte  = "#B78398",
  Spermatid     = "#D95F52"
)

CELLTYPE_MARKERS <- list(
  Sertoli       = c("SOX9", "AMH", "CLDN11", "WT1", "FSHR", "GATA4", "VIM", "INHA", "KITLG"),
  Leydig        = c("CYP11A1", "CYP17A1", "STAR", "HSD3B1", "INSL3", "NR5A1", "LHCGR"),
  PMC           = c("ACTA2", "MYH11", "MYOCD", "DES", "CNN1", "TAGLN"),
  Endothelial   = c("PECAM1", "CD34", "VWF", "ENG", "CDH5", "KDR"),
  Macrophage    = c("CD68", "CD163", "MRC1", "CSF1R", "LST1", "TYROBP"),
  T_cell        = c("CD3D", "CD3E", "CD3G", "TRAC", "CD8A", "CD4"),
  Spermatogonia = c("MAGEA4", "DAZL", "STRA8", "ID4", "FGFR3", "GFRA1", "UTF1"),
  Spermatocyte  = c("SYCP1", "SYCP3", "MLH3", "SPO11", "MEIOC", "HORMAD1"),
  Spermatid     = c("TNP1", "TNP2", "PRM1", "PRM2", "ACR", "ACRV1", "SPEM1")
)

MATURE_MARKER_CONFIG <- list(
  Sertoli = list(
    immature = c("AMH", "JUN", "FOS", "FOSB", "EGR1", "EGR3", "NR4A1", "BTG2", "DUSP1"),
    mature = c("CLDN11", "TJP1", "GATA4", "SOX9", "WT1", "HOPX", "INHA", "INHBB", "KITLG", "FSHR")
  ),
  Leydig = list(
    immature = c("TCF21", "PDGFRA", "DLK1", "NR2F2", "LIFR"),
    mature = c("STAR", "CYP11A1", "CYP17A1", "HSD3B1", "INSL3", "LHCGR", "NR5A1")
  ),
  PMC = list(
    immature = c("PDGFRA", "LUM", "DCN", "COL1A1", "COL3A1"),
    mature = c("ACTA2", "MYH11", "MYOCD", "CNN1", "TAGLN", "DES")
  ),
  Endothelial = list(
    immature = c("KDR", "SOX17", "ESAM", "EMCN", "PLVAP"),
    mature = c("PECAM1", "VWF", "CDH5", "ENG", "KLF2", "KLF4")
  ),
  Macrophage = list(
    immature = c("LYZ", "S100A8", "S100A9", "IL1B", "CXCL8"),
    mature = c("CD68", "CD163", "MRC1", "CSF1R", "MAFB", "SPI1")
  ),
  T_cell = list(
    immature = c("TCF7", "LEF1", "IL7R", "CCR7"),
    mature = c("CD3D", "CD3E", "TRAC", "TBX21", "GZMK", "GZMB")
  ),
  Spermatogonia = list(
    immature = c("ID4", "FGFR3", "GFRA1", "UTF1", "ZBTB16"),
    mature = c("MAGEA4", "DAZL", "STRA8", "DMRT1", "KIT")
  ),
  Spermatocyte = list(
    immature = c("STRA8", "DMC1", "SPO11", "SYCP1"),
    mature = c("SYCP3", "HORMAD1", "MEIOC", "MLH3", "PIWIL1")
  ),
  Spermatid = list(
    immature = c("ACR", "ACRV1", "SPATA18", "TNP1"),
    mature = c("TNP2", "PRM1", "PRM2", "ODF1", "AKAP4")
  )
)

# ------------------------------------------------------------------------------
# 1. Path discovery and output folders
# ------------------------------------------------------------------------------

rel_exists <- function(base, rel_paths) {
  any(file.exists(file.path(base, rel_paths)))
}

is_absolute_path <- function(path) {
  grepl("^[A-Za-z]:[/\\\\]", path) || startsWith(path, "/") || startsWith(path, "\\\\")
}

get_script_dir <- function() {
  frames <- sys.frames()
  ofiles <- vapply(frames, function(x) {
    if (!is.null(x$ofile)) as.character(x$ofile) else NA_character_
  }, character(1))
  ofiles <- ofiles[!is.na(ofiles) & nzchar(ofiles)]
  if (length(ofiles) > 0) {
    return(dirname(normalizePath(tail(ofiles, 1), winslash = "/", mustWork = FALSE)))
  }
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = FALSE)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active_file <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (nzchar(active_file)) {
      return(dirname(normalizePath(active_file, winslash = "/", mustWork = FALSE)))
    }
  }
  NA_character_
}

parent_dirs <- function(path) {
  if (is.na(path) || !nzchar(path) || !dir.exists(path)) return(character())
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  out <- character()
  repeat {
    out <- c(out, path)
    parent <- dirname(path)
    if (identical(parent, path)) break
    path <- parent
  }
  out
}

detect_base_dir <- function() {
  part3_rel <- c(
    file.path(PART3_DIR_NAME, "sertoli", "rds", "seurat_with_pseudotime_PAGA_DPT.rds"),
    file.path(LEGACY_PART3_DIR_NAME, "sertoli", "rds", "seurat_with_pseudotime_PAGA_DPT.rds")
  )
  disease_rel <- c(
    "NOA_results/iNOA_merged_final.rds",
    "NOA_results/iNOA_merged.rds",
    "iNOA_merged_final.rds",
    "iNOA_merged.rds",
    "NOA_results/AZFa_Del_1_processed.rds",
    "NOA_results/AZFa_Del_1_final.rds",
    "AZFa_Del_1_processed.rds",
    "AZFa_Del_1_final.rds"
  )
  explicit_dirs <- dirname(c(USER_DISEASE_RDS, USER_iNOA_RDS, USER_AZFA_RDS))
  explicit_dirs <- explicit_dirs[nzchar(explicit_dirs) & explicit_dirs != "."]
  seeds <- unique(c(USER_BASE_DIR, explicit_dirs, getwd(), get_script_dir()))
  seeds <- seeds[!is.na(seeds) & nzchar(seeds)]
  seeds <- seeds[dir.exists(seeds)]
  candidates <- unique(unlist(lapply(seeds, parent_dirs), use.names = FALSE))
  for (base in candidates) {
    if (rel_exists(base, part3_rel) && rel_exists(base, disease_rel)) return(base)
  }
  for (base in candidates) {
    if (dir.exists(file.path(base, PART3_DIR_NAME)) || dir.exists(file.path(base, LEGACY_PART3_DIR_NAME))) return(base)
  }
  stop(
    "Could not detect project base directory. Set working directory to the ",
    "project root or set Sys.setenv(PART4_BASE_DIR = '/path/to/project/root')."
  )
}

BASE_DIR <- normalizePath(detect_base_dir(), winslash = "/", mustWork = TRUE)
PART3_DIR <- file.path(BASE_DIR, PART3_DIR_NAME)
if (!dir.exists(PART3_DIR) && dir.exists(file.path(BASE_DIR, LEGACY_PART3_DIR_NAME))) {
  PART3_DIR <- file.path(BASE_DIR, LEGACY_PART3_DIR_NAME)
}
OUTDIR <- if (is_absolute_path(OUTDIR_NAME)) OUTDIR_NAME else file.path(BASE_DIR, OUTDIR_NAME)
TABLE_DIR <- file.path(OUTDIR, "tables")
RDS_DIR <- file.path(OUTDIR, "rds")
FIG_DIR <- file.path(OUTDIR, "figures")
MAIN_FIG_DIR <- file.path(FIG_DIR, "main")
SUPP_FIG_DIR <- file.path(FIG_DIR, "supplementary")
REPORT_DIR <- file.path(OUTDIR, "reports")

dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RDS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(MAIN_FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(SUPP_FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(REPORT_DIR, showWarnings = FALSE, recursive = TRUE)

message("Part 4 base directory: ", BASE_DIR)
message("Part 4 output directory: ", OUTDIR)

# ------------------------------------------------------------------------------
# 2. General helpers
# ------------------------------------------------------------------------------

safe_celltype_id <- function(x) {
  if (x %in% names(CELLTYPE_DIR_IDS)) return(unname(CELLTYPE_DIR_IDS[[x]]))
  tolower(gsub("[^A-Za-z0-9]+", "_", x))
}

first_existing_col <- function(candidates, df) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) NA_character_ else hit[[1]]
}

to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

fmt_num <- function(x, digits = 3) {
  ifelse(is.finite(x), formatC(x, digits = digits, format = "f"), "NA")
}

theme_cell <- function(base_size = 8.8) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(color = PAL$black),
      plot.title = element_text(face = "bold", size = base_size + 2.2, hjust = 0),
      plot.subtitle = element_text(size = base_size, color = PAL$grey_dark, hjust = 0),
      axis.title = element_text(face = "bold", size = base_size + 0.4),
      axis.text = element_text(color = PAL$black, size = base_size),
      axis.line = element_line(linewidth = 0.32, color = PAL$black),
      axis.ticks = element_line(linewidth = 0.28, color = PAL$black),
      legend.title = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size),
      legend.key.size = grid::unit(3.6, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size + 0.4),
      panel.background = element_rect(fill = PAL$white, color = NA),
      plot.background = element_rect(fill = PAL$white, color = NA),
      plot.margin = margin(5, 6, 5, 5)
    )
}

save_fig <- function(name, plot, width, height, section = c("main", "supplementary")) {
  section <- match.arg(section)
  out_dir <- if (section == "main") MAIN_FIG_DIR else SUPP_FIG_DIR
  out_file <- file.path(out_dir, name)
  tryCatch(
    ggsave(out_file, plot = plot, width = width, height = height, device = cairo_pdf, bg = "white"),
    error = function(e) ggsave(out_file, plot = plot, width = width, height = height, bg = "white")
  )
  message("Saved ", section, " figure: ", normalizePath(out_file, winslash = "/", mustWork = FALSE))
  invisible(out_file)
}

safe_fwrite <- function(x, file) {
  fwrite(x, file)
  message("Saved table: ", normalizePath(file, winslash = "/", mustWork = FALSE))
}

join_layers_safe <- function(obj, assay = ASSAY_USE) {
  if (!assay %in% names(obj@assays)) return(obj)
  join_fun <- NULL
  if (exists("JoinLayers", where = asNamespace("Seurat"), mode = "function")) {
    join_fun <- get("JoinLayers", envir = asNamespace("Seurat"))
  } else if (exists("JoinLayers", where = asNamespace("SeuratObject"), mode = "function")) {
    join_fun <- get("JoinLayers", envir = asNamespace("SeuratObject"))
  }
  if (is.null(join_fun)) return(obj)
  tryCatch(join_fun(obj, assay = assay), error = function(e) obj)
}

get_assay_data_safe <- function(obj, assay = ASSAY_USE, slot = SLOT_USE) {
  assay_obj <- obj[[assay]]
  layer_names <- tryCatch(SeuratObject::Layers(assay_obj), error = function(e) character(0))
  if (length(layer_names) > 0) {
    exact <- layer_names[layer_names == slot]
    if (length(exact) == 1) return(SeuratObject::LayerData(assay_obj, layer = exact))
    prefixed <- layer_names[startsWith(layer_names, paste0(slot, "."))]
    if (length(prefixed) > 0) {
      mats <- lapply(prefixed, function(ly) SeuratObject::LayerData(assay_obj, layer = ly))
      genes <- Reduce(union, lapply(mats, rownames))
      mats <- lapply(mats, function(m) {
        missing <- setdiff(genes, rownames(m))
        if (length(missing) > 0) {
          zero <- Matrix::Matrix(0, nrow = length(missing), ncol = ncol(m), sparse = TRUE)
          rownames(zero) <- missing
          colnames(zero) <- colnames(m)
          m <- rbind(m, zero)
        }
        m[genes, , drop = FALSE]
      })
      return(do.call(cbind, mats))
    }
  }
  if ("layer" %in% names(formals(Seurat::GetAssayData))) {
    tryCatch(
      Seurat::GetAssayData(obj, assay = assay, layer = slot),
      error = function(e) Seurat::GetAssayData(obj, assay = assay, slot = slot)
    )
  } else {
    Seurat::GetAssayData(obj, assay = assay, slot = slot)
  }
}

ensure_normalized <- function(obj, assay = ASSAY_USE) {
  DefaultAssay(obj) <- assay
  obj <- join_layers_safe(obj, assay)
  data_ok <- tryCatch({
    m <- get_assay_data_safe(obj, assay, "data")
    nrow(m) > 0 && ncol(m) > 0 && length(m@x) > 0
  }, error = function(e) FALSE)
  if (!data_ok) obj <- NormalizeData(obj, assay = assay, verbose = FALSE)
  obj
}

scale_rows_by_reference <- function(ref_expr, query_expr, genes) {
  ref_mat <- as.matrix(ref_expr[genes, , drop = FALSE])
  query_mat <- as.matrix(query_expr[genes, , drop = FALSE])
  mu <- rowMeans(ref_mat, na.rm = TRUE)
  s <- apply(ref_mat, 1, sd, na.rm = TRUE)
  keep <- is.finite(mu) & is.finite(s) & s > 0
  ref_z <- sweep(ref_mat[keep, , drop = FALSE], 1, mu[keep], "-")
  ref_z <- sweep(ref_z, 1, s[keep], "/")
  query_z <- sweep(query_mat[keep, , drop = FALSE], 1, mu[keep], "-")
  query_z <- sweep(query_z, 1, s[keep], "/")
  ref_z[ref_z > 10] <- 10
  ref_z[ref_z < -10] <- -10
  query_z[query_z > 10] <- 10
  query_z[query_z < -10] <- -10
  list(ref_z = ref_z, query_z = query_z, genes = rownames(ref_z))
}

select_balanced_cells <- function(meta, max_cells, age_col = AGE_COL) {
  if (nrow(meta) <= max_cells) return(rownames(meta))
  if (age_col %in% colnames(meta)) {
    split_idx <- split(rownames(meta), as.character(meta[[age_col]]))
    n_each <- ceiling(max_cells / length(split_idx))
    cells <- unlist(lapply(split_idx, function(v) sample(v, min(length(v), n_each))), use.names = FALSE)
    if (length(cells) > max_cells) cells <- sample(cells, max_cells)
    return(cells)
  }
  sample(rownames(meta), max_cells)
}

run_reference_pca <- function(ref_z, query_z, n_pcs = N_PCS_PROJECTION) {
  n_pcs <- min(n_pcs, nrow(ref_z) - 1, ncol(ref_z) - 1)
  if (n_pcs < 2) stop("Too few genes/cells for PCA projection.")
  ref_t <- t(ref_z)
  if (HAS_IRLBA && nrow(ref_t) > n_pcs + 1 && ncol(ref_t) > n_pcs + 1) {
    pca <- irlba::prcomp_irlba(ref_t, n = n_pcs, center = FALSE, scale. = FALSE)
  } else {
    pca <- prcomp(ref_t, center = FALSE, scale. = FALSE, rank. = n_pcs)
  }
  ref_pca <- pca$x[, 1:n_pcs, drop = FALSE]
  query_pca <- t(query_z) %*% pca$rotation[, 1:n_pcs, drop = FALSE]
  colnames(ref_pca) <- paste0("PC", seq_len(ncol(ref_pca)))
  colnames(query_pca) <- colnames(ref_pca)
  list(ref_pca = ref_pca, query_pca = as.matrix(query_pca), rotation = pca$rotation[, 1:n_pcs, drop = FALSE])
}

find_knn <- function(query, ref, k = K_NEIGHBORS) {
  k <- min(k, nrow(ref))
  if (HAS_FNN) {
    out <- FNN::get.knnx(data = ref, query = query, k = k)
    return(list(idx = out$nn.index, dist = out$nn.dist))
  }
  if (HAS_RANN) {
    out <- RANN::nn2(data = ref, query = query, k = k)
    return(list(idx = out$nn.idx, dist = out$nn.dists))
  }
  idx <- matrix(NA_integer_, nrow = nrow(query), ncol = k)
  dist <- matrix(NA_real_, nrow = nrow(query), ncol = k)
  for (i in seq_len(nrow(query))) {
    d <- sqrt(rowSums(sweep(ref, 2, query[i, ], "-")^2))
    ord <- order(d)[seq_len(k)]
    idx[i, ] <- ord
    dist[i, ] <- d[ord]
  }
  list(idx = idx, dist = dist)
}

weights_from_distance <- function(dist_row) {
  d <- as.numeric(dist_row)
  d[!is.finite(d)] <- max(d[is.finite(d)], na.rm = TRUE)
  sigma <- median(d[d > 0], na.rm = TRUE)
  if (!is.finite(sigma) || sigma <= 0) sigma <- max(d, na.rm = TRUE)
  if (!is.finite(sigma) || sigma <= 0) return(rep(1 / length(d), length(d)))
  w <- exp(-(d^2) / (2 * sigma^2))
  if (sum(w) <= 0 || !is.finite(sum(w))) return(rep(1 / length(d), length(d)))
  w / sum(w)
}

weighted_mean_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w)
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

weighted_sd_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w)
  if (sum(ok) < 2) return(NA_real_)
  m <- weighted_mean_safe(x[ok], w[ok])
  sqrt(sum(w[ok] * (x[ok] - m)^2) / sum(w[ok]))
}

weight_entropy <- function(w) {
  w <- w[is.finite(w) & w > 0]
  if (length(w) == 0) return(NA_real_)
  -sum(w * log(w)) / log(length(w))
}

wasserstein_1d <- function(x, y, grid_n = WASSERSTEIN_GRID_N) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) == 0 || length(y) == 0) return(NA_real_)
  probs <- seq(0, 1, length.out = grid_n)
  qx <- quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 8)
  qy <- quantile(y, probs = probs, na.rm = TRUE, names = FALSE, type = 8)
  mean(abs(qx - qy), na.rm = TRUE)
}

nearest_reference_age <- function(pt, age_summary) {
  vapply(pt, function(z) {
    if (!is.finite(z)) return(NA_real_)
    age_summary$age_numeric[which.min(abs(age_summary$median_pseudotime - z))]
  }, numeric(1))
}

infer_reference_age <- function(pt, age_summary) {
  age_summary <- age_summary %>%
    arrange(median_pseudotime) %>%
    distinct(median_pseudotime, .keep_all = TRUE)
  x <- age_summary$median_pseudotime
  y <- age_summary$age_numeric
  vapply(pt, function(z) {
    if (!is.finite(z)) return(NA_real_)
    if (length(x) < 2 || length(unique(x)) < 2) {
      return(y[which.min(abs(x - z))])
    }
    approx(x = x, y = y, xout = z, rule = 2, ties = "ordered")$y
  }, numeric(1))
}

score_marker_set <- function(expr, genes) {
  rn <- rownames(expr)
  lookup <- setNames(rn, toupper(rn))
  matched <- unique(unname(lookup[toupper(genes)]))
  matched <- matched[!is.na(matched)]
  if (length(matched) == 0) return(rep(NA_real_, ncol(expr)))
  Matrix::colMeans(expr[matched, , drop = FALSE])
}

zscore <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

# ------------------------------------------------------------------------------
# 3. Load normal references and disease objects
# ------------------------------------------------------------------------------

reference_file_for_celltype <- function(celltype) {
  file.path(
    PART3_DIR,
    safe_celltype_id(celltype),
    "rds",
    "seurat_with_pseudotime_PAGA_DPT.rds"
  )
}

load_reference <- function(celltype) {
  f <- reference_file_for_celltype(celltype)
  if (!file.exists(f)) {
    warning("Missing Part 3 reference for ", celltype, ": ", f)
    return(NULL)
  }
  obj <- readRDS(f)
  DefaultAssay(obj) <- ASSAY_USE
  obj <- ensure_normalized(obj, ASSAY_USE)
  if (!"pseudotime_final" %in% colnames(obj@meta.data)) {
    warning("Reference lacks pseudotime_final for ", celltype)
    return(NULL)
  }
  if (!AGE_COL %in% colnames(obj@meta.data)) {
    warning("Reference lacks age_numeric for ", celltype)
    return(NULL)
  }
  obj$part4_reference_celltype <- celltype
  obj
}

find_first_file <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) NA_character_ else normalizePath(hit[[1]], winslash = "/", mustWork = TRUE)
}

disease_candidates <- list(
  combined = c(
    USER_DISEASE_RDS,
    file.path(BASE_DIR, "NOA_results", "combined_disease_final.rds"),
    file.path(BASE_DIR, "combined_disease_final.rds"),
    file.path(BASE_DIR, "part4_disease_projection_results", "rds", "combined_final_with_projections.rds")
  ),
  iNOA = c(
    USER_iNOA_RDS,
    file.path(BASE_DIR, "NOA_results", "iNOA_merged_final.rds"),
    file.path(BASE_DIR, "NOA_results", "iNOA_merged.rds"),
    file.path(BASE_DIR, "iNOA_merged_final.rds"),
    file.path(BASE_DIR, "iNOA_merged.rds")
  ),
  AZFa_Del = c(
    USER_AZFA_RDS,
    file.path(BASE_DIR, "NOA_results", "AZFa_Del_1_processed.rds"),
    file.path(BASE_DIR, "NOA_results", "AZFa_Del_1_final.rds"),
    file.path(BASE_DIR, "AZFa_Del_1_processed.rds"),
    file.path(BASE_DIR, "AZFa_Del_1_final.rds")
  )
)

standardize_disease_metadata <- function(obj, default_group) {
  DefaultAssay(obj) <- ASSAY_USE
  obj <- ensure_normalized(obj, ASSAY_USE)
  meta <- obj@meta.data
  sample_col <- first_existing_col(SAMPLE_COL_CANDIDATES, meta)
  disease_col <- first_existing_col(DISEASE_COL_CANDIDATES, meta)
  if (is.na(sample_col)) {
    obj$sample_id <- if ("orig.ident" %in% colnames(meta)) as.character(obj$orig.ident) else default_group
  } else {
    obj$sample_id <- as.character(meta[[sample_col]])
  }
  if (is.na(disease_col)) {
    obj$disease_group <- default_group
  } else {
    obj$disease_group <- as.character(meta[[disease_col]])
  }
  obj$disease_group <- dplyr::case_when(
    grepl("azfa", obj$disease_group, ignore.case = TRUE) ~ "AZFa_Del",
    grepl("inoa|noa", obj$disease_group, ignore.case = TRUE) ~ "iNOA",
    TRUE ~ as.character(obj$disease_group)
  )
  obj
}

load_disease_objects <- function() {
  combined_file <- find_first_file(disease_candidates$combined[nzchar(disease_candidates$combined)])
  if (!is.na(combined_file)) {
    message("Loading combined disease object: ", combined_file)
    obj <- readRDS(combined_file)
    obj <- standardize_disease_metadata(obj, "Disease")
    return(obj)
  }
  
  objs <- list()
  for (nm in c("iNOA", "AZFa_Del")) {
    f <- find_first_file(disease_candidates[[nm]][nzchar(disease_candidates[[nm]])])
    if (!is.na(f)) {
      message("Loading ", nm, " object: ", f)
      obj <- readRDS(f)
      obj <- standardize_disease_metadata(obj, nm)
      objs[[nm]] <- obj
    }
  }
  if (length(objs) == 0) {
    stop("No disease Seurat object found. Set PART4_DISEASE_RDS or PART4_iNOA_RDS/PART4_AZFA_RDS.")
  }
  if (length(objs) == 1) return(objs[[1]])
  merge(objs[[1]], y = objs[-1], add.cell.ids = names(objs), project = "Disease_projection")
}

annotate_disease_celltypes_if_needed <- function(obj) {
  meta <- obj@meta.data
  
  if (nzchar(USER_CELLTYPE_COL) && USER_CELLTYPE_COL %in% colnames(meta) && !FORCE_MARKER_ANNOTATION) {
    message("Using user-specified disease cell-type column: ", USER_CELLTYPE_COL)
    obj$cell_type_part4 <- as.character(meta[[USER_CELLTYPE_COL]])
    return(list(
      obj = obj,
      score_qc = data.frame(),
      annotation_source = paste0("user_column:", USER_CELLTYPE_COL)
    ))
  }
  
  ct_col <- if (FORCE_MARKER_ANNOTATION) NA_character_ else first_existing_col(CELLTYPE_COL_CANDIDATES, meta)
  if (!is.na(ct_col)) {
    current <- as.character(meta[[ct_col]])
    recognized <- current %in% CELLTYPES_TO_PROJECT
    n_recognized_types <- length(unique(current[recognized]))
    enough_labels <- sum(recognized, na.rm = TRUE) >= MIN_DISEASE_CELLS
    enough_celltype_coverage <- n_recognized_types >= 3
    if (enough_labels && enough_celltype_coverage) {
      message("Using detected disease cell-type column: ", ct_col)
      obj$cell_type_part4 <- current
      return(list(
        obj = obj,
        score_qc = data.frame(),
        annotation_source = paste0("detected_column:", ct_col)
      ))
    }
  }
  
  message("No complete disease cell-type labels found. Running marker-score fallback annotation.")
  expr <- get_assay_data_safe(obj, ASSAY_USE, SLOT_USE)
  score_mat <- sapply(CELLTYPES_TO_PROJECT, function(ct) score_marker_set(expr, CELLTYPE_MARKERS[[ct]]))
  score_mat <- as.matrix(score_mat)
  rownames(score_mat) <- colnames(obj)
  
  z_mat <- apply(score_mat, 2, zscore)
  if (is.null(dim(z_mat))) z_mat <- matrix(z_mat, ncol = 1)
  colnames(z_mat) <- colnames(score_mat)
  rownames(z_mat) <- rownames(score_mat)
  
  top_idx <- apply(z_mat, 1, function(x) {
    if (all(!is.finite(x))) return(NA_integer_)
    which.max(ifelse(is.finite(x), x, -Inf))
  })
  top_score <- rep(NA_real_, nrow(z_mat))
  valid_top <- is.finite(top_idx)
  top_score[valid_top] <- z_mat[cbind(which(valid_top), top_idx[valid_top])]
  z_tmp <- z_mat
  z_tmp[cbind(which(valid_top), top_idx[valid_top])] <- NA_real_
  second_score <- apply(z_tmp, 1, function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) NA_real_ else max(x)
  })
  score_margin <- top_score - second_score
  label <- rep("Unknown", nrow(z_mat))
  label[valid_top] <- colnames(z_mat)[top_idx[valid_top]]
  label[!is.finite(top_score) | score_margin < 0.15] <- "Unknown"
  obj$cell_type_part4 <- label
  
  score_qc <- data.frame(
    cell = rownames(z_mat),
    sample_id = obj$sample_id,
    disease_group = obj$disease_group,
    assigned_cell_type = label,
    top_score = top_score,
    second_score = second_score,
    score_margin = score_margin,
    stringsAsFactors = FALSE
  )
  safe_fwrite(score_qc, file.path(TABLE_DIR, "disease_celltype_marker_score_qc.csv"))
  list(
    obj = obj,
    score_qc = score_qc,
    annotation_source = "marker_score_fallback"
  )
}

refs <- lapply(CELLTYPES_TO_PROJECT, load_reference)
names(refs) <- CELLTYPES_TO_PROJECT
refs <- refs[!vapply(refs, is.null, logical(1))]
if (length(refs) == 0) stop("No usable Part 3 references were found.")

disease_obj <- load_disease_objects()
anno <- annotate_disease_celltypes_if_needed(disease_obj)
disease_obj <- anno$obj
annotation_score_qc <- anno$score_qc
annotation_source <- anno$annotation_source

safe_fwrite(
  data.frame(
    annotation_source = annotation_source,
    user_celltype_col = USER_CELLTYPE_COL,
    force_marker_annotation = FORCE_MARKER_ANNOTATION,
    stringsAsFactors = FALSE
  ),
  file.path(TABLE_DIR, "part4_disease_celltype_annotation_source.csv")
)

disease_celltype_counts <- disease_obj@meta.data %>%
  mutate(
    sample_id = as.character(sample_id),
    disease_group = as.character(disease_group),
    cell_type_part4 = as.character(cell_type_part4)
  ) %>%
  count(disease_group, sample_id, cell_type_part4, name = "n_cells") %>%
  group_by(disease_group, sample_id) %>%
  mutate(fraction = n_cells / sum(n_cells)) %>%
  ungroup() %>%
  arrange(disease_group, sample_id, desc(n_cells))
safe_fwrite(disease_celltype_counts, file.path(TABLE_DIR, "disease_celltype_counts_for_projection.csv"))

safe_fwrite(
  data.frame(
    color_name = names(PAL),
    hex = unlist(PAL, use.names = FALSE),
    stringsAsFactors = FALSE
  ),
  file.path(TABLE_DIR, "part4_cell_style_palette.csv")
)

# ------------------------------------------------------------------------------
# 4. Projection core
# ------------------------------------------------------------------------------

choose_projection_genes <- function(ref, disease_ct, celltype) {
  ref_expr <- get_assay_data_safe(ref, ASSAY_USE, SLOT_USE)
  dis_expr <- get_assay_data_safe(disease_ct, ASSAY_USE, SLOT_USE)
  common <- intersect(rownames(ref_expr), rownames(dis_expr))
  vf <- tryCatch(VariableFeatures(ref), error = function(e) character(0))
  marker_genes <- unique(c(CELLTYPE_MARKERS[[celltype]], unlist(MATURE_MARKER_CONFIG[[celltype]])))
  marker_genes <- intersect(marker_genes, common)
  if (length(vf) == 0) {
    gene_vars <- Matrix::rowMeans(ref_expr[common, , drop = FALSE]^2) -
      Matrix::rowMeans(ref_expr[common, , drop = FALSE])^2
    vf <- names(sort(gene_vars, decreasing = TRUE))
  }
  genes <- unique(c(marker_genes, intersect(vf, common)))
  genes <- genes[genes %in% common]
  if (length(genes) > N_FEATURES_PROJECTION) {
    keep_marker <- intersect(marker_genes, genes)
    keep_vf <- setdiff(genes, keep_marker)
    genes <- unique(c(keep_marker, head(keep_vf, N_FEATURES_PROJECTION - length(keep_marker))))
  }
  genes
}

project_one_celltype <- function(celltype, ref, disease_obj) {
  disease_label <- as.character(disease_obj$cell_type_part4)
  disease_cells <- colnames(disease_obj)[which(!is.na(disease_label) & disease_label == celltype)]
  if (length(disease_cells) < MIN_DISEASE_CELLS) {
    warning("Skipping ", celltype, ": too few disease cells (", length(disease_cells), ").")
    return(NULL)
  }
  if (ncol(ref) < MIN_REF_CELLS) {
    warning("Skipping ", celltype, ": too few reference cells (", ncol(ref), ").")
    return(NULL)
  }
  
  message("Projecting cell type: ", celltype)
  disease_ct <- subset(disease_obj, cells = disease_cells)
  ref_cells_use <- select_balanced_cells(ref@meta.data, MAX_REF_CELLS_FOR_PCA, AGE_COL)
  ref_train <- subset(ref, cells = ref_cells_use)
  
  genes <- choose_projection_genes(ref_train, disease_ct, celltype)
  if (length(genes) < 100) {
    warning("Skipping ", celltype, ": too few shared projection genes (", length(genes), ").")
    return(NULL)
  }
  
  ref_expr <- get_assay_data_safe(ref_train, ASSAY_USE, SLOT_USE)
  dis_expr <- get_assay_data_safe(disease_ct, ASSAY_USE, SLOT_USE)
  scaled <- scale_rows_by_reference(ref_expr, dis_expr, genes)
  pca <- run_reference_pca(scaled$ref_z, scaled$query_z, N_PCS_PROJECTION)
  
  knn <- find_knn(pca$query_pca, pca$ref_pca, K_NEIGHBORS)
  ref_meta <- ref_train@meta.data
  ref_cells <- rownames(ref_meta)
  ref_pt <- to_num(ref_meta$pseudotime_final)
  ref_age <- to_num(ref_meta[[AGE_COL]])
  ref_maturation <- if ("maturation_score" %in% colnames(ref_meta)) {
    to_num(ref_meta$maturation_score)
  } else {
    rep(NA_real_, nrow(ref_meta))
  }
  
  ref_self <- find_knn(pca$ref_pca, pca$ref_pca, min(K_NEIGHBORS + 1, nrow(pca$ref_pca)))
  self_dist <- ref_self$dist[, -1, drop = FALSE]
  distance_threshold <- unname(quantile(rowMeans(self_dist, na.rm = TRUE),
                                        PROJECTION_DISTANCE_QC_QUANTILE,
                                        na.rm = TRUE))
  
  rows <- vector("list", nrow(pca$query_pca))
  neighbor_rows <- vector("list", nrow(pca$query_pca))
  for (i in seq_len(nrow(pca$query_pca))) {
    idx <- knn$idx[i, ]
    dist <- knn$dist[i, ]
    w <- weights_from_distance(dist)
    n_pt <- ref_pt[idx]
    n_age <- ref_age[idx]
    n_mat <- ref_maturation[idx]
    neighbor_cells <- ref_cells[idx]
    projected_pt <- weighted_mean_safe(n_pt, w)
    projected_age <- weighted_mean_safe(n_age, w)
    projected_mat <- weighted_mean_safe(n_mat, w)
    mean_dist <- mean(dist, na.rm = TRUE)
    quality <- case_when(
      is.finite(mean_dist) && mean_dist <= distance_threshold ~ "high",
      is.finite(mean_dist) && mean_dist <= distance_threshold * 1.25 ~ "moderate",
      TRUE ~ "low"
    )
    rows[[i]] <- data.frame(
      cell = rownames(pca$query_pca)[i],
      cell_type = celltype,
      disease_group = as.character(disease_ct$disease_group[i]),
      sample_id = as.character(disease_ct$sample_id[i]),
      projected_pseudotime = projected_pt,
      projected_maturation_score = projected_mat,
      weighted_reference_age = projected_age,
      nearest_reference_cell = neighbor_cells[[1]],
      nearest_reference_pseudotime = n_pt[[1]],
      nearest_reference_age = n_age[[1]],
      nearest_distance = dist[[1]],
      mean_knn_distance = mean_dist,
      kth_distance = dist[[length(dist)]],
      neighbor_age_sd = weighted_sd_safe(n_age, w),
      neighbor_pseudotime_sd = weighted_sd_safe(n_pt, w),
      weight_entropy = weight_entropy(w),
      projection_quality = quality,
      n_projection_genes = length(scaled$genes),
      n_reference_cells_used = ncol(scaled$ref_z),
      stringsAsFactors = FALSE
    )
    neighbor_rows[[i]] <- data.frame(
      cell = rownames(pca$query_pca)[i],
      cell_type = celltype,
      disease_group = as.character(disease_ct$disease_group[i]),
      sample_id = as.character(disease_ct$sample_id[i]),
      neighbor_rank = seq_along(idx),
      reference_cell = neighbor_cells,
      reference_age = n_age,
      reference_pseudotime = n_pt,
      distance = dist,
      weight = w,
      stringsAsFactors = FALSE
    )
  }
  
  projected <- bind_rows(rows)
  age_summary <- data.frame(
    age_numeric = ref_age,
    pseudotime_final = ref_pt,
    maturation_score = ref_maturation
  ) %>%
    filter(is.finite(age_numeric), is.finite(pseudotime_final)) %>%
    group_by(age_numeric) %>%
    summarise(
      n_reference_cells = n(),
      median_pseudotime = median(pseudotime_final, na.rm = TRUE),
      q25 = quantile(pseudotime_final, 0.25, na.rm = TRUE),
      q75 = quantile(pseudotime_final, 0.75, na.rm = TRUE),
      median_maturation_score = median(maturation_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(age_numeric)
  
  projected$inferred_normal_equiv_age <- infer_reference_age(projected$projected_pseudotime, age_summary)
  projected$nearest_normal_age <- nearest_reference_age(projected$projected_pseudotime, age_summary)
  mature_age <- max(age_summary$age_numeric, na.rm = TRUE)
  mature_ref <- age_summary %>% filter(age_numeric == mature_age) %>% dplyr::slice(1)
  projected$mature_reference_age <- mature_age
  projected$mature_reference_pseudotime <- mature_ref$median_pseudotime[[1]]
  projected$pseudotime_deficit_vs_mature <- mature_ref$median_pseudotime[[1]] - projected$projected_pseudotime
  projected$developmental_age_deficit_vs_mature <- mature_age - projected$inferred_normal_equiv_age
  
  reference_table <- data.frame(
    cell = ref_cells,
    cell_type = celltype,
    age_numeric = ref_age,
    pseudotime_final = ref_pt,
    maturation_score = ref_maturation,
    stringsAsFactors = FALSE
  )
  
  quality <- projected %>%
    group_by(cell_type, disease_group) %>%
    summarise(
      n_projected_cells = n(),
      n_samples = n_distinct(sample_id),
      median_mean_knn_distance = median(mean_knn_distance, na.rm = TRUE),
      high_quality_fraction = mean(projection_quality == "high", na.rm = TRUE),
      low_quality_fraction = mean(projection_quality == "low", na.rm = TRUE),
      median_neighbor_age_sd = median(neighbor_age_sd, na.rm = TRUE),
      n_projection_genes = median(n_projection_genes, na.rm = TRUE),
      distance_threshold = distance_threshold,
      .groups = "drop"
    )
  
  list(
    projected = projected,
    neighbors = bind_rows(neighbor_rows),
    reference = reference_table,
    age_summary = age_summary %>% mutate(cell_type = celltype),
    quality = quality
  )
}

projection_results <- lapply(names(refs), function(ct) {
  project_one_celltype(ct, refs[[ct]], disease_obj)
})
names(projection_results) <- names(refs)
projection_results <- projection_results[!vapply(projection_results, is.null, logical(1))]
if (length(projection_results) == 0) stop("No cell type was projected successfully.")

projected_all <- bind_rows(lapply(projection_results, `[[`, "projected"))
neighbors_all <- bind_rows(lapply(projection_results, `[[`, "neighbors"))
reference_all <- bind_rows(lapply(projection_results, `[[`, "reference"))
reference_age_all <- bind_rows(lapply(projection_results, `[[`, "age_summary"))
quality_all <- bind_rows(lapply(projection_results, `[[`, "quality"))

# Attach projected values back to the disease object for downstream Part 5/6/7.
for (cc in c(
  "cell_type", "projected_pseudotime", "projected_maturation_score",
  "weighted_reference_age", "nearest_normal_age", "inferred_normal_equiv_age",
  "pseudotime_deficit_vs_mature", "developmental_age_deficit_vs_mature",
  "projection_quality", "mean_knn_distance", "neighbor_age_sd"
)) {
  disease_obj@meta.data[[paste0("part4_", cc)]] <- NA
}

match_idx <- match(projected_all$cell, rownames(disease_obj@meta.data))
ok <- !is.na(match_idx)
for (cc in c(
  "cell_type", "projected_pseudotime", "projected_maturation_score",
  "weighted_reference_age", "nearest_normal_age", "inferred_normal_equiv_age",
  "pseudotime_deficit_vs_mature", "developmental_age_deficit_vs_mature",
  "projection_quality", "mean_knn_distance", "neighbor_age_sd"
)) {
  disease_obj@meta.data[match_idx[ok], paste0("part4_", cc)] <- projected_all[[cc]][ok]
}

safe_fwrite(projected_all, file.path(TABLE_DIR, "disease_projected_pseudotime_values.csv"))
safe_fwrite(neighbors_all, file.path(TABLE_DIR, "disease_projection_knn_neighbors.csv"))
safe_fwrite(reference_all, file.path(TABLE_DIR, "normal_reference_cell_level_values_for_projection.csv"))
safe_fwrite(reference_age_all, file.path(TABLE_DIR, "normal_reference_age_summary_for_projection.csv"))
safe_fwrite(quality_all, file.path(TABLE_DIR, "disease_projection_quality_by_celltype.csv"))
saveRDS(disease_obj, file.path(RDS_DIR, "combined_final_with_projections.rds"))

# ------------------------------------------------------------------------------
# 5. Sample-level and group-level summaries
# ------------------------------------------------------------------------------

summarise_projection_by_sample <- function(projected, reference_age) {
  mature_refs <- reference_age %>%
    group_by(cell_type) %>%
    filter(age_numeric == max(age_numeric, na.rm = TRUE)) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    dplyr::select(cell_type,
                  ref_mature_age = age_numeric,
                  ref_mature_pseudotime = median_pseudotime,
                  ref_mature_maturation_score = median_maturation_score)
  
  projected %>%
    left_join(mature_refs, by = "cell_type") %>%
    group_by(cell_type, disease_group, sample_id) %>%
    summarise(
      n_cells = n(),
      n_high_quality_cells = sum(projection_quality == "high", na.rm = TRUE),
      high_quality_fraction = mean(projection_quality == "high", na.rm = TRUE),
      median_projected_pseudotime = median(projected_pseudotime, na.rm = TRUE),
      mean_projected_pseudotime = mean(projected_pseudotime, na.rm = TRUE),
      q25_projected_pseudotime = quantile(projected_pseudotime, 0.25, na.rm = TRUE),
      q75_projected_pseudotime = quantile(projected_pseudotime, 0.75, na.rm = TRUE),
      median_projected_maturation_score = median(projected_maturation_score, na.rm = TRUE),
      median_weighted_reference_age = median(weighted_reference_age, na.rm = TRUE),
      median_inferred_normal_equiv_age = median(inferred_normal_equiv_age, na.rm = TRUE),
      median_nearest_normal_age = median(nearest_normal_age, na.rm = TRUE),
      median_mean_knn_distance = median(mean_knn_distance, na.rm = TRUE),
      median_neighbor_age_sd = median(neighbor_age_sd, na.rm = TRUE),
      mature_age = dplyr::first(ref_mature_age),
      mature_reference_pseudotime = dplyr::first(ref_mature_pseudotime),
      mature_reference_maturation_score = dplyr::first(ref_mature_maturation_score),
      median_pseudotime_deficit_vs_mature =
        dplyr::first(ref_mature_pseudotime) - median(projected_pseudotime, na.rm = TRUE),
      median_maturation_deficit_vs_mature =
        dplyr::first(ref_mature_maturation_score) - median(projected_maturation_score, na.rm = TRUE),
      developmental_age_deficit_vs_mature =
        dplyr::first(ref_mature_age) - median(inferred_normal_equiv_age, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(cell_type) %>%
    mutate(
      immature_reference_pseudotime = min(reference_age$median_pseudotime[reference_age$cell_type == dplyr::first(cell_type)], na.rm = TRUE),
      normalized_maturity_position =
        (median_projected_pseudotime - immature_reference_pseudotime) /
        (mature_reference_pseudotime - immature_reference_pseudotime)
    ) %>%
    ungroup()
}

sample_summary <- summarise_projection_by_sample(projected_all, reference_age_all)

sample_wasserstein <- projected_all %>%
  group_by(cell_type, disease_group, sample_id) %>%
  summarise(
    wasserstein_to_mature_reference = {
      ct <- dplyr::first(cell_type)
      mature_age <- max(reference_age_all$age_numeric[reference_age_all$cell_type == ct], na.rm = TRUE)
      mature_ref_cells <- reference_all %>%
        filter(cell_type == ct, age_numeric == mature_age) %>%
        pull(pseudotime_final)
      wasserstein_1d(projected_pseudotime, mature_ref_cells)
    },
    .groups = "drop"
  )

sample_summary <- sample_summary %>%
  left_join(sample_wasserstein, by = c("cell_type", "disease_group", "sample_id")) %>%
  arrange(cell_type, disease_group, sample_id)

group_summary <- sample_summary %>%
  group_by(cell_type, disease_group) %>%
  summarise(
    n_samples = n_distinct(sample_id),
    high_quality_fraction = {
      ok <- is.finite(high_quality_fraction) & is.finite(n_cells) & n_cells > 0
      if (!any(ok)) NA_real_ else sum(high_quality_fraction[ok] * n_cells[ok]) / sum(n_cells[ok])
    },
    n_cells = sum(n_cells, na.rm = TRUE),
    median_projected_pseudotime = median(median_projected_pseudotime, na.rm = TRUE),
    median_projected_maturation_score = median(median_projected_maturation_score, na.rm = TRUE),
    median_inferred_normal_equiv_age = median(median_inferred_normal_equiv_age, na.rm = TRUE),
    median_pseudotime_deficit_vs_mature = median(median_pseudotime_deficit_vs_mature, na.rm = TRUE),
    median_maturation_deficit_vs_mature = median(median_maturation_deficit_vs_mature, na.rm = TRUE),
    median_developmental_age_deficit_vs_mature = median(developmental_age_deficit_vs_mature, na.rm = TRUE),
    mean_wasserstein_to_mature_reference = mean(wasserstein_to_mature_reference, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(cell_type, disease_group)

priority_table <- sample_summary %>%
  filter(disease_group == "iNOA") %>%
  group_by(cell_type) %>%
  summarise(
    n_inoa_samples = n_distinct(sample_id),
    n_inoa_cells = sum(n_cells, na.rm = TRUE),
    median_inoa_pseudotime_deficit = median(median_pseudotime_deficit_vs_mature, na.rm = TRUE),
    median_inoa_developmental_age_deficit =
      median(developmental_age_deficit_vs_mature, na.rm = TRUE),
    fraction_inoa_samples_positive =
      mean(median_pseudotime_deficit_vs_mature > 0, na.rm = TRUE),
    median_high_quality_fraction = median(high_quality_fraction, na.rm = TRUE),
    priority_score =
      median_inoa_pseudotime_deficit * fraction_inoa_samples_positive * median_high_quality_fraction,
    .groups = "drop"
  ) %>%
  arrange(desc(priority_score), desc(median_inoa_pseudotime_deficit))

safe_fwrite(sample_summary, file.path(TABLE_DIR, "disease_projection_summary_by_sample.csv"))
safe_fwrite(group_summary, file.path(TABLE_DIR, "disease_projection_summary_by_group.csv"))
safe_fwrite(priority_table, file.path(TABLE_DIR, "part4_multicelltype_iNOA_developmental_deficit_priority.csv"))

# ------------------------------------------------------------------------------
# 6. Figures
# ------------------------------------------------------------------------------

reference_axis <- reference_age_all %>%
  mutate(
    age_numeric = to_num(age_numeric),
    cell_type = factor(cell_type, levels = CELLTYPES_TO_PROJECT)
  )

sample_plot <- sample_summary %>%
  mutate(
    cell_type = factor(cell_type, levels = CELLTYPES_TO_PROJECT),
    disease_group = factor(disease_group, levels = c("AZFa_Del", "iNOA", "Disease", "Unknown"))
  )

projected_plot <- projected_all
if (nrow(projected_plot) > MAX_DISEASE_CELLS_FOR_PLOTS) {
  cells_per_plot_group <- ceiling(
    MAX_DISEASE_CELLS_FOR_PLOTS /
      max(1, dplyr::n_distinct(projected_all$cell_type, projected_all$disease_group, projected_all$sample_id))
  )
  projected_plot <- projected_plot %>%
    group_by(cell_type, disease_group, sample_id) %>%
    dplyr::group_modify(~ dplyr::slice_sample(.x, n = min(nrow(.x), cells_per_plot_group))) %>%
    ungroup()
}

# Fig5A: multi-cell-type priority screen.
if (nrow(priority_table) > 0) {
  p_priority <- priority_table %>%
    mutate(
      cell_type = factor(cell_type, levels = rev(cell_type)),
      is_sertoli = as.character(cell_type) == "Sertoli"
    ) %>%
    ggplot(aes(x = cell_type, y = median_inoa_pseudotime_deficit, fill = is_sertoli)) +
    geom_col(width = 0.72, color = PAL$black, linewidth = 0.25) +
    geom_point(aes(y = priority_score), shape = 21, size = 2.7,
               fill = PAL$white, color = PAL$black, stroke = 0.35) +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = PAL$pink_dark, "FALSE" = PAL$blue_light), guide = "none") +
    labs(
      x = NULL,
      y = "Median iNOA pseudotime deficit vs mature reference",
      title = "Multi-cell-type projection reveals broad developmental shifts in iNOA",
      subtitle = "Bars show sample-median deficits; open dots show quality-weighted scores; Sertoli is highlighted for support-function follow-up"
    ) +
    theme_cell(9)
  save_fig("Fig5A_multicelltype_iNOA_developmental_deficit_priority.pdf",
           p_priority, width = 6.4, height = 4.4, section = "main")
}

# Fig5B: Sertoli projection onto normal age axis.
sertoli_ref_axis <- reference_axis %>% filter(cell_type == "Sertoli")
sertoli_sample <- sample_plot %>% filter(cell_type == "Sertoli")
if (nrow(sertoli_ref_axis) > 0 && nrow(sertoli_sample) > 0) {
  p_sertoli_axis <- ggplot(sertoli_ref_axis, aes(x = age_numeric, y = median_pseudotime)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = PAL$blue_light, alpha = 0.28) +
    geom_line(color = PAL$blue_dark, linewidth = 0.85) +
    geom_point(color = PAL$blue_dark, fill = PAL$white, shape = 21, size = 2.4, stroke = 0.35) +
    geom_point(
      data = sertoli_sample,
      aes(x = median_inferred_normal_equiv_age, y = median_projected_pseudotime,
          fill = disease_group, shape = disease_group),
      size = 3.2, color = PAL$black, stroke = 0.35,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = disease_palette, na.value = PAL$grey_dark) +
    scale_shape_manual(values = c("AZFa_Del" = 22, "iNOA" = 24, "Disease" = 21, "Unknown" = 21)) +
    labs(
      x = "Normal-equivalent developmental age",
      y = "Sertoli developmental pseudotime",
      fill = "Disease",
      shape = "Disease",
      title = "Disease Sertoli cells project to earlier normal developmental states",
      subtitle = "Blue line shows normal Sertoli age medians; disease symbols show sample medians"
    ) +
    theme_cell(9)
  if (HAS_GGREPEL) {
    p_sertoli_axis <- p_sertoli_axis +
      ggrepel::geom_text_repel(
        data = sertoli_sample,
        aes(x = median_inferred_normal_equiv_age, y = median_projected_pseudotime, label = sample_id),
        size = 2.6, color = PAL$black, min.segment.length = 0,
        box.padding = 0.25, point.padding = 0.18,
        inherit.aes = FALSE
      )
  }
  save_fig("Fig5B_Sertoli_projection_on_normal_axis.pdf",
           p_sertoli_axis, width = 5.6, height = 4.4, section = "main")
}

# Fig5C: Sertoli sample-level deficit.
if (nrow(sertoli_sample) > 0) {
  p_sertoli_deficit <- sertoli_sample %>%
    arrange(disease_group, desc(median_pseudotime_deficit_vs_mature)) %>%
    mutate(sample_id = factor(sample_id, levels = unique(sample_id))) %>%
    ggplot(aes(x = sample_id, y = median_pseudotime_deficit_vs_mature, fill = disease_group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = PAL$grey_dark, linewidth = 0.35) +
    geom_col(width = 0.68, color = PAL$black, linewidth = 0.25) +
    geom_text(aes(label = sprintf("%.2f", median_pseudotime_deficit_vs_mature)),
              vjust = -0.35, size = 2.7, color = PAL$black) +
    scale_fill_manual(values = disease_palette, na.value = PAL$grey_dark) +
    labs(
      x = NULL,
      y = "Pseudotime deficit vs mature normal Sertoli",
      fill = "Disease",
      title = "iNOA Sertoli cells show a consistent maturation deficit",
      subtitle = "Positive values indicate projection earlier than mature normal Sertoli cells"
    ) +
    theme_cell(9) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  save_fig("Fig5C_Sertoli_pseudotime_deficit_by_sample.pdf",
           p_sertoli_deficit, width = 5.8, height = 4.1, section = "main")
}

# Fig5D: projection workflow schematic for the main text.
workflow_nodes <- data.frame(
  step = c("Normal reference", "Disease mapping", "Sample-level readout", "Interpretation"),
  label = c(
    "Normal cell-type\nreference trajectories",
    "KNN projection in\nreference PCA space",
    "Projected pseudotime\nand maturity deficit",
    "Early-state shift in\niNOA Sertoli cells"
  ),
  x = c(1, 2.45, 3.9, 5.35),
  y = c(1, 1, 1, 1),
  fill = c(PAL$blue_light, PAL$lavender, PAL$pink_light, PAL$pink_mid),
  stringsAsFactors = FALSE
)
workflow_edges <- data.frame(
  x = workflow_nodes$x[-nrow(workflow_nodes)] + 0.48,
  xend = workflow_nodes$x[-1] - 0.48,
  y = 1,
  yend = 1
)
p_workflow <- ggplot() +
  geom_segment(
    data = workflow_edges,
    aes(x = x, xend = xend, y = y, yend = yend),
    arrow = grid::arrow(length = grid::unit(2.2, "mm"), type = "closed"),
    linewidth = 0.35,
    color = PAL$grey_dark
  ) +
  geom_label(
    data = workflow_nodes,
    aes(x = x, y = y, label = label, fill = step),
    label.size = 0.28,
    label.r = grid::unit(1.8, "mm"),
    size = 3.2,
    lineheight = 0.95,
    color = PAL$black,
    label.padding = grid::unit(3.6, "mm")
  ) +
  scale_fill_manual(values = setNames(workflow_nodes$fill, workflow_nodes$step), guide = "none") +
  coord_cartesian(xlim = c(0.4, 5.95), ylim = c(0.58, 1.42), expand = FALSE) +
  labs(
    title = "Disease-cell projection strategy",
    subtitle = "Disease cells are interpreted by their position on normal cell-type-specific developmental references"
  ) +
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold", size = 11, hjust = 0, color = PAL$black),
    plot.subtitle = element_text(size = 8.8, hjust = 0, color = PAL$grey_dark),
    plot.margin = margin(5, 6, 5, 5),
    plot.background = element_rect(fill = PAL$white, color = NA)
  )
save_fig("Fig5D_projection_workflow_schematic.pdf",
         p_workflow, width = 7.2, height = 1.7, section = "main")

# Supplementary Sertoli projection quality / maturation coupling QC.
if (nrow(sertoli_sample) > 0) {
  sertoli_sample_q <- sertoli_sample %>%
    mutate(high_quality_fraction_plot = pmin(pmax(high_quality_fraction, 0), 1))
  has_maturation_axis <- sum(is.finite(sertoli_ref_axis$median_maturation_score)) >= 2 &&
    sum(is.finite(sertoli_sample_q$median_projected_maturation_score)) >= 2
  
  if (has_maturation_axis) {
    p_quality_maturation <- ggplot(
      sertoli_sample_q,
      aes(x = median_projected_pseudotime, y = median_projected_maturation_score,
          fill = disease_group, size = high_quality_fraction_plot)
    ) +
      geom_path(
        data = sertoli_ref_axis %>% filter(is.finite(median_maturation_score)),
        aes(x = median_pseudotime, y = median_maturation_score),
        inherit.aes = FALSE,
        color = PAL$blue_dark, linewidth = 0.8, alpha = 0.85
      ) +
      geom_point(
        data = sertoli_ref_axis %>% filter(is.finite(median_maturation_score)),
        aes(x = median_pseudotime, y = median_maturation_score),
        inherit.aes = FALSE,
        shape = 21, fill = PAL$white, color = PAL$blue_dark, size = 2.3, stroke = 0.35
      ) +
      geom_point(shape = 21, color = PAL$black, stroke = 0.35, alpha = 0.95) +
      labs(
        x = "Median projected pseudotime",
        y = "Median projected maturation score",
        title = "Projected pseudotime remains coupled to maturation programs",
        subtitle = "Reference age medians are shown as the blue trajectory"
      )
  } else {
    p_quality_maturation <- ggplot(
      sertoli_sample_q,
      aes(x = median_projected_pseudotime, y = high_quality_fraction_plot,
          fill = disease_group, size = n_cells)
    ) +
      geom_hline(yintercept = c(0.5, 0.8), linetype = "dotted",
                 color = PAL$grey_dark, linewidth = 0.3) +
      geom_point(shape = 21, color = PAL$black, stroke = 0.35, alpha = 0.95) +
      scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
      labs(
        x = "Median projected pseudotime",
        y = "High-quality projection fraction",
        title = "Sertoli projection quality supports the developmental shift",
        subtitle = "Fallback view used because maturation-score values are incomplete"
      )
  }
  
  p_quality_maturation <- p_quality_maturation +
    scale_fill_manual(values = disease_palette, na.value = PAL$grey_dark) +
    scale_size_continuous(range = c(2.4, 5.2), name = ifelse(has_maturation_axis, "High-quality\nfraction", "Cells")) +
    labs(fill = "Disease") +
    theme_cell(9)
  if (HAS_GGREPEL) {
    p_quality_maturation <- p_quality_maturation +
      ggrepel::geom_text_repel(
        aes(label = sample_id),
        size = 2.6, color = PAL$black, min.segment.length = 0,
        box.padding = 0.25, point.padding = 0.18
      )
  }
  save_fig("FigS5G_Sertoli_projection_quality_qc.pdf",
           p_quality_maturation, width = 5.8, height = 4.4, section = "supplementary")
}

# FigS5A: all cell types on their normal axes.
p_all_axes <- ggplot(reference_axis, aes(x = age_numeric, y = median_pseudotime)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = PAL$blue_light, alpha = 0.22) +
  geom_line(color = PAL$blue_dark, linewidth = 0.65) +
  geom_point(color = PAL$blue_dark, fill = PAL$white, shape = 21, size = 1.7, stroke = 0.25) +
  geom_point(
    data = sample_plot,
    aes(x = median_inferred_normal_equiv_age, y = median_projected_pseudotime,
        fill = disease_group, shape = disease_group),
    color = PAL$black, size = 2.1, stroke = 0.28,
    inherit.aes = FALSE
  ) +
  facet_wrap(~cell_type, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = disease_palette, na.value = PAL$grey_dark) +
  scale_shape_manual(values = c("AZFa_Del" = 22, "iNOA" = 24, "Disease" = 21, "Unknown" = 21)) +
  labs(
    x = "Normal-equivalent developmental age",
    y = "Developmental pseudotime",
    fill = "Disease",
    shape = "Disease",
    title = "Disease-cell projections across normal cell-type developmental references"
  ) +
  theme_cell(8.4)
save_fig("FigS5A_all_celltypes_projection_axis.pdf",
         p_all_axes, width = 9.8, height = 7.2, section = "supplementary")

# FigS5B: cell-level density.
density_df <- bind_rows(
  reference_all %>%
    group_by(cell_type) %>%
    mutate(ref_age_label = paste0("Normal age ", age_numeric)) %>%
    ungroup() %>%
    transmute(cell_type, group = ref_age_label, disease_group = "Normal",
              pseudotime = pseudotime_final),
  projected_plot %>%
    transmute(cell_type, group = disease_group, disease_group,
              pseudotime = projected_pseudotime)
) %>%
  filter(is.finite(pseudotime))

p_density <- density_df %>%
  mutate(
    group_class = ifelse(disease_group == "Normal", "Normal reference", as.character(disease_group)),
    group_class = factor(group_class, levels = c("Normal reference", "AZFa_Del", "iNOA", "Disease", "Unknown"))
  ) %>%
  ggplot(aes(x = pseudotime, fill = group_class, color = group_class)) +
  geom_density(alpha = 0.22, linewidth = 0.35, adjust = 1.1) +
  facet_wrap(~cell_type, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = c("Normal reference" = PAL$blue_light, disease_palette), na.value = PAL$grey_dark) +
  scale_color_manual(values = c("Normal reference" = PAL$blue_dark, disease_palette), na.value = PAL$grey_dark) +
  labs(
    x = "Projected / reference pseudotime",
    y = "Density",
    fill = NULL,
    color = NULL,
    title = "Cell-level projected pseudotime distributions",
    subtitle = "Use as distributional support; sample-level summaries remain the primary statistical unit"
  ) +
  theme_cell(8.4)
save_fig("FigS5B_cell_level_projected_pseudotime_density.pdf",
         p_density, width = 10.4, height = 7.4, section = "supplementary")

# FigS5C: quality by cell type.
p_quality <- quality_all %>%
  mutate(cell_type = factor(cell_type, levels = CELLTYPES_TO_PROJECT)) %>%
  ggplot(aes(x = cell_type, y = high_quality_fraction, fill = disease_group)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62,
           color = PAL$black, linewidth = 0.22) +
  geom_hline(yintercept = c(0.5, 0.8), linetype = "dotted", color = PAL$grey_dark, linewidth = 0.28) +
  scale_fill_manual(values = disease_palette, na.value = PAL$grey_dark) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    x = NULL,
    y = "High-quality projection fraction",
    fill = "Disease",
    title = "Projection quality across cell types"
  ) +
  theme_cell(8.6) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_fig("FigS5C_projection_quality_by_celltype.pdf",
         p_quality, width = 8.4, height = 4.6, section = "supplementary")

# FigS5D: nearest reference age composition.
age_comp <- projected_all %>%
  group_by(cell_type, disease_group, sample_id, nearest_normal_age) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cell_type, disease_group, sample_id) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

age_levels <- sort(unique(age_comp$nearest_normal_age))
age_cols <- setNames(colorRampPalette(c(PAL$blue_light, PAL$pink_dark))(length(age_levels)), as.character(age_levels))

p_age_comp <- age_comp %>%
  mutate(
    cell_type = factor(cell_type, levels = CELLTYPES_TO_PROJECT),
    sample_id = factor(sample_id, levels = unique(sample_id)),
    nearest_normal_age = factor(nearest_normal_age, levels = age_levels)
  ) %>%
  ggplot(aes(x = sample_id, y = prop, fill = nearest_normal_age)) +
  geom_col(width = 0.72, color = PAL$white, linewidth = 0.12) +
  facet_wrap(~cell_type, ncol = 3, scales = "free_x") +
  scale_fill_manual(values = age_cols, name = "Nearest\nnormal age") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = NULL,
    y = "Cell fraction",
    title = "Nearest normal developmental age composition"
  ) +
  theme_cell(8.2) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
save_fig("FigS5D_nearest_reference_age_composition.pdf",
         p_age_comp, width = 10.4, height = 7.4, section = "supplementary")

# FigS5E: disease cell-type composition used for projection.
composition_plot <- disease_celltype_counts %>%
  filter(cell_type_part4 != "Unknown") %>%
  mutate(
    sample_id = factor(sample_id, levels = unique(sample_id)),
    cell_type_part4 = factor(cell_type_part4, levels = CELLTYPES_TO_PROJECT)
  ) %>%
  ggplot(aes(x = sample_id, y = fraction, fill = cell_type_part4)) +
  geom_col(width = 0.72, color = PAL$white, linewidth = 0.12) +
  scale_fill_manual(values = celltype_palette, na.value = PAL$grey_mid, name = "Cell type") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = NULL,
    y = "Cell fraction",
    title = "Disease cell-type composition used for projection",
    subtitle = paste0("Annotation source: ", annotation_source)
  ) +
  theme_cell(8.5) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_fig("FigS5E_disease_celltype_composition_for_projection.pdf",
         composition_plot, width = 7.6, height = 4.6, section = "supplementary")

# FigS5F: marker-score fallback QC, only when fallback annotation was used.
if (nrow(annotation_score_qc) > 0) {
  p_anno <- annotation_score_qc %>%
    filter(assigned_cell_type != "Unknown") %>%
    ggplot(aes(x = assigned_cell_type, y = score_margin, fill = assigned_cell_type)) +
    geom_violin(scale = "width", alpha = 0.72, color = NA) +
    geom_boxplot(width = 0.16, outlier.shape = NA, color = PAL$black, linewidth = 0.28) +
    scale_fill_manual(values = celltype_palette, guide = "none") +
    labs(
      x = NULL,
      y = "Top marker-score margin",
      title = "Fallback disease cell-type annotation score margin",
      subtitle = "Higher margins indicate cleaner marker-based assignment"
    ) +
    theme_cell(8.5) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  save_fig("FigS5F_celltype_annotation_score_qc.pdf",
           p_anno, width = 8.6, height = 4.6, section = "supplementary")
}

# ------------------------------------------------------------------------------
# 7. Manifest and report
# ------------------------------------------------------------------------------

figure_manifest <- data.frame(
  figure = c(
    "Fig5A_multicelltype_iNOA_developmental_deficit_priority.pdf",
    "Fig5B_Sertoli_projection_on_normal_axis.pdf",
    "Fig5C_Sertoli_pseudotime_deficit_by_sample.pdf",
    "Fig5D_Sertoli_projection_quality_and_maturation_coupling.pdf",
    "FigS5A_all_celltypes_projection_axis.pdf",
    "FigS5B_cell_level_projected_pseudotime_density.pdf",
    "FigS5C_projection_quality_by_celltype.pdf",
    "FigS5D_nearest_reference_age_composition.pdf",
    "FigS5E_celltype_annotation_score_qc.pdf"
  ),
  placement = c(
    rep("Main text", 4),
    rep("Supplementary", 5)
  ),
  file = c(
    file.path("figures/main", "Fig5A_multicelltype_iNOA_developmental_deficit_priority.pdf"),
    file.path("figures/main", "Fig5B_Sertoli_projection_on_normal_axis.pdf"),
    file.path("figures/main", "Fig5C_Sertoli_pseudotime_deficit_by_sample.pdf"),
    file.path("figures/main", "Fig5D_Sertoli_projection_quality_and_maturation_coupling.pdf"),
    file.path("figures/supplementary", "FigS5A_all_celltypes_projection_axis.pdf"),
    file.path("figures/supplementary", "FigS5B_cell_level_projected_pseudotime_density.pdf"),
    file.path("figures/supplementary", "FigS5C_projection_quality_by_celltype.pdf"),
    file.path("figures/supplementary", "FigS5D_nearest_reference_age_composition.pdf"),
    file.path("figures/supplementary", "FigS5E_celltype_annotation_score_qc.pdf")
  ),
  purpose = c(
    "Prioritize which cell types show the strongest iNOA developmental deficit.",
    "Place disease Sertoli sample medians on the normal Sertoli developmental axis.",
    "Show sample-level Sertoli pseudotime deficit relative to mature normal reference.",
    "Show maturation-score coupling and projection quality for Sertoli samples.",
    "Display all projected cell types on their normal developmental axes.",
    "Show cell-level projected pseudotime distribution support.",
    "Document projection quality across cell types and diseases.",
    "Show nearest normal age composition for projected disease cells.",
    "QC for fallback marker-score cell-type annotation when used."
  ),
  stringsAsFactors = FALSE
)
safe_fwrite(figure_manifest, file.path(TABLE_DIR, "part4_figure_manifest.csv"))

report_lines <- c(
  "Part 4 optimized disease projection report",
  "==========================================",
  paste0("Base directory: ", BASE_DIR),
  paste0("Projected cell types: ", paste(unique(projected_all$cell_type), collapse = ", ")),
  paste0("Projected disease cells: ", nrow(projected_all)),
  "",
  "Primary interpretation:",
  "- Use sample-level projected pseudotime and maturation deficit as the main evidence.",
  "- Use cell-level distributions only as supplementary support because cells from the same patient are not independent samples.",
  "- A positive pseudotime deficit means the disease sample projects earlier than mature normal cells of the same cell type.",
  "",
  "Main text figures:",
  paste0("- ", figure_manifest$figure[figure_manifest$placement == "Main text"]),
  "",
  "Supplementary figures:",
  paste0("- ", figure_manifest$figure[figure_manifest$placement == "Supplementary"]),
  "",
  "Downstream compatibility:",
  "- tables/disease_projected_pseudotime_values.csv is compatible with Part 5.",
  "- tables/disease_projection_summary_by_sample.csv and by_group.csv are compatible with Part 6/7.",
  "- rds/combined_final_with_projections.rds stores part4_* metadata for downstream mechanism exploration."
)
writeLines(report_lines, file.path(REPORT_DIR, "part4_projection_report.txt"))

message("Part 4 optimized disease projection finished: ", normalizePath(OUTDIR, winslash = "/"))

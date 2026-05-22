#!/usr/bin/env Rscript

# ==============================================================================
# 03_developmental_reference_construction.R
#
# Purpose:
#   Build normal postnatal developmental references for Sertoli and other major
#   human testicular cell types using PAGA/DPT pseudotime, age ordering,
#   marker-based maturation refinement, and cross-cell-type trajectory-quality
#   summaries.
#
# Scope:
#   This script constructs normal developmental reference axes only. Disease
#   cells are not used to define these trajectories; disease-cell projection is
#   handled in the next workflow step.
#
# Expected inputs:
#   <NORMAL_CELLTYPE_RDS_DIR>/sertoli_final_annotated.rds
#   <NORMAL_CELLTYPE_RDS_DIR>/Leydig_final_annotated.rds
#   <NORMAL_CELLTYPE_RDS_DIR>/PMC_final_annotated.rds
#   ...
#
# Main outputs:
#   results/03_developmental_reference_construction/tables/
#   results/03_developmental_reference_construction/figures/
#   results/03_developmental_reference_construction/rds/
#
# Default paths and Python settings can be overridden with environment variables:
#   PROJECT_DIR, NORMAL_CELLTYPE_RDS_DIR, OUTPUT_DIR, RETICULATE_PYTHON,
#   PART3_CONDA_ENV
# ==============================================================================

# IMPORTANT for Windows/RStudio/reticulate:
# Set the Python environment before loading reticulate or packages that may
# initialize Python. The environment should contain scanpy. Do not hard-code a
# personal Python path in the repository; set RETICULATE_PYTHON or PART3_CONDA_ENV
# before running this script when needed.
if (nzchar(Sys.getenv("RETICULATE_PYTHON", unset = ""))) {
  Sys.setenv(RETICULATE_PYTHON = Sys.getenv("RETICULATE_PYTHON"))
}
Sys.setenv(NUMBA_DISABLE_INTEL_SVML = Sys.getenv("NUMBA_DISABLE_INTEL_SVML", unset = "1"))
Sys.setenv(NUMBA_THREADING_LAYER = Sys.getenv("NUMBA_THREADING_LAYER", unset = "workqueue"))

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(readr)
  library(reticulate)
  library(ggplot2)
  library(ggridges)
  library(patchwork)
  library(scales)
})

# ------------------------------------------------------------------------------
# 0. User configuration
# ------------------------------------------------------------------------------

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = getwd())
INPUT_DIR <- Sys.getenv(
  "NORMAL_CELLTYPE_RDS_DIR",
  unset = file.path(PROJECT_DIR, "data", "processed", "normal_celltype_objects")
)
OUTDIR <- Sys.getenv(
  "OUTPUT_DIR",
  unset = file.path(PROJECT_DIR, "results", "03_developmental_reference_construction")
)
TABLE_DIR <- file.path(OUTDIR, "tables")
RDS_DIR <- file.path(OUTDIR, "rds")
FIG_DIR <- file.path(OUTDIR, "figures")
MAIN_FIG_DIR <- file.path(FIG_DIR, "main")
SUPP_FIG_DIR <- file.path(FIG_DIR, "supplementary")
PAGA_FIG_DIR <- file.path(SUPP_FIG_DIR, "paga_networks")
MARKER_TREND_FIG_DIR <- file.path(SUPP_FIG_DIR, "marker_trends")

dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RDS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(MAIN_FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(SUPP_FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PAGA_FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(MARKER_TREND_FIG_DIR, showWarnings = FALSE, recursive = TRUE)

CELLTYPES_TO_RUN <- c(
  "Sertoli", "Leydig", "PMC", "Endothelial", "Macrophage",
  "T_cell", "Spermatogonia", "Spermatocyte", "Spermatid"
)

INPUT_MAP <- c(
  Sertoli       = file.path(INPUT_DIR, "sertoli_final_annotated.rds"),
  Leydig        = file.path(INPUT_DIR, "Leydig_final_annotated.rds"),
  PMC           = file.path(INPUT_DIR, "PMC_final_annotated.rds"),
  Endothelial   = file.path(INPUT_DIR, "Endothelial_final_annotated.rds"),
  Macrophage    = file.path(INPUT_DIR, "Macrophage_final_annotated.rds"),
  T_cell        = file.path(INPUT_DIR, "T_cell_final_annotated.rds"),
  Spermatogonia = file.path(INPUT_DIR, "Spermatogonia_final_annotated.rds"),
  Spermatocyte  = file.path(INPUT_DIR, "Spermatocyte_final_annotated.rds"),
  Spermatid     = file.path(INPUT_DIR, "Spermatid_final_annotated.rds")
)

AGE_COL <- "age_numeric"
DONOR_COL <- "sample"
CLUSTER_COL <- "seurat_clusters"
ASSAY_USE <- "RNA"
SLOT_USE <- "data"

# PAGA/DPT.
N_PCS <- 30
N_NEIGHBORS_PAGA <- 15
RESOLUTION_PAGA <- 0.5
N_DCS <- 15
PAGA_EDGE_MIN_WEIGHT <- 0.01
PAGA_MAX_EDGES <- 80

# Embedding used only for visualization. PAGA/DPT itself is still computed from
# the PCA neighbor graph.
VIS_REDUCTION <- "tsne"
VIS_FALLBACK_REDUCTION <- "umap"
MARKER_TREND_N_BINS <- 20

# Pseudotime windows are retained only as trajectory intervals for summary plots.
WINDOW_OVERLAP <- 0.20
MIN_CELLS_PER_STAGE <- 100

# Cache behavior.
FORCE_RECOMPUTE_PSEUDOTIME <- FALSE
FORCE_RECOMPUTE_WINDOWS <- FALSE

# Python setup. If already configured in RStudio, these can stay unset.
RETICULATE_PYTHON <- Sys.getenv("RETICULATE_PYTHON", unset = NA_character_)
CONDA_ENV <- Sys.getenv("PART3_CONDA_ENV", unset = NA_character_)

# Marker sets for pseudotime orientation/refinement.
MARKER_CONFIG <- list(
  Sertoli = list(
    immature = c("AMH", "JUN", "FOS", "FOSB", "EGR1", "EGR3", "NR4A1", "BTG2", "DUSP1", "RHOB", "GADD45A"),
    mature   = c("CLDN11", "TJP1", "GATA4", "SOX9", "WT1", "HOPX", "INHA", "INHBB", "KITLG", "FSHR")
  ),
  Leydig = list(
    immature = c("TCF21", "PDGFRA", "DLK1", "NR2F2", "LIFR"),
    mature   = c("STAR", "CYP11A1", "CYP17A1", "HSD3B1", "INSL3", "LHCGR", "NR5A1")
  ),
  PMC = list(
    immature = c("PDGFRA", "LUM", "DCN", "COL1A1", "COL3A1"),
    mature   = c("ACTA2", "MYH11", "MYOCD", "CNN1", "TAGLN", "DES")
  ),
  Endothelial = list(
    immature = c("KDR", "SOX17", "ESAM", "EMCN", "PLVAP"),
    mature   = c("PECAM1", "VWF", "CDH5", "ENG", "KLF2", "KLF4")
  ),
  Macrophage = list(
    immature = c("LYZ", "S100A8", "S100A9", "IL1B", "CXCL8"),
    mature   = c("CD68", "CD163", "MRC1", "CSF1R", "MAFB", "SPI1")
  ),
  T_cell = list(
    immature = c("TCF7", "LEF1", "IL7R", "CCR7"),
    mature   = c("CD3D", "CD3E", "TRAC", "TBX21", "GZMK", "GZMB")
  ),
  Spermatogonia = list(
    immature = c("ID4", "FGFR3", "GFRA1", "UTF1", "ZBTB16"),
    mature   = c("MAGEA4", "DAZL", "STRA8", "DMRT1", "KIT")
  ),
  Spermatocyte = list(
    immature = c("STRA8", "DMC1", "SPO11", "SYCP1"),
    mature   = c("SYCP3", "HORMAD1", "MEIOC", "MLH3", "PIWIL1")
  ),
  Spermatid = list(
    immature = c("ACR", "ACRV1", "SPATA18", "TNP1"),
    mature   = c("TNP2", "PRM1", "PRM2", "ODF1", "AKAP4")
  )
)

celltype_order <- CELLTYPES_TO_RUN
FIG4A_MAIN_CELLTYPES <- c("Sertoli", "Leydig", "PMC", "Endothelial", "Spermatogonia")

celltype_palette <- c(
  Sertoli       = "#2C8C8C",
  Leydig        = "#C08A4B",
  PMC           = "#9A9087",
  Endothelial   = "#8CB7C9",
  Macrophage    = "#B0186B",
  T_cell        = "#E6B7C6",
  Spermatogonia = "#BFDDE3",
  Spermatocyte  = "#B78398",
  Spermatid     = "#D95F52"
)

age_palette <- c(
  "2" = "#8CB7C9",
  "5" = "#BFDDE3",
  "8" = "#A7C7A0",
  "11" = "#C7A15A",
  "17" = "#BB6A8A"
)

set.seed(1)

# ------------------------------------------------------------------------------
# 1. Helper functions
# ------------------------------------------------------------------------------

safe_celltype_id <- function(x) tolower(gsub("[^A-Za-z0-9]+", "_", x))

ct_dirs <- function(celltype) {
  id <- safe_celltype_id(celltype)
  list(
    base = file.path(OUTDIR, id),
    tables = file.path(OUTDIR, id, "tables"),
    rds = file.path(OUTDIR, id, "rds")
  )
}

ensure_ct_dirs <- function(celltype) {
  dirs <- ct_dirs(celltype)
  dir.create(dirs$tables, showWarnings = FALSE, recursive = TRUE)
  dir.create(dirs$rds, showWarnings = FALSE, recursive = TRUE)
  dirs
}

get_assay_data_safe <- function(obj, assay = "RNA", slot = "data") {
  if ("layer" %in% names(formals(Seurat::GetAssayData))) {
    tryCatch(
      Seurat::GetAssayData(obj, assay = assay, layer = slot),
      error = function(e) Seurat::GetAssayData(obj, assay = assay, slot = slot)
    )
  } else {
    Seurat::GetAssayData(obj, assay = assay, slot = slot)
  }
}

min_max_scale <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) return(rep(0.5, length(x)))
  (x - rng[1]) / diff(rng)
}

zscore <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

score_markers <- function(obj, genes) {
  expr <- get_assay_data_safe(obj, assay = ASSAY_USE, slot = SLOT_USE)
  genes <- intersect(genes, rownames(expr))
  if (length(genes) == 0) return(rep(NA_real_, ncol(obj)))
  as.numeric(Matrix::colMeans(expr[genes, , drop = FALSE]))
}

ensure_visual_reduction <- function(obj) {
  vis <- tolower(VIS_REDUCTION)
  fallback <- tolower(VIS_FALLBACK_REDUCTION)
  
  if (vis %in% names(obj@reductions)) return(obj)
  
  if (vis == "tsne") {
    obj <- tryCatch(
      RunTSNE(obj, dims = 1:N_PCS, reduction = "pca", reduction.name = "tsne", verbose = FALSE),
      error = function(e) {
        warning("Could not compute t-SNE for visualization: ", e$message)
        obj
      }
    )
    if ("tsne" %in% names(obj@reductions)) return(obj)
  }
  
  if (!(fallback %in% names(obj@reductions)) && fallback == "umap") {
    obj <- tryCatch(
      RunUMAP(obj, dims = 1:N_PCS, reduction = "pca", reduction.name = "umap", verbose = FALSE),
      error = function(e) {
        warning("Could not compute UMAP fallback for visualization: ", e$message)
        obj
      }
    )
  }
  
  obj
}

setup_python <- function() {
  if (!is.na(RETICULATE_PYTHON) && nzchar(RETICULATE_PYTHON)) {
    Sys.setenv(RETICULATE_PYTHON = RETICULATE_PYTHON)
  }
  if (!is.na(CONDA_ENV) && nzchar(CONDA_ENV)) {
    reticulate::use_condaenv(CONDA_ENV, required = FALSE)
  }
  list(sc = reticulate::import("scanpy"))
}

prepare_object <- function(celltype) {
  input_file <- INPUT_MAP[[celltype]]
  if (is.na(input_file) || !file.exists(input_file)) {
    warning("Missing input RDS for ", celltype, ": ", input_file)
    return(NULL)
  }
  
  obj <- readRDS(input_file)
  DefaultAssay(obj) <- ASSAY_USE
  
  missing_cols <- setdiff(c(AGE_COL, DONOR_COL, CLUSTER_COL), colnames(obj@meta.data))
  if (length(missing_cols) > 0) {
    stop(celltype, " is missing metadata columns: ", paste(missing_cols, collapse = ", "))
  }
  
  if (!"pca" %in% names(obj@reductions)) {
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
    obj <- ScaleData(obj, verbose = FALSE)
    obj <- RunPCA(obj, npcs = max(50, N_PCS), verbose = FALSE)
  }
  
  obj <- ensure_visual_reduction(obj)
  obj$part3_cell_type <- celltype
  obj
}

compute_paga_dpt <- function(obj, celltype) {
  py_mod <- setup_python()
  sc <- py_mod$sc
  
  meta <- obj@meta.data
  pca_coords <- Embeddings(obj, "pca")[, 1:min(N_PCS, ncol(Embeddings(obj, "pca"))), drop = FALSE]
  obs <- data.frame(
    cell = rownames(meta),
    age_numeric = meta[[AGE_COL]],
    sample = meta[[DONOR_COL]],
    seurat_cluster = meta[[CLUSTER_COL]],
    row.names = rownames(meta),
    stringsAsFactors = FALSE
  )
  
  adata <- sc$AnnData(X = pca_coords, obs = obs)
  adata$obsm$update(list(X_pca = pca_coords))
  sc$pp$neighbors(adata, n_neighbors = as.integer(N_NEIGHBORS_PAGA), n_pcs = as.integer(ncol(pca_coords)))
  sc$tl$leiden(adata, resolution = RESOLUTION_PAGA)
  sc$tl$paga(adata, groups = "leiden")
  
  py <- reticulate::py
  py$adata <- adata
  reticulate::py_run_string("
import numpy as np
paga_conn_array = adata.uns['paga']['connectivities'].toarray()
leiden_labels = adata.obs['leiden'].values.astype(str)
")
  paga_conn <- py$paga_conn_array
  paga_groups <- py$leiden_labels
  
  root_idx <- select_root_cell(obj, celltype, paga_groups, pca_coords)
  py$root_cell_idx <- root_idx
  reticulate::py_run_string("adata.uns['iroot'] = int(root_cell_idx) - 1")
  
  # Explicit diffmap avoids scanpy's fallback warning before DPT.
  sc$tl$diffmap(adata, n_comps = as.integer(N_DCS))
  sc$tl$dpt(adata, n_dcs = as.integer(N_DCS))
  raw <- as.numeric(reticulate::py_to_r(adata$obs$dpt_pseudotime))
  
  list(
    pseudotime_raw = raw,
    paga_groups = paste0("C", paga_groups),
    paga_connectivities = paga_conn,
    root_cell = colnames(obj)[root_idx]
  )
}

select_root_cell <- function(obj, celltype, paga_groups, pca_coords) {
  meta <- obj@meta.data
  age_num <- suppressWarnings(as.numeric(as.character(meta[[AGE_COL]])))
  min_age <- min(age_num, na.rm = TRUE)
  youngest <- which(age_num == min_age)
  candidates <- youngest
  
  cfg <- MARKER_CONFIG[[celltype]]
  if (!is.null(cfg) && length(cfg$immature) > 0 && length(youngest) >= 10) {
    score <- score_markers(obj, cfg$immature)
    score_y <- score[youngest]
    if (sum(is.finite(score_y)) >= 10) {
      candidates <- youngest[score_y >= quantile(score_y, 0.8, na.rm = TRUE)]
    }
  }
  
  if (length(candidates) == 0) candidates <- youngest
  if (length(candidates) == 1) return(candidates[1])
  
  cluster_distribution <- table(paga_groups[candidates])
  main_cluster <- names(cluster_distribution)[which.max(cluster_distribution)]
  cells_in_cluster <- which(paga_groups == main_cluster)
  use_dims <- 1:min(20, ncol(pca_coords))
  centroid <- colMeans(pca_coords[cells_in_cluster, use_dims, drop = FALSE])
  distances <- apply(
    pca_coords[candidates, use_dims, drop = FALSE],
    1,
    function(x) sqrt(sum((x - centroid)^2))
  )
  candidates[which.min(distances)]
}

age_ordered_pseudotime <- function(raw_pt, age) {
  age_num <- suppressWarnings(as.numeric(as.character(age)))
  ages <- sort(unique(age_num[is.finite(age_num)]))
  out <- rep(NA_real_, length(raw_pt))
  
  for (i in seq_along(ages)) {
    idx <- which(age_num == ages[i] & is.finite(raw_pt))
    if (length(idx) == 0) next
    local <- min_max_scale(raw_pt[idx])
    out[idx] <- (i - 1) / length(ages) + local * (0.8 / length(ages))
  }
  
  min_max_scale(out)
}

marker_refine <- function(obj, celltype, base_pt) {
  cfg <- MARKER_CONFIG[[celltype]]
  if (is.null(cfg)) return(base_pt)
  immature <- score_markers(obj, cfg$immature)
  mature <- score_markers(obj, cfg$mature)
  marker_pt <- min_max_scale(zscore(mature) - zscore(immature))
  min_max_scale(0.7 * base_pt + 0.3 * marker_pt)
}

write_pseudotime_outputs <- function(obj, celltype, paga_res, dirs) {
  meta <- obj@meta.data
  age_num <- suppressWarnings(as.numeric(as.character(meta[[AGE_COL]])))
  pt_raw <- min_max_scale(paga_res$pseudotime_raw)
  pt_age <- age_ordered_pseudotime(pt_raw, age_num)
  pt_final <- marker_refine(obj, celltype, pt_age)
  
  obj$pseudotime_raw <- pt_raw
  obj$pseudotime_age_ordered <- pt_age
  obj$pseudotime_final <- pt_final
  obj$paga_cluster <- paga_res$paga_groups
  
  cfg <- MARKER_CONFIG[[celltype]]
  immature_score <- if (!is.null(cfg)) score_markers(obj, cfg$immature) else NA_real_
  mature_score <- if (!is.null(cfg)) score_markers(obj, cfg$mature) else NA_real_
  
  pt_table <- data.frame(
    cell_type = celltype,
    cell = colnames(obj),
    pseudotime_final = pt_final,
    pseudotime_raw = pt_raw,
    pseudotime_age_ordered = pt_age,
    age_numeric = age_num,
    sample = meta[[DONOR_COL]],
    seurat_cluster = meta[[CLUSTER_COL]],
    paga_cluster = paga_res$paga_groups,
    immature_score = immature_score,
    mature_score = mature_score,
    maturation_score = zscore(mature_score) - zscore(immature_score),
    stringsAsFactors = FALSE
  )
  
  age_summary <- pt_table %>%
    group_by(cell_type, age_numeric) %>%
    summarise(
      n_cells = n(),
      median_pseudotime = median(pseudotime_final, na.rm = TRUE),
      mean_pseudotime = mean(pseudotime_final, na.rm = TRUE),
      q25 = quantile(pseudotime_final, 0.25, na.rm = TRUE),
      q75 = quantile(pseudotime_final, 0.75, na.rm = TRUE),
      median_maturation_score = median(maturation_score, na.rm = TRUE),
      mean_maturation_score = mean(maturation_score, na.rm = TRUE),
      .groups = "drop"
    )
  
  quality <- data.frame(
    cell_type = celltype,
    metric = c("spearman_raw_age", "spearman_age_ordered_age", "spearman_final_age"),
    value = c(
      suppressWarnings(cor(pt_raw, age_num, method = "spearman", use = "complete.obs")),
      suppressWarnings(cor(pt_age, age_num, method = "spearman", use = "complete.obs")),
      suppressWarnings(cor(pt_final, age_num, method = "spearman", use = "complete.obs"))
    )
  )
  
  root_info <- data.frame(
    cell_type = celltype,
    root_cell = paga_res$root_cell,
    stringsAsFactors = FALSE
  )
  
  fwrite(pt_table, file.path(dirs$tables, "pseudotime_values.csv"))
  fwrite(age_summary, file.path(dirs$tables, "pseudotime_age_summary.csv"))
  fwrite(quality, file.path(dirs$tables, "pseudotime_quality_metrics.csv"))
  fwrite(root_info, file.path(dirs$tables, "pseudotime_root_cell.csv"))
  saveRDS(obj, file.path(dirs$rds, "seurat_with_pseudotime_PAGA_DPT.rds"))
  saveRDS(paga_res$paga_connectivities, file.path(dirs$rds, "paga_connectivities.rds"))
  pt_table
}

find_kde_intersection <- function(pt_a, pt_b, n_grid = 4096) {
  pt_a <- pt_a[is.finite(pt_a)]
  pt_b <- pt_b[is.finite(pt_b)]
  if (length(pt_a) < 50 || length(pt_b) < 50) return(NA_real_)
  
  da <- density(pt_a, n = n_grid)
  db <- density(pt_b, n = n_grid)
  x_min <- max(min(da$x), min(db$x))
  x_max <- min(max(da$x), max(db$x))
  grid <- seq(x_min, x_max, length.out = n_grid)
  ya <- approx(da$x, da$y, xout = grid, rule = 2)$y
  yb <- approx(db$x, db$y, xout = grid, rule = 2)$y
  idx <- which(sign((ya - yb)[-1]) * sign((ya - yb)[-length(ya)]) <= 0)
  
  if (length(idx) == 0) return((median(pt_a) + median(pt_b)) / 2)
  cand <- grid[idx]
  cand[which.min(abs(cand - (median(pt_a) + median(pt_b)) / 2))]
}

define_windows <- function(pt_table, dirs) {
  pt <- pt_table$pseudotime_final
  stage <- as.character(pt_table$age_numeric)
  med <- tapply(pt, stage, median, na.rm = TRUE)
  stages <- names(sort(med))
  bounds <- min(pt, na.rm = TRUE)
  
  for (i in seq_len(length(stages) - 1)) {
    pt1 <- pt[stage == stages[i]]
    pt2 <- pt[stage == stages[i + 1]]
    b <- if (sum(is.finite(pt1)) < MIN_CELLS_PER_STAGE || sum(is.finite(pt2)) < MIN_CELLS_PER_STAGE) {
      (median(pt1, na.rm = TRUE) + median(pt2, na.rm = TRUE)) / 2
    } else {
      find_kde_intersection(pt1, pt2)
    }
    if (!is.finite(b)) b <- (median(pt1, na.rm = TRUE) + median(pt2, na.rm = TRUE)) / 2
    bounds <- c(bounds, b)
  }
  
  bounds <- sort(unique(c(bounds, max(pt, na.rm = TRUE))))
  segments <- data.frame(
    seg_id = seq_len(length(bounds) - 1),
    start = bounds[-length(bounds)],
    end = bounds[-1]
  )
  segments$width <- segments$end - segments$start
  half_overlap <- WINDOW_OVERLAP / 2
  
  windows <- segments %>%
    mutate(
      win_id = sprintf("W%02d", seg_id),
      win_start = pmax(min(pt, na.rm = TRUE), start - half_overlap * width),
      win_end = pmin(max(pt, na.rm = TRUE), end + half_overlap * width),
      win_center = (win_start + win_end) / 2
    ) %>%
    dplyr::select(win_id, win_start, win_end, win_center)
  
  fwrite(data.frame(boundary = bounds), file.path(dirs$tables, "pseudotime_kde_boundaries.csv"))
  fwrite(windows, file.path(dirs$tables, "pseudotime_windows.csv"))
  windows
}

summarise_window_cells <- function(pt_table, windows, dirs, celltype) {
  rows <- lapply(seq_len(nrow(windows)), function(i) {
    in_window <- is.finite(pt_table$pseudotime_final) &
      pt_table$pseudotime_final >= windows$win_start[i] &
      pt_table$pseudotime_final <= windows$win_end[i]
    
    data.frame(
      cell_type = celltype,
      win_id = windows$win_id[i],
      win_start = windows$win_start[i],
      win_end = windows$win_end[i],
      win_center = windows$win_center[i],
      n_cells = sum(in_window),
      stringsAsFactors = FALSE
    )
  })
  
  out <- bind_rows(rows)
  fwrite(out, file.path(dirs$tables, "pseudotime_window_cell_counts.csv"))
  out
}

run_one_celltype <- function(celltype) {
  message("\n==============================")
  message("Part 3 multi-celltype pseudotime analysis: ", celltype)
  message("==============================")
  
  dirs <- ensure_ct_dirs(celltype)
  pt_file <- file.path(dirs$tables, "pseudotime_values.csv")
  obj_file <- file.path(dirs$rds, "seurat_with_pseudotime_PAGA_DPT.rds")
  
  if (!FORCE_RECOMPUTE_PSEUDOTIME && file.exists(pt_file) && file.exists(obj_file)) {
    pt_table <- fread(pt_file)
  } else {
    obj <- prepare_object(celltype)
    if (is.null(obj)) return(NULL)
    paga_res <- compute_paga_dpt(obj, celltype)
    pt_table <- write_pseudotime_outputs(obj, celltype, paga_res, dirs)
  }
  
  win_file <- file.path(dirs$tables, "pseudotime_windows.csv")
  if (!FORCE_RECOMPUTE_WINDOWS && file.exists(win_file)) {
    windows <- fread(win_file)
  } else {
    windows <- define_windows(pt_table, dirs)
  }
  
  window_counts <- summarise_window_cells(pt_table, windows, dirs, celltype)
  
  data.frame(
    cell_type = celltype,
    n_cells = nrow(pt_table),
    n_ages = length(unique(pt_table$age_numeric)),
    spearman_final_age = suppressWarnings(cor(
      pt_table$pseudotime_final,
      pt_table$age_numeric,
      method = "spearman",
      use = "complete.obs"
    )),
    n_windows = nrow(windows),
    min_window_cells = min(window_counts$n_cells, na.rm = TRUE),
    median_window_cells = median(window_counts$n_cells, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

collect_table <- function(filename) {
  rows <- lapply(CELLTYPES_TO_RUN, function(ct) {
    f <- file.path(ct_dirs(ct)$tables, filename)
    if (!file.exists(f)) return(NULL)
    fread(f)
  })
  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

marker_config_table <- function() {
  rows <- lapply(names(MARKER_CONFIG), function(ct) {
    cfg <- MARKER_CONFIG[[ct]]
    bind_rows(
      data.frame(cell_type = ct, marker_class = "immature", gene = cfg$immature, stringsAsFactors = FALSE),
      data.frame(cell_type = ct, marker_class = "mature", gene = cfg$mature, stringsAsFactors = FALSE)
    )
  })
  bind_rows(rows)
}

# ------------------------------------------------------------------------------
# 2. Main analysis
# ------------------------------------------------------------------------------

summary_rows <- list()
for (ct in CELLTYPES_TO_RUN) {
  summary_rows[[ct]] <- tryCatch(
    run_one_celltype(ct),
    error = function(e) {
      warning("Part 3 failed for ", ct, ": ", e$message)
      data.frame(
        cell_type = ct,
        n_cells = NA_integer_,
        n_ages = NA_integer_,
        spearman_final_age = NA_real_,
        n_windows = NA_integer_,
        min_window_cells = NA_integer_,
        median_window_cells = NA_real_,
        stringsAsFactors = FALSE
      )
    }
  )
}

summary_all <- bind_rows(summary_rows)
fwrite(summary_all, file.path(TABLE_DIR, "part3_multicelltype_analysis_summary.csv"))
fwrite(marker_config_table(), file.path(TABLE_DIR, "part3_multicelltype_marker_panel.csv"))

aggregate_files <- c(
  "pseudotime_values.csv",
  "pseudotime_age_summary.csv",
  "pseudotime_quality_metrics.csv",
  "pseudotime_windows.csv",
  "pseudotime_window_cell_counts.csv",
  "pseudotime_root_cell.csv"
)

for (f in aggregate_files) {
  agg <- collect_table(f)
  if (!is.null(agg) && nrow(agg) > 0) {
    fwrite(agg, file.path(TABLE_DIR, paste0("all_celltypes_", f)))
  }
}

manifest <- data.frame(
  item = c(
    "Main Figure schematic",
    "Main Figure 4A",
    "Main Figure 4B",
    "Supplementary Figure S4A",
    "Supplementary Figure S4B",
    "Supplementary Figure S4C",
    "Supplementary Figure S4D",
    "Supplementary Table S10",
    "Supplementary Table S11",
    "Supplementary Table S12",
    "Supplementary Table S13",
    "Supplementary Table S14",
    "Supplementary Table S15"
  ),
  source_file = c(
    "manual schematic based on this analysis workflow",
    "all_celltypes_pseudotime_values.csv",
    "part3_multicelltype_analysis_summary.csv / all_celltypes_pseudotime_window_cell_counts.csv",
    "all_celltypes_pseudotime_age_summary.csv",
    "paga_connectivities.rds / seurat_with_pseudotime_PAGA_DPT.rds",
    "all_celltypes_marker_trends.csv",
    "all_celltypes_pseudotime_values.csv",
    "part3_multicelltype_analysis_summary.csv",
    "all_celltypes_pseudotime_age_summary.csv",
    "all_celltypes_pseudotime_values.csv",
    "all_celltypes_marker_trends.csv",
    "part3_multicelltype_marker_panel.csv",
    "part3_multicelltype_paga_network_plot_files.csv"
  ),
  role = c(
    "Explain normal-reference trajectory construction and later disease projection logic",
    "Compare pseudotime organization in the best-supported developmental reference cell types",
    "Assess raw trajectory alignment, final age-calibrated reference, and trajectory interval coverage",
    "Show age-level pseudotime and marker-based maturation changes",
    "Show per-cell-type PAGA trajectory topology and pseudotime embedding",
    "Validate trajectory direction using known immature and mature marker genes",
    "Show the complete all-cell-type pseudotime density panel",
    "Overall trajectory quality and feasibility table",
    "Age-level pseudotime and maturation summary",
    "Cell-level pseudotime and maturation score table",
    "Binned marker expression trends along pseudotime",
    "List marker genes used to orient and validate trajectories",
    "List generated PAGA trajectory network plots"
  )
)
fwrite(manifest, file.path(TABLE_DIR, "part3_multicelltype_figure_table_manifest.csv"))

# ------------------------------------------------------------------------------
# 3. Visualization
# ------------------------------------------------------------------------------

theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 2),
      plot.subtitle = element_text(color = "grey35", size = base_size),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.35),
      strip.background = element_rect(fill = "grey94", color = NA),
      strip.text = element_text(face = "bold", color = "black"),
      legend.title = element_text(face = "bold")
    )
}

theme_cell_style <- function(base_size = 8.5) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 1.8),
      plot.subtitle = element_blank(),
      axis.title = element_text(face = "bold", size = base_size),
      axis.text = element_text(color = "black", size = base_size - 0.5),
      axis.line = element_line(linewidth = 0.35, color = "black"),
      axis.ticks = element_line(linewidth = 0.3, color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", color = "black", size = base_size + 0.5),
      panel.spacing = grid::unit(0.72, "lines"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = base_size - 0.5),
      legend.text = element_text(size = base_size - 0.8),
      plot.margin = margin(5.5, 6, 5.5, 5.5)
    )
}

save_main <- function(name, plot, width, height) {
  ggsave(file.path(MAIN_FIG_DIR, name), plot = plot, width = width, height = height, device = cairo_pdf)
}

save_supp <- function(name, plot, width, height) {
  ggsave(file.path(SUPP_FIG_DIR, name), plot = plot, width = width, height = height, device = cairo_pdf)
}

natural_cluster_order <- function(x) {
  x <- unique(as.character(x))
  numeric_id <- suppressWarnings(as.numeric(sub("^C", "", x)))
  x[order(is.na(numeric_id), numeric_id, x)]
}

build_paga_edges <- function(conn, cluster_levels) {
  conn <- as.matrix(conn)
  cluster_levels <- natural_cluster_order(cluster_levels)
  n_use <- min(nrow(conn), ncol(conn), length(cluster_levels))
  if (!is.finite(n_use) || n_use < 2) {
    return(data.frame(from = character(), to = character(), weight = numeric()))
  }
  
  conn <- conn[seq_len(n_use), seq_len(n_use), drop = FALSE]
  cluster_levels <- cluster_levels[seq_len(n_use)]
  idx <- which(conn > PAGA_EDGE_MIN_WEIGHT & upper.tri(conn), arr.ind = TRUE)
  if (nrow(idx) == 0) {
    return(data.frame(from = character(), to = character(), weight = numeric()))
  }
  
  edges <- data.frame(
    from = cluster_levels[idx[, 1]],
    to = cluster_levels[idx[, 2]],
    weight = conn[idx],
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(weight))
  
  if (!is.null(PAGA_MAX_EDGES) && nrow(edges) > PAGA_MAX_EDGES) {
    edges <- edges[seq_len(PAGA_MAX_EDGES), , drop = FALSE]
  }
  edges
}

select_visual_reduction <- function(obj) {
  preferred <- tolower(VIS_REDUCTION)
  fallback <- tolower(VIS_FALLBACK_REDUCTION)
  reductions <- names(obj@reductions)
  
  if (preferred %in% reductions) return(preferred)
  if (fallback %in% reductions) return(fallback)
  character(0)
}

visual_reduction_label <- function(reduction) {
  dplyr::case_when(
    reduction == "tsne" ~ "t-SNE",
    reduction == "umap" ~ "UMAP",
    TRUE ~ reduction
  )
}

plot_celltype_paga_network <- function(celltype) {
  dirs <- ct_dirs(celltype)
  obj_file <- file.path(dirs$rds, "seurat_with_pseudotime_PAGA_DPT.rds")
  conn_file <- file.path(dirs$rds, "paga_connectivities.rds")
  if (!file.exists(obj_file) || !file.exists(conn_file)) return(NULL)
  
  obj <- readRDS(obj_file)
  if (!(tolower(VIS_REDUCTION) %in% names(obj@reductions))) {
    obj <- ensure_visual_reduction(obj)
    saveRDS(obj, obj_file)
  }
  reduction_use <- select_visual_reduction(obj)
  if (length(reduction_use) == 0) return(NULL)
  reduction_label <- visual_reduction_label(reduction_use)
  
  conn <- readRDS(conn_file)
  emb <- Embeddings(obj, reduction_use)[, 1:2, drop = FALSE]
  meta <- obj@meta.data
  cell_df <- data.frame(
    cell = colnames(obj),
    emb_1 = emb[, 1],
    emb_2 = emb[, 2],
    age_numeric = meta[[AGE_COL]],
    pseudotime_final = meta$pseudotime_final,
    paga_cluster = meta$paga_cluster,
    stringsAsFactors = FALSE
  )
  
  cell_df <- cell_df %>%
    filter(!is.na(paga_cluster), is.finite(emb_1), is.finite(emb_2))
  if (nrow(cell_df) == 0) return(NULL)
  
  nodes <- cell_df %>%
    group_by(paga_cluster) %>%
    summarise(
      x = median(emb_1, na.rm = TRUE),
      y = median(emb_2, na.rm = TRUE),
      n_cells = n(),
      mean_pseudotime = mean(pseudotime_final, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(paga_cluster = as.character(paga_cluster))
  
  cluster_levels <- natural_cluster_order(nodes$paga_cluster)
  nodes$paga_cluster <- factor(nodes$paga_cluster, levels = cluster_levels)
  nodes_join <- nodes %>% mutate(paga_cluster = as.character(paga_cluster))
  
  edges <- build_paga_edges(conn, cluster_levels)
  from_nodes <- nodes_join %>% transmute(from = paga_cluster, x = x, y = y)
  to_nodes <- nodes_join %>% transmute(to = paga_cluster, xend = x, yend = y)
  edge_df <- edges %>%
    left_join(from_nodes, by = "from") %>%
    left_join(to_nodes, by = "to") %>%
    filter(is.finite(x), is.finite(y), is.finite(xend), is.finite(yend))
  
  p_embed <- ggplot(cell_df, aes(x = emb_1, y = emb_2, color = pseudotime_final)) +
    geom_point(size = 0.18, alpha = 0.55) +
    scale_color_gradientn(colors = c("#91B0CB", "white", "#BB6A8A"), name = "Pseudotime") +
    coord_equal() +
    labs(
      x = NULL,
      y = NULL,
      title = paste0(celltype, " cells on ", reduction_label)
    ) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      legend.position = "right"
    )
  
  p_graph <- ggplot()
  if (nrow(edge_df) > 0) {
    p_graph <- p_graph +
      geom_segment(
        data = edge_df,
        aes(x = x, y = y, xend = xend, yend = yend, linewidth = weight),
        color = "grey35",
        alpha = 0.65,
        lineend = "round"
      ) +
      scale_linewidth(range = c(0.25, 2.4), name = "PAGA\nweight")
  }
  
  p_graph <- p_graph +
    geom_point(
      data = nodes,
      aes(x = x, y = y, size = n_cells, fill = mean_pseudotime),
      shape = 21,
      color = "black",
      stroke = 0.25,
      alpha = 0.95
    ) +
    geom_text(data = nodes, aes(x = x, y = y, label = paga_cluster), size = 2.8, vjust = -1.05) +
    scale_size_area(max_size = 8, name = "Cells") +
    scale_fill_gradientn(colors = c("#91B0CB", "white", "#BB6A8A"), name = "Mean\npseudotime") +
    coord_equal() +
    labs(
      x = NULL,
      y = NULL,
      title = paste0(celltype, " PAGA trajectory network")
    ) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      legend.position = "right"
    )
  
  plot_out <- p_embed | p_graph
  out_file <- file.path(PAGA_FIG_DIR, paste0("PAGA_", safe_celltype_id(celltype), "_trajectory_network.pdf"))
  ggsave(out_file, plot = plot_out, width = 10.5, height = 4.8, device = cairo_pdf)
  data.frame(cell_type = celltype, file = out_file, stringsAsFactors = FALSE)
}

summarise_marker_trends <- function(celltype, n_bins = MARKER_TREND_N_BINS) {
  dirs <- ct_dirs(celltype)
  obj_file <- file.path(dirs$rds, "seurat_with_pseudotime_PAGA_DPT.rds")
  if (!file.exists(obj_file)) return(NULL)
  
  cfg <- MARKER_CONFIG[[celltype]]
  if (is.null(cfg)) return(NULL)
  
  marker_df <- bind_rows(
    data.frame(marker_class = "immature", gene = cfg$immature, stringsAsFactors = FALSE),
    data.frame(marker_class = "mature", gene = cfg$mature, stringsAsFactors = FALSE)
  ) %>%
    distinct(marker_class, gene)
  
  obj <- readRDS(obj_file)
  if (!"pseudotime_final" %in% colnames(obj@meta.data)) return(NULL)
  
  expr <- get_assay_data_safe(obj, assay = ASSAY_USE, slot = SLOT_USE)
  marker_df <- marker_df %>% filter(gene %in% rownames(expr))
  if (nrow(marker_df) == 0) return(NULL)
  
  pt <- obj$pseudotime_final
  names(pt) <- colnames(obj)
  ok <- is.finite(pt)
  if (sum(ok) < 50) return(NULL)
  
  pt_ok <- pt[ok]
  breaks <- unique(quantile(pt_ok, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE))
  if (length(breaks) < 3) {
    breaks <- seq(min(pt_ok, na.rm = TRUE), max(pt_ok, na.rm = TRUE), length.out = min(n_bins, length(pt_ok)) + 1)
    breaks <- unique(breaks)
  }
  if (length(breaks) < 3) return(NULL)
  
  bin <- cut(pt_ok, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  cells_ok <- names(pt_ok)
  expr_sub <- as.matrix(expr[marker_df$gene, cells_ok, drop = FALSE])
  
  rows <- lapply(sort(unique(bin[is.finite(bin)])), function(b) {
    cells_b <- cells_ok[bin == b]
    if (length(cells_b) == 0) return(NULL)
    vals <- Matrix::rowMeans(expr_sub[, cells_b, drop = FALSE])
    data.frame(
      cell_type = celltype,
      pseudotime_bin = b,
      bin_center = median(pt[cells_b], na.rm = TRUE),
      n_cells = length(cells_b),
      gene = names(vals),
      mean_expression = as.numeric(vals),
      stringsAsFactors = FALSE
    )
  })
  
  out <- bind_rows(rows) %>%
    left_join(marker_df, by = "gene") %>%
    select(cell_type, marker_class, gene, pseudotime_bin, bin_center, n_cells, mean_expression)
  
  if (nrow(out) > 0) {
    fwrite(out, file.path(dirs$tables, "marker_trends.csv"))
  }
  out
}

plot_marker_trends <- function(marker_trends, celltype) {
  sub <- marker_trends %>% filter(.data$cell_type == celltype)
  if (nrow(sub) == 0) return(NULL)
  
  gene_levels <- sub %>%
    distinct(marker_class, gene) %>%
    arrange(marker_class, gene) %>%
    pull(gene)
  
  sub <- sub %>%
    mutate(
      gene = factor(gene, levels = gene_levels),
      marker_class = factor(marker_class, levels = c("immature", "mature"))
    )
  
  p <- ggplot(sub, aes(x = bin_center, y = mean_expression, color = marker_class, group = gene)) +
    geom_line(linewidth = 0.65, alpha = 0.9) +
    geom_point(size = 0.9, alpha = 0.9) +
    facet_wrap(~gene, scales = "free_y", ncol = 4) +
    scale_color_manual(values = c(immature = "#91B0CB", mature = "#BB6A8A"), name = "Marker class") +
    labs(
      x = "Pseudotime",
      y = "Mean expression",
      title = paste0(celltype, " trajectory marker validation")
    ) +
    theme_pub(9) +
    theme(
      strip.text = element_text(size = 8),
      legend.position = "bottom"
    )
  
  n_genes <- dplyr::n_distinct(sub$gene)
  out_file <- file.path(MARKER_TREND_FIG_DIR, paste0("Marker_trends_", safe_celltype_id(celltype), ".pdf"))
  ggsave(out_file, plot = p, width = 10.5, height = max(4.8, ceiling(n_genes / 4) * 1.8), device = cairo_pdf)
  data.frame(cell_type = celltype, file = out_file, stringsAsFactors = FALSE)
}

summary_all <- fread(file.path(TABLE_DIR, "part3_multicelltype_analysis_summary.csv")) %>%
  mutate(cell_type = factor(cell_type, levels = celltype_order))

raw_rho <- fread(file.path(TABLE_DIR, "all_celltypes_pseudotime_quality_metrics.csv")) %>%
  filter(metric == "spearman_raw_age") %>%
  transmute(
    cell_type = factor(cell_type, levels = celltype_order),
    spearman_raw_age = suppressWarnings(as.numeric(value))
  )

summary_all <- summary_all %>%
  left_join(raw_rho, by = "cell_type") %>%
  mutate(
    trajectory_support = case_when(
      n_ages < 3 ~ "Limited age coverage",
      min_window_cells < 50 ~ "Sparse intervals",
      min_window_cells < 100 ~ "Moderate coverage",
      TRUE ~ "Good coverage"
    ),
    trajectory_support = factor(
      trajectory_support,
      levels = c("Good coverage", "Moderate coverage", "Sparse intervals", "Limited age coverage")
    )
  )

pt_all <- fread(file.path(TABLE_DIR, "all_celltypes_pseudotime_values.csv")) %>%
  mutate(
    cell_type = factor(cell_type, levels = celltype_order),
    age_numeric = as.numeric(age_numeric),
    age_factor = factor(age_numeric, levels = sort(unique(age_numeric)))
  )

window_counts <- fread(file.path(TABLE_DIR, "all_celltypes_pseudotime_window_cell_counts.csv")) %>%
  mutate(cell_type = factor(cell_type, levels = celltype_order))

age_summary <- fread(file.path(TABLE_DIR, "all_celltypes_pseudotime_age_summary.csv")) %>%
  mutate(
    cell_type = factor(cell_type, levels = celltype_order),
    age_numeric = as.numeric(age_numeric)
  )

age_levels_chr <- as.character(sort(unique(pt_all$age_numeric)))
age_colors_use <- age_palette
if (!all(age_levels_chr %in% names(age_colors_use))) {
  age_colors_use <- setNames(colorRampPalette(c("#91B0CB", "#BB6A8A"))(length(age_levels_chr)), age_levels_chr)
}

p_quality_data <- summary_all %>% filter(is.finite(spearman_raw_age))
p_quality_missing <- summary_all %>% filter(!is.finite(spearman_raw_age))
if (nrow(p_quality_missing) > 0) {
  message(
    "Skipping cell types with missing raw DPT-age Spearman rho in Fig4B: ",
    paste(as.character(p_quality_missing$cell_type), collapse = ", ")
  )
}

paga_plot_list <- Filter(Negate(is.null), lapply(CELLTYPES_TO_RUN, plot_celltype_paga_network))
paga_plot_files <- if (length(paga_plot_list) > 0) rbindlist(paga_plot_list, use.names = TRUE, fill = TRUE) else NULL
if (!is.null(paga_plot_files) && nrow(paga_plot_files) > 0) {
  fwrite(paga_plot_files, file.path(TABLE_DIR, "part3_multicelltype_paga_network_plot_files.csv"))
}

marker_trend_list <- lapply(CELLTYPES_TO_RUN, summarise_marker_trends)
marker_trend_list <- Filter(Negate(is.null), marker_trend_list)
marker_trends <- if (length(marker_trend_list) > 0) rbindlist(marker_trend_list, use.names = TRUE, fill = TRUE) else NULL
if (!is.null(marker_trends) && nrow(marker_trends) > 0) {
  fwrite(marker_trends, file.path(TABLE_DIR, "all_celltypes_marker_trends.csv"))
  marker_plot_list <- Filter(Negate(is.null), lapply(CELLTYPES_TO_RUN, function(ct) plot_marker_trends(marker_trends, ct)))
  marker_plot_files <- if (length(marker_plot_list) > 0) rbindlist(marker_plot_list, use.names = TRUE, fill = TRUE) else NULL
  if (!is.null(marker_plot_files) && nrow(marker_plot_files) > 0) {
    fwrite(marker_plot_files, file.path(TABLE_DIR, "part3_multicelltype_marker_trend_plot_files.csv"))
  }
}

pt_medians <- pt_all %>%
  group_by(cell_type, age_factor) %>%
  summarise(
    median_pseudotime = median(pseudotime_final, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  )

make_fig4a_density <- function(plot_data, median_data, ncol, title = NULL) {
  ggplot(plot_data, aes(x = pseudotime_final, y = age_factor, fill = age_factor)) +
    ggridges::geom_density_ridges(
      scale = 0.92,
      rel_min_height = 0.015,
      alpha = 0.88,
      color = "white",
      linewidth = 0.18
    ) +
    geom_point(
      data = median_data,
      aes(x = median_pseudotime, y = age_factor),
      inherit.aes = FALSE,
      shape = 21,
      size = 1.35,
      stroke = 0.28,
      color = "black",
      fill = "white"
    ) +
    facet_wrap(~cell_type, ncol = ncol, scales = "free_y") +
    scale_fill_manual(values = age_colors_use, name = "Age") +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.5, 1),
      labels = c("0", "0.5", "1"),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    coord_cartesian(clip = "on") +
    labs(
      x = "Developmental pseudotime",
      y = "Age",
      title = title
    ) +
    theme_cell_style(8.4) +
    theme(
      plot.title = element_text(size = 10.5, face = "bold"),
      legend.position = "bottom",
      legend.key.width = grid::unit(0.55, "cm"),
      legend.margin = margin(t = -2, b = -3),
      panel.spacing.x = grid::unit(0.58, "lines"),
      panel.spacing.y = grid::unit(0.7, "lines"),
      strip.text = element_text(size = 8.8, face = "bold"),
      axis.title.x = element_text(margin = margin(t = 4)),
      axis.title.y = element_text(margin = margin(r = 4))
    )
}

fig4a_main_data <- pt_all %>%
  filter(as.character(cell_type) %in% FIG4A_MAIN_CELLTYPES) %>%
  mutate(cell_type = factor(as.character(cell_type), levels = FIG4A_MAIN_CELLTYPES))

fig4a_main_medians <- pt_medians %>%
  filter(as.character(cell_type) %in% FIG4A_MAIN_CELLTYPES) %>%
  mutate(cell_type = factor(as.character(cell_type), levels = FIG4A_MAIN_CELLTYPES))

# Figure 4A: Cell-style main panel with the best-supported developmental reference cell types.
p_density <- make_fig4a_density(
  fig4a_main_data,
  fig4a_main_medians,
  ncol = length(FIG4A_MAIN_CELLTYPES),
  title = "Developmental pseudotime reference in normal testicular cells"
)

save_main("Fig4A_multicelltype_pseudotime_density.pdf", p_density, width = 9.4, height = 2.95)

# Full all-cell-type version for supplementary review.
p_density_all <- make_fig4a_density(
  pt_all,
  pt_medians,
  ncol = 3,
  title = "Age-resolved developmental pseudotime across all analyzed cell types"
)

save_supp("FigS4D_all_celltypes_pseudotime_density.pdf", p_density_all, width = 10.8, height = 8.3)

# Figure 4B: raw trajectory alignment and interval coverage.
support_palette <- c(
  "Good coverage" = "#2C8C8C",
  "Moderate coverage" = "#C7A15A",
  "Sparse intervals" = "#BB6A8A",
  "Limited age coverage" = "grey70"
)

p_quality <- ggplot(
  p_quality_data,
  aes(x = spearman_raw_age, y = reorder(cell_type, spearman_raw_age), fill = trajectory_support)
) +
  geom_vline(xintercept = 0, color = "grey75", linewidth = 0.35) +
  geom_col(width = 0.68, color = "white", linewidth = 0.25) +
  geom_point(aes(x = spearman_final_age), shape = 21, fill = "white", color = "black", size = 2.2, stroke = 0.35) +
  scale_fill_manual(values = support_palette, name = "Coverage") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(
    x = "Spearman rho",
    y = NULL,
    title = "Trajectory alignment with developmental age",
    subtitle = "Bars: raw PAGA/DPT pseudotime; open circles: final age-calibrated reference"
  ) +
  theme_pub(9.5) +
  theme(legend.position = "bottom")

coverage_data <- summary_all %>%
  filter(is.finite(min_window_cells), is.finite(median_window_cells), min_window_cells > 0, median_window_cells > 0)

p_window_cells <- ggplot(
  coverage_data,
  aes(y = reorder(cell_type, median_window_cells), color = trajectory_support)
) +
  geom_segment(aes(x = min_window_cells, xend = median_window_cells, yend = reorder(cell_type, median_window_cells)),
               linewidth = 1.05, alpha = 0.72, lineend = "round") +
  geom_point(aes(x = min_window_cells), shape = 21, fill = "white", size = 2.1, stroke = 0.35) +
  geom_point(aes(x = median_window_cells), size = 2.3) +
  scale_color_manual(values = support_palette, name = "Coverage") +
  scale_x_log10(
    breaks = c(5, 10, 50, 100, 500, 1000, 5000),
    labels = scales::comma
  ) +
  labs(
    x = "Cells per pseudotime interval (log10)",
    y = NULL,
    title = "Trajectory interval coverage",
    subtitle = "Open circle: minimum; filled circle: median"
  ) +
  theme_pub(9.5) +
  theme(legend.position = "none")

fig4b <- p_quality | p_window_cells
save_main("Fig4B_multicelltype_pseudotime_quality_and_coverage.pdf", fig4b, width = 12.4, height = 5.2)

# Supplementary: median pseudotime and maturation score by age.
p_age_pt <- ggplot(age_summary, aes(x = age_numeric, y = median_pseudotime, color = cell_type, group = cell_type)) +
  geom_line(linewidth = 0.75, alpha = 0.85) +
  geom_point(size = 1.8, alpha = 0.9) +
  scale_color_manual(values = celltype_palette, name = "Cell type") +
  labs(
    x = "Age",
    y = "Median pseudotime",
    title = "Median pseudotime by age"
  ) +
  theme_pub(10)

p_age_maturation <- ggplot(age_summary, aes(x = age_numeric, y = median_maturation_score, color = cell_type, group = cell_type)) +
  geom_hline(yintercept = 0, color = "grey82", linewidth = 0.35) +
  geom_line(linewidth = 0.75, alpha = 0.85) +
  geom_point(size = 1.8, alpha = 0.9) +
  scale_color_manual(values = celltype_palette, name = "Cell type") +
  labs(
    x = "Age",
    y = "Median maturation score",
    title = "Marker-based maturation score by age"
  ) +
  theme_pub(10)

save_supp("FigS4A_multicelltype_age_pseudotime_and_maturation.pdf", p_age_pt | p_age_maturation, width = 12, height = 4.8)

placement <- data.frame(
  placement = c(
    "Main Figure schematic",
    "Main Figure 4A",
    "Main Figure 4B",
    "Supplementary Figure S4A",
    "Supplementary Figure S4B",
    "Supplementary Figure S4C",
    "Supplementary Figure S4D",
    "Supplementary Table S10",
    "Supplementary Table S11",
    "Supplementary Table S12",
    "Supplementary Table S13",
    "Supplementary Table S14",
    "Supplementary Table S15"
  ),
  file = c(
    "manual schematic: normal reference trajectory construction and disease projection concept",
    "figures/main/Fig4A_multicelltype_pseudotime_density.pdf",
    "figures/main/Fig4B_multicelltype_pseudotime_quality_and_coverage.pdf",
    "figures/supplementary/FigS4A_multicelltype_age_pseudotime_and_maturation.pdf",
    "figures/supplementary/paga_networks/PAGA_*_trajectory_network.pdf",
    "figures/supplementary/marker_trends/Marker_trends_*.pdf",
    "figures/supplementary/FigS4D_all_celltypes_pseudotime_density.pdf",
    "tables/part3_multicelltype_analysis_summary.csv",
    "tables/all_celltypes_pseudotime_age_summary.csv",
    "tables/all_celltypes_pseudotime_values.csv",
    "tables/all_celltypes_marker_trends.csv",
    "tables/part3_multicelltype_marker_panel.csv",
    "tables/part3_multicelltype_paga_network_plot_files.csv"
  ),
  rationale = c(
    "Introduce the Part 3 logic: normal trajectories define the reference used for later disease projection.",
    "Compare developmental pseudotime organization in the best-supported reference cell types.",
    "Identify which cell types have raw trajectory age alignment and sufficient interval coverage.",
    "Show age-level changes in pseudotime and marker-based maturation score.",
    "Show per-cell-type PAGA trajectory topology and pseudotime embedding.",
    "Validate trajectory direction with known immature and mature marker genes.",
    "Show the complete all-cell-type pseudotime density panel.",
    "Overall trajectory quality and coverage summary.",
    "Age-level pseudotime and maturation summary.",
    "Cell-level pseudotime and maturation score table.",
    "Binned marker expression trends along pseudotime.",
    "Marker panel used for pseudotime orientation and validation.",
    "List generated PAGA trajectory network plots."
  )
)

fwrite(placement, file.path(TABLE_DIR, "part3_multicelltype_main_supplementary_figure_table_plan.csv"))

message("Part 3 multi-celltype pseudotime analysis finished: ", normalizePath(OUTDIR, winslash = "/"))

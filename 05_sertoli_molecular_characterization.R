#!/usr/bin/env Rscript

# ==============================================================================
# 05_sertoli_molecular_characterization.R
#
# Purpose:
#   Use the Part 4 disease projection results to focus on Sertoli cells and test
#   whether iNOA-associated developmental immaturity is accompanied by molecular
#   changes in stress-response, senescence-associated, metabolic and Sertoli
#   support-related programs. The useful non-redundant parts of the previous
#   Part 7 script are folded in here as an extension: multicelltype
#   prioritization and expanded molecular-program summaries.
#
# Main interpretation:
#   Part 4 indicates broad disease-associated developmental shifts across cell
#   types, but Sertoli cells are prioritized here because their early-state
#   projection is consistent across iNOA samples and directly relates to the
#   somatic support niche required for spermatogenesis.
#
# This script intentionally starts from Part 4 outputs. It does not rerun
# pseudotime projection and does not run inferCSN. The analysis is framed as
# molecular characterization/supporting evidence rather than an independent
# validation experiment because Part 4 and Part 5 use related expression data.
# The folded Part 7 extension should be interpreted as molecular-program
# exploration, not as a separate causal or independent-validation layer.
#
# Main outputs:
#   results/05_sertoli_molecular_characterization/tables/
#   results/05_sertoli_molecular_characterization/figures/
#
# Default paths can be overridden with environment variables:
#   PROJECT_DIR, PART5_BASE_DIR, PART5_PART3_DIR, PART5_PART4_DIR,
#   PART5_OUTPUT_DIR, PART5_iNOA_RDS, PART5_AZFA_RDS
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

HAS_GGREPEL <- requireNamespace("ggrepel", quietly = TRUE)
HAS_CLUSTERPROFILER <- requireNamespace("clusterProfiler", quietly = TRUE)
HAS_ORGDB <- requireNamespace("org.Hs.eg.db", quietly = TRUE)
HAS_ENRICHPLOT <- requireNamespace("enrichplot", quietly = TRUE)

# ------------------------------------------------------------------------------
# 0. User configuration
# ------------------------------------------------------------------------------

is_absolute_path <- function(path) {
  grepl("^[A-Za-z]:[/\\\\]", path) || startsWith(path, "/") || startsWith(path, "\\\\")
}

resolve_dir <- function(path, base_dir, legacy_path = NULL) {
  out <- if (is_absolute_path(path)) path else file.path(base_dir, path)
  if (!dir.exists(out) && !is.null(legacy_path)) {
    legacy <- if (is_absolute_path(legacy_path)) legacy_path else file.path(base_dir, legacy_path)
    if (dir.exists(legacy)) out <- legacy
  }
  out
}

BASE_DIR <- Sys.getenv("PART5_BASE_DIR", unset = Sys.getenv("PROJECT_DIR", unset = getwd()))
BASE_DIR <- normalizePath(BASE_DIR, winslash = "/", mustWork = FALSE)

PART3_DIR <- resolve_dir(
  Sys.getenv("PART5_PART3_DIR", unset = file.path("results", "03_developmental_reference_construction")),
  BASE_DIR,
  legacy_path = "part3_multicelltype_pseudotime_results"
)
PART4_DIR <- resolve_dir(
  Sys.getenv("PART5_PART4_DIR", unset = file.path("results", "04_disease_cell_projection")),
  BASE_DIR,
  legacy_path = "part4_disease_projection_results"
)

OUTDIR_SETTING <- Sys.getenv("PART5_OUTPUT_DIR", unset = file.path("results", "05_sertoli_molecular_characterization"))
OUTDIR <- if (is_absolute_path(OUTDIR_SETTING)) OUTDIR_SETTING else file.path(BASE_DIR, OUTDIR_SETTING)
TABLE_DIR <- file.path(OUTDIR, "tables")
FIG_DIR <- file.path(OUTDIR, "figures")
MAIN_FIG_DIR <- file.path(FIG_DIR, "main")
SUPP_FIG_DIR <- file.path(FIG_DIR, "supplementary")

dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(MAIN_FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(SUPP_FIG_DIR, showWarnings = FALSE, recursive = TRUE)

ASSAY_USE <- "RNA"
SLOT_USE <- "data"
AGE_COL <- "age_numeric"
SERTOLI_CELLTYPE <- "Sertoli"

PART4_PROJECTED_VALUES <- file.path(
  PART4_DIR,
  "tables",
  "disease_projected_pseudotime_values.csv"
)

PART4_QUALITY_TABLE <- file.path(
  PART4_DIR,
  "tables",
  "disease_projection_quality_by_celltype.csv"
)

PART4_SAMPLE_SUMMARY <- file.path(
  PART4_DIR,
  "tables",
  "disease_projection_summary_by_sample.csv"
)

DISEASE_OBJECT_CANDIDATES <- list(
  iNOA = c(
    Sys.getenv("PART5_iNOA_RDS", unset = ""),
    file.path(BASE_DIR, "NOA_results", "iNOA_merged_final.rds"),
    file.path(BASE_DIR, "iNOA_merged_final.rds")
  ),
  AZFa_Del = c(
    Sys.getenv("PART5_AZFA_RDS", unset = ""),
    file.path(BASE_DIR, "NOA_results", "AZFa_Del_1_processed.rds"),
    file.path(BASE_DIR, "AZFa_Del_1_processed.rds")
  )
)

LOW_PSEUDOTIME_QUANTILE <- 0.33
HIGH_PSEUDOTIME_QUANTILE <- 0.67
MIN_CELLS_PER_GROUP <- 20
MIN_PCT_DEG <- 0.05
LOG2FC_CUTOFF <- 0.25
FDR_CUTOFF <- 0.05
MAX_CELLS_PER_GROUP_FOR_DEG <- 1200
USE_PART4_QUALITY_FILTER <- TRUE
PART4_QUALITY_KEEP <- c("high", "moderate")
GROUP_INOA_WITHIN_SAMPLE <- TRUE
MIN_CELLS_PER_SAMPLE_SUBGROUP <- 30
SAMPLE_CONSISTENCY_MIN_FRACTION <- 2 / 3
VOLCANO_Y_CAP <- 160
VOLCANO_IMMATURE_LABEL_GENES <- c(
  "DDIT3", "HSPA5", "HSPA6", "CXCL2", "CTGF",
  "EGR1", "EGR2", "EGR3", "DUSP2", "DUSP5",
  "IER3", "THBS1"
)
VOLCANO_MATURE_LABEL_GENES <- c(
  "HOPX", "CALM1", "CALM2", "PEG10", "RAB31", "GSTM3", "FTH1"
)
set.seed(20260429)

analysis_group_order <- c(
  "Normal immature reference",
  "Normal mature reference",
  "AZFa_Del Sertoli",
  "iNOA immature-like Sertoli",
  "iNOA intermediate Sertoli",
  "iNOA mature-like Sertoli"
)

analysis_group_palette <- c(
  "Normal immature reference" = "#8CB7C9",
  "Normal mature reference" = "#2C8C8C",
  "AZFa_Del Sertoli" = "#D95F52",
  "iNOA immature-like Sertoli" = "#B0186B",
  "iNOA intermediate Sertoli" = "#C08A4B",
  "iNOA mature-like Sertoli" = "#6E6BB8"
)

SENESCENCE_MARKERS <- list(
  cell_cycle_arrest = c("CDKN1A", "CDKN2A", "CDKN1B", "CDKN2B", "TP53", "RB1"),
  sasp_factors = c("IL6", "CXCL8", "IL8", "CXCL1", "CXCL2", "CCL2", "MMP1", "MMP3", "ICAM1"),
  dna_damage = c("ATM", "ATR", "CHEK1", "CHEK2", "H2AFX", "MDC1"),
  autophagy_senescence = c("BECN1", "ATG5", "ATG7", "MAP1LC3B", "SQSTM1")
)

STRESS_MARKERS <- list(
  oxidative_stress = c("SOD1", "SOD2", "CAT", "GPX1", "GPX4", "PRDX1", "NQO1", "HMOX1"),
  er_stress = c("HSPA5", "DDIT3", "ATF4", "ATF6", "XBP1", "ERN1"),
  heat_shock = c("HSP90AA1", "HSPA1A", "HSPA8", "HSPB1", "DNAJB1"),
  hypoxia = c("HIF1A", "EPAS1", "VEGFA", "SLC2A1", "LDHA")
)

FUNCTIONAL_GENES <- list(
  spermatogenesis_support = c("GDNF", "SCF", "KITLG", "BMP4", "FGF2", "NGF"),
  blood_testis_barrier = c("OCLN", "CLDN11", "TJP1", "CDH2", "CTNNB1"),
  hormone_signaling = c("AR", "FSHR", "ESR1", "NR5A1", "CYP19A1"),
  structural_identity = c("SOX9", "WT1", "AMH", "GATA4", "DHH", "DMRT1"),
  metabolic_support = c("LDHA", "LDHB", "G6PD", "PFKP", "HK2")
)

MAIN_GENE_SETS <- list(
  Senescence = unique(unlist(SENESCENCE_MARKERS)),
  `Stress response` = unique(unlist(STRESS_MARKERS)),
  `Sertoli function` = unique(unlist(FUNCTIONAL_GENES))
)

EXTENDED_GENE_SETS <- list(
  `Immature developmental state` = c("AMH", "SOX9", "WT1", "GATA4", "DHH", "DMRT1"),
  `Sertoli support function` = c("GDNF", "KITLG", "SCF", "BMP4", "FGF2", "CLDN11", "OCLN", "TJP1", "CDH2", "AR", "FSHR", "INHBB"),
  `Stress response` = unique(unlist(STRESS_MARKERS)),
  `Senescence SASP` = c("CDKN1A", "CDKN2A", "CDKN1B", "TP53", "RB1", "IL6", "IL8", "CXCL1", "CXCL2", "CCL2", "MMP1", "MMP3", "ICAM1"),
  `Metabolic support` = c("LDHA", "LDHB", "G6PD", "PFKP", "HK2", "SLC2A1", "SLC16A1", "SLC16A3")
)

SUBCATEGORY_GENE_SETS <- c(
  setNames(SENESCENCE_MARKERS, paste0("Senescence: ", names(SENESCENCE_MARKERS))),
  setNames(STRESS_MARKERS, paste0("Stress: ", names(STRESS_MARKERS))),
  setNames(FUNCTIONAL_GENES, paste0("Function: ", names(FUNCTIONAL_GENES)))
)

REPRESENTATIVE_GENES <- unique(c(
  "CDKN1A", "CDKN2A", "IL6", "CXCL8", "HMOX1", "HSPA5",
  "KITLG", "CLDN11", "SOX9", "AMH", "FSHR", "GDNF"
))

# ------------------------------------------------------------------------------
# 1. Helper functions
# ------------------------------------------------------------------------------

safe_celltype_id <- function(x) tolower(gsub("[^A-Za-z0-9]+", "_", x))

ct_dirs_part3 <- function(celltype) {
  id <- safe_celltype_id(celltype)
  list(
    tables = file.path(PART3_DIR, id, "tables"),
    rds = file.path(PART3_DIR, id, "rds")
  )
}

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

theme_cell <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(color = "black"),
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = element_text(size = base_size, color = "grey35", hjust = 0),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.35),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", color = "black"),
      legend.title = element_text(face = "bold"),
      legend.key.height = grid::unit(0.35, "cm"),
      legend.key.width = grid::unit(0.35, "cm")
    )
}

short_group_labels <- c(
  "Normal immature reference" = "Normal\nimmature",
  "Normal mature reference" = "Normal\nmature",
  "AZFa_Del Sertoli" = "AZFa_Del\nSertoli",
  "iNOA immature-like Sertoli" = "iNOA\nimmature-like",
  "iNOA intermediate Sertoli" = "iNOA\nintermediate",
  "iNOA mature-like Sertoli" = "iNOA\nmature-like"
)

short_program_labels <- c(
  Senescence = "Senescence",
  `Stress response` = "Stress\nresponse",
  `Sertoli function` = "Sertoli\nfunction",
  `Immature developmental state` = "Immature\ndevelopment",
  `Sertoli support function` = "Sertoli\nsupport",
  `Senescence SASP` = "Senescence\nSASP",
  `Metabolic support` = "Metabolic\nsupport"
)

wrap_text <- function(x, width = 34) {
  vapply(x, function(xx) paste(strwrap(xx, width = width), collapse = "\n"), character(1))
}

parse_gene_ratio <- function(x) {
  vapply(strsplit(as.character(x), "/", fixed = TRUE), function(z) {
    if (length(z) != 2) return(NA_real_)
    as.numeric(z[1]) / as.numeric(z[2])
  }, numeric(1))
}

to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

first_existing_col <- function(candidates, df) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) NA_character_ else hit[[1]]
}

save_main <- function(name, plot, width, height) {
  ggsave(file.path(MAIN_FIG_DIR, name), plot = plot, width = width, height = height, device = cairo_pdf)
}

save_supp <- function(name, plot, width, height) {
  ggsave(file.path(SUPP_FIG_DIR, name), plot = plot, width = width, height = height, device = cairo_pdf)
}

get_assay_data_safe <- function(obj, assay = "RNA", slot = "data") {
  assay_obj <- obj[[assay]]
  layer_names <- tryCatch(SeuratObject::Layers(assay_obj), error = function(e) character(0))
  if (length(layer_names) > 0) {
    layer_use <- layer_names[layer_names == slot | startsWith(layer_names, paste0(slot, "."))]
    if (length(layer_use) == 1) {
      return(SeuratObject::LayerData(assay_obj, layer = layer_use))
    }
    if (length(layer_use) > 1) {
      mats <- lapply(layer_use, function(ly) SeuratObject::LayerData(assay_obj, layer = ly))
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

load_sertoli_reference <- function() {
  f <- file.path(ct_dirs_part3(SERTOLI_CELLTYPE)$rds, "seurat_with_pseudotime_PAGA_DPT.rds")
  if (!file.exists(f)) stop("Missing Part 3 Sertoli reference: ", f)
  obj <- readRDS(f)
  DefaultAssay(obj) <- ASSAY_USE
  obj <- join_layers_safe(obj, ASSAY_USE)
  if (!"pseudotime_final" %in% colnames(obj@meta.data)) {
    stop("Part 3 Sertoli reference lacks pseudotime_final.")
  }
  
  pt_file <- file.path(ct_dirs_part3(SERTOLI_CELLTYPE)$tables, "pseudotime_values.csv")
  if (file.exists(pt_file)) {
    pt_table <- fread(pt_file)
    add_cols <- intersect(
      c("immature_score", "mature_score", "maturation_score"),
      colnames(pt_table)
    )
    if (length(add_cols) > 0 && "cell" %in% colnames(pt_table)) {
      idx <- match(colnames(obj), pt_table$cell)
      for (cc in add_cols) obj[[cc]] <- pt_table[[cc]][idx]
    }
  }
  
  obj
}

load_part4_sertoli_projection <- function() {
  if (!file.exists(PART4_PROJECTED_VALUES)) {
    stop("Missing Part 4 projected values table: ", PART4_PROJECTED_VALUES)
  }
  projected <- fread(PART4_PROJECTED_VALUES)
  required <- c("cell", "cell_type", "disease_group", "projected_pseudotime")
  missing <- setdiff(required, colnames(projected))
  if (length(missing) > 0) {
    stop("Part 4 projected values table lacks required columns: ", paste(missing, collapse = ", "))
  }
  if (!"sample_id" %in% colnames(projected)) projected$sample_id <- projected$disease_group
  if (!"projection_quality" %in% colnames(projected)) {
    projected$projection_quality <- "not_recorded"
    warning("Part 4 projected values table lacks projection_quality; Part 5 quality filtering is disabled.")
  }
  if (!"mean_knn_distance" %in% colnames(projected)) projected$mean_knn_distance <- NA_real_
  if (!"nearest_normal_age" %in% colnames(projected)) projected$nearest_normal_age <- NA_real_
  if (!"inferred_normal_equiv_age" %in% colnames(projected)) projected$inferred_normal_equiv_age <- NA_real_
  
  projected <- projected %>%
    filter(cell_type == SERTOLI_CELLTYPE, is.finite(projected_pseudotime)) %>%
    mutate(
      projection_quality = tolower(as.character(projection_quality)),
      part5_quality_keep = !USE_PART4_QUALITY_FILTER |
        projection_quality %in% tolower(PART4_QUALITY_KEEP) |
        projection_quality == "not_recorded"
    ) %>%
    as.data.frame()
  
  quality_summary <- projected %>%
    group_by(disease_group, sample_id, projection_quality) %>%
    summarise(n_cells = n(), .groups = "drop") %>%
    group_by(disease_group, sample_id) %>%
    mutate(fraction = n_cells / sum(n_cells)) %>%
    ungroup()
  fwrite(quality_summary, file.path(TABLE_DIR, "part5_sertoli_projection_quality_by_sample.csv"))
  
  retained <- projected %>% filter(part5_quality_keep)
  filter_summary <- data.frame(
    cell_type = SERTOLI_CELLTYPE,
    use_part4_quality_filter = USE_PART4_QUALITY_FILTER,
    retained_quality_classes = paste(PART4_QUALITY_KEEP, collapse = ";"),
    n_cells_before_filter = nrow(projected),
    n_cells_after_filter = nrow(retained),
    retained_fraction = ifelse(nrow(projected) == 0, NA_real_, nrow(retained) / nrow(projected)),
    stringsAsFactors = FALSE
  )
  fwrite(filter_summary, file.path(TABLE_DIR, "part5_sertoli_projection_quality_filter_summary.csv"))
  
  if (nrow(retained) < MIN_CELLS_PER_GROUP) {
    warning("Too few high/moderate-quality projected Sertoli cells after filtering; falling back to all projected Sertoli cells.")
    retained <- projected
  }
  if (nrow(retained) < MIN_CELLS_PER_GROUP) {
    stop("Too few projected disease Sertoli cells in Part 4 table.")
  }
  retained
}

load_disease_objects <- function() {
  objs <- list()
  for (nm in names(DISEASE_OBJECT_CANDIDATES)) {
    candidates <- DISEASE_OBJECT_CANDIDATES[[nm]]
    f <- candidates[file.exists(candidates)][1]
    if (is.na(f) || length(f) == 0) {
      warning("Missing disease object for ", nm, ". Tried: ", paste(candidates, collapse = " ; "))
      next
    }
    message("Loading disease object for ", nm, ": ", normalizePath(f, winslash = "/"))
    obj <- readRDS(f)
    DefaultAssay(obj) <- ASSAY_USE
    obj <- join_layers_safe(obj, ASSAY_USE)
    obj$disease_group <- nm
    if (!"sample_id" %in% colnames(obj@meta.data)) obj$sample_id <- nm
    objs[[nm]] <- obj
  }
  if (length(objs) == 0) stop("No disease objects were loaded.")
  objs
}

combine_disease_objects <- function(objs) {
  if (length(objs) == 1) return(join_layers_safe(objs[[1]], ASSAY_USE))
  combined <- merge(objs[[1]], y = objs[-1], add.cell.ids = names(objs), project = "Disease_combined")
  join_layers_safe(combined, ASSAY_USE)
}

strip_known_disease_prefix <- function(x) {
  sub("^(iNOA|AZFa_Del)_", "", x)
}

match_projection_cells <- function(projected, object_cells) {
  exact <- match(projected$cell, object_cells)
  matched <- object_cells[exact]
  
  missing <- is.na(matched)
  if (any(missing)) {
    projected_stripped <- strip_known_disease_prefix(projected$cell[missing])
    object_stripped <- strip_known_disease_prefix(object_cells)
    fallback <- match(projected_stripped, object_stripped)
    matched[missing] <- object_cells[fallback]
  }
  
  projected$matched_cell <- matched
  projected <- projected[!is.na(projected$matched_cell), , drop = FALSE]
  projected <- projected[!duplicated(projected$matched_cell), , drop = FALSE]
  projected
}

attach_projection_metadata <- function(disease_obj, projected) {
  projected <- match_projection_cells(projected, colnames(disease_obj))
  if (nrow(projected) < MIN_CELLS_PER_GROUP) {
    stop("Could not match enough Part 4 Sertoli projected cells to loaded disease objects.")
  }
  
  meta_cols <- setdiff(colnames(projected), "matched_cell")
  for (cc in meta_cols) {
    new_col <- paste0("part4_", cc)
    disease_obj@meta.data[[new_col]] <- NA
    disease_obj@meta.data[projected$matched_cell, new_col] <- projected[[cc]]
  }
  
  subset(disease_obj, cells = projected$matched_cell)
}

assign_normal_reference_groups <- function(ref) {
  ref_age <- suppressWarnings(as.numeric(as.character(ref@meta.data[[AGE_COL]])))
  ref_pt <- ref$pseudotime_final
  
  if (any(is.finite(ref_age))) {
    mature_age <- max(ref_age[is.finite(ref_age)], na.rm = TRUE)
    immature_age <- min(ref_age[is.finite(ref_age)], na.rm = TRUE)
    group <- ifelse(
      ref_age == mature_age,
      "Normal mature reference",
      ifelse(ref_age == immature_age, "Normal immature reference", "Normal other reference")
    )
  } else {
    q_low <- quantile(ref_pt, LOW_PSEUDOTIME_QUANTILE, na.rm = TRUE)
    q_high <- quantile(ref_pt, HIGH_PSEUDOTIME_QUANTILE, na.rm = TRUE)
    group <- ifelse(
      ref_pt >= q_high,
      "Normal mature reference",
      ifelse(ref_pt <= q_low, "Normal immature reference", "Normal other reference")
    )
  }
  
  ref$analysis_group <- group
  ref
}

assign_disease_sertoli_groups <- function(disease_sertoli) {
  meta <- disease_sertoli@meta.data
  pt <- suppressWarnings(as.numeric(meta$part4_projected_pseudotime))
  disease_group <- as.character(meta$disease_group)
  sample_id <- as.character(meta$sample_id)
  if (all(is.na(sample_id)) || all(!nzchar(sample_id))) sample_id <- disease_group
  inoa_pt <- pt[disease_group == "iNOA" & is.finite(pt)]
  
  if (length(inoa_pt) < MIN_CELLS_PER_GROUP) {
    stop("Too few iNOA Sertoli cells for pseudotime subgrouping.")
  }
  
  group <- rep(NA_character_, nrow(meta))
  group[disease_group == "AZFa_Del"] <- "AZFa_Del Sertoli"
  threshold_rows <- list()
  
  if (GROUP_INOA_WITHIN_SAMPLE) {
    for (sid in unique(sample_id[disease_group == "iNOA"])) {
      idx <- which(disease_group == "iNOA" & sample_id == sid & is.finite(pt))
      if (length(idx) >= MIN_CELLS_PER_SAMPLE_SUBGROUP * 3) {
        q_low <- unname(quantile(pt[idx], LOW_PSEUDOTIME_QUANTILE, na.rm = TRUE))
        q_high <- unname(quantile(pt[idx], HIGH_PSEUDOTIME_QUANTILE, na.rm = TRUE))
        method <- "within_sample_quantiles"
      } else {
        q_low <- unname(quantile(inoa_pt, LOW_PSEUDOTIME_QUANTILE, na.rm = TRUE))
        q_high <- unname(quantile(inoa_pt, HIGH_PSEUDOTIME_QUANTILE, na.rm = TRUE))
        method <- "pooled_quantiles_fallback"
      }
      group[idx[pt[idx] <= q_low]] <- "iNOA immature-like Sertoli"
      group[idx[pt[idx] > q_low & pt[idx] < q_high]] <- "iNOA intermediate Sertoli"
      group[idx[pt[idx] >= q_high]] <- "iNOA mature-like Sertoli"
      threshold_rows[[length(threshold_rows) + 1]] <- data.frame(
        cell_type = SERTOLI_CELLTYPE,
        disease_group = "iNOA",
        sample_id = sid,
        grouping_method = method,
        low_pseudotime_quantile = LOW_PSEUDOTIME_QUANTILE,
        high_pseudotime_quantile = HIGH_PSEUDOTIME_QUANTILE,
        low_pseudotime_cutoff = q_low,
        high_pseudotime_cutoff = q_high,
        n_cells_used_for_thresholds = length(idx),
        stringsAsFactors = FALSE
      )
    }
  } else {
    q_low <- unname(quantile(inoa_pt, LOW_PSEUDOTIME_QUANTILE, na.rm = TRUE))
    q_high <- unname(quantile(inoa_pt, HIGH_PSEUDOTIME_QUANTILE, na.rm = TRUE))
    idx <- which(disease_group == "iNOA" & is.finite(pt))
    group[idx[pt[idx] <= q_low]] <- "iNOA immature-like Sertoli"
    group[idx[pt[idx] > q_low & pt[idx] < q_high]] <- "iNOA intermediate Sertoli"
    group[idx[pt[idx] >= q_high]] <- "iNOA mature-like Sertoli"
    threshold_rows[[1]] <- data.frame(
      cell_type = SERTOLI_CELLTYPE,
      disease_group = "iNOA",
      sample_id = "pooled_iNOA",
      grouping_method = "pooled_quantiles",
      low_pseudotime_quantile = LOW_PSEUDOTIME_QUANTILE,
      high_pseudotime_quantile = HIGH_PSEUDOTIME_QUANTILE,
      low_pseudotime_cutoff = q_low,
      high_pseudotime_cutoff = q_high,
      n_cells_used_for_thresholds = length(idx),
      stringsAsFactors = FALSE
    )
  }
  
  disease_sertoli$analysis_group <- group
  disease_sertoli$pseudotime_subgroup <- case_when(
    group == "iNOA immature-like Sertoli" ~ "iNOA low projected pseudotime",
    group == "iNOA intermediate Sertoli" ~ "iNOA intermediate projected pseudotime",
    group == "iNOA mature-like Sertoli" ~ "iNOA high projected pseudotime",
    group == "AZFa_Del Sertoli" ~ "AZFa_Del Sertoli",
    TRUE ~ NA_character_
  )
  
  thresholds <- bind_rows(threshold_rows)
  fwrite(thresholds, file.path(TABLE_DIR, "part5_iNOA_sertoli_pseudotime_subgroup_thresholds.csv"))
  
  disease_sertoli
}

build_combined_marker_expression <- function(ref, disease_sertoli, marker_genes) {
  normal_expr <- get_assay_data_safe(ref, ASSAY_USE, SLOT_USE)
  disease_expr <- get_assay_data_safe(disease_sertoli, ASSAY_USE, SLOT_USE)
  genes <- intersect(marker_genes, intersect(rownames(normal_expr), rownames(disease_expr)))
  if (length(genes) < 5) stop("Too few marker genes available in both normal and disease matrices.")
  
  normal_mat <- normal_expr[genes, , drop = FALSE]
  disease_mat <- disease_expr[genes, , drop = FALSE]
  colnames(normal_mat) <- paste0("Normal__", colnames(normal_mat))
  colnames(disease_mat) <- paste0("Disease__", colnames(disease_mat))
  
  expr <- cbind(normal_mat, disease_mat)
  
  normal_meta <- data.frame(
    analysis_cell = colnames(normal_mat),
    source = "Normal",
    disease_group = "Normal",
    sample_id = if ("sample_id" %in% colnames(ref@meta.data)) as.character(ref$sample_id) else "Normal",
    analysis_group = as.character(ref$analysis_group),
    projected_pseudotime = ref$pseudotime_final,
    stringsAsFactors = FALSE
  )
  
  disease_meta <- data.frame(
    analysis_cell = colnames(disease_mat),
    source = "Disease",
    disease_group = as.character(disease_sertoli$disease_group),
    sample_id = as.character(disease_sertoli$sample_id),
    analysis_group = as.character(disease_sertoli$analysis_group),
    projected_pseudotime = suppressWarnings(as.numeric(disease_sertoli$part4_projected_pseudotime)),
    stringsAsFactors = FALSE
  )
  
  meta <- bind_rows(normal_meta, disease_meta) %>%
    mutate(analysis_group = factor(analysis_group, levels = analysis_group_order))
  
  list(expr = expr, meta = meta, genes = genes)
}

compute_gene_set_scores <- function(expr, meta, gene_sets) {
  genes_use <- unique(unlist(gene_sets))
  genes_use <- intersect(genes_use, rownames(expr))
  mat <- as.matrix(expr[genes_use, , drop = FALSE])
  gene_sd <- apply(mat, 1, sd, na.rm = TRUE)
  keep <- is.finite(gene_sd) & gene_sd > 0
  mat <- mat[keep, , drop = FALSE]
  mat_z <- t(scale(t(mat)))
  mat_z[!is.finite(mat_z)] <- 0
  
  score_df <- meta
  for (nm in names(gene_sets)) {
    genes <- intersect(gene_sets[[nm]], rownames(mat_z))
    score_df[[nm]] <- if (length(genes) == 0) NA_real_ else colMeans(mat_z[genes, , drop = FALSE])
  }
  
  score_df %>%
    pivot_longer(
      cols = all_of(names(gene_sets)),
      names_to = "score_set",
      values_to = "score"
    )
}

summarise_scores <- function(score_long) {
  score_long %>%
    filter(!is.na(analysis_group), is.finite(score)) %>%
    group_by(score_set, analysis_group) %>%
    summarise(
      n_cells = n(),
      mean_score = mean(score, na.rm = TRUE),
      median_score = median(score, na.rm = TRUE),
      q25_score = quantile(score, 0.25, na.rm = TRUE),
      q75_score = quantile(score, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
}

summarise_scores_by_sample <- function(score_long) {
  score_long %>%
    filter(!is.na(analysis_group), is.finite(score)) %>%
    mutate(sample_id = ifelse(is.na(sample_id) | !nzchar(sample_id), disease_group, sample_id)) %>%
    group_by(score_set, source, disease_group, sample_id, analysis_group) %>%
    summarise(
      n_cells = n(),
      mean_score = mean(score, na.rm = TRUE),
      median_score = median(score, na.rm = TRUE),
      q25_score = quantile(score, 0.25, na.rm = TRUE),
      q75_score = quantile(score, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
}

summarise_sample_score_effects <- function(sample_summary, baseline_group = "Normal mature reference") {
  baseline <- sample_summary %>%
    filter(analysis_group == baseline_group) %>%
    group_by(score_set) %>%
    summarise(
      baseline_median_of_samples = median(median_score, na.rm = TRUE),
      baseline_n_samples = n_distinct(sample_id),
      .groups = "drop"
    )
  
  sample_summary %>%
    left_join(baseline, by = "score_set") %>%
    filter(analysis_group != baseline_group) %>%
    mutate(delta_vs_mature_reference_sample_median = median_score - baseline_median_of_samples) %>%
    group_by(score_set, analysis_group) %>%
    summarise(
      n_samples = n_distinct(sample_id),
      n_cells = sum(n_cells, na.rm = TRUE),
      median_of_sample_medians = median(median_score, na.rm = TRUE),
      q25_sample_median = quantile(median_score, 0.25, na.rm = TRUE),
      q75_sample_median = quantile(median_score, 0.75, na.rm = TRUE),
      baseline_median_of_samples = dplyr::first(baseline_median_of_samples),
      median_delta_vs_mature_reference = median(delta_vs_mature_reference_sample_median, na.rm = TRUE),
      fraction_samples_above_mature_reference =
        mean(delta_vs_mature_reference_sample_median > 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(score_set, analysis_group)
}

compare_score_groups <- function(score_long, baseline_group = "Normal mature reference") {
  rows <- list()
  for (ss in unique(score_long$score_set)) {
    baseline <- score_long %>%
      filter(.data$score_set == ss, .data$analysis_group == baseline_group) %>%
      pull(score)
    baseline <- baseline[is.finite(baseline)]
    
    for (grp in setdiff(unique(as.character(score_long$analysis_group)), baseline_group)) {
      vals <- score_long %>%
        filter(.data$score_set == ss, .data$analysis_group == grp) %>%
        pull(score)
      vals <- vals[is.finite(vals)]
      if (length(baseline) < 2 || length(vals) < 2) next
      p <- suppressWarnings(wilcox.test(vals, baseline)$p.value)
      rows[[length(rows) + 1]] <- data.frame(
        score_set = ss,
        baseline_group = baseline_group,
        comparison_group = grp,
        n_baseline = length(baseline),
        n_comparison = length(vals),
        median_baseline = median(baseline),
        median_comparison = median(vals),
        delta_median = median(vals) - median(baseline),
        p_value = p,
        stringsAsFactors = FALSE
      )
    }
  }
  
  bind_rows(rows) %>%
    mutate(p_adj_bh = p.adjust(p_value, method = "BH")) %>%
    arrange(score_set, p_adj_bh)
}

run_multicelltype_priority_extension <- function() {
  if (!file.exists(PART4_SAMPLE_SUMMARY)) {
    message("Part 5 extension skipped multicelltype prioritization: missing ", PART4_SAMPLE_SUMMARY)
    return(data.frame())
  }
  
  part4_sample <- fread(PART4_SAMPLE_SUMMARY)
  required_cols <- c("cell_type", "disease_group", "sample_id")
  deficit_col <- first_existing_col(
    c(
      "median_pseudotime_deficit_vs_mature",
      "developmental_age_deficit_vs_mature",
      "median_maturation_deficit_vs_mature"
    ),
    part4_sample
  )
  
  if (!all(required_cols %in% colnames(part4_sample)) || is.na(deficit_col)) {
    warning("Part 4 sample summary lacks columns needed for multicelltype prioritization.")
    return(data.frame())
  }
  
  priority_table <- part4_sample %>%
    mutate(
      disease_group = as.character(disease_group),
      sample_id = as.character(sample_id),
      deficit = to_num(.data[[deficit_col]])
    ) %>%
    filter(disease_group == "iNOA", is.finite(deficit)) %>%
    group_by(cell_type) %>%
    summarise(
      n_inoa_samples = n_distinct(sample_id),
      median_inoa_deficit = median(deficit, na.rm = TRUE),
      min_inoa_deficit = min(deficit, na.rm = TRUE),
      max_inoa_deficit = max(deficit, na.rm = TRUE),
      fraction_positive = mean(deficit > 0, na.rm = TRUE),
      priority_score = median_inoa_deficit * fraction_positive,
      .groups = "drop"
    ) %>%
    arrange(desc(priority_score), desc(median_inoa_deficit))
  
  fwrite(priority_table, file.path(TABLE_DIR, "part5_multicelltype_iNOA_developmental_deficit_priority.csv"))
  
  if (nrow(priority_table) > 0) {
    priority_plot_data <- priority_table %>%
      mutate(
        cell_type = factor(cell_type, levels = rev(cell_type)),
        is_sertoli = tolower(as.character(cell_type)) == tolower(SERTOLI_CELLTYPE)
      )
    
    p_priority <- ggplot(priority_plot_data, aes(x = cell_type, y = median_inoa_deficit, fill = is_sertoli)) +
      geom_col(width = 0.68, color = "white", linewidth = 0.22) +
      geom_linerange(aes(ymin = min_inoa_deficit, ymax = max_inoa_deficit),
                     linewidth = 0.38, color = "grey35") +
      geom_point(shape = 21, size = 1.7, color = "black", stroke = 0.22) +
      coord_flip() +
      scale_fill_manual(values = c("TRUE" = "#B0186B", "FALSE" = "#8CB7C9"), guide = "none") +
      labs(
        x = NULL,
        y = "Median iNOA developmental deficit",
        title = "Multicelltype screen supports Sertoli-focused mechanism analysis",
        subtitle = "Bars show iNOA sample median deficit; ticks show sample range"
      ) +
      theme_pub(10)
    
    save_supp(
      "FigS6G_multicelltype_iNOA_developmental_deficit_priority.pdf",
      p_priority,
      width = 7.0,
      height = max(4.2, 0.28 * nrow(priority_plot_data) + 2.0)
    )
  }
  
  priority_table
}

compare_gene_expression_matrix <- function(expr, cells1, cells2, group1_name, group2_name,
                                           genes = rownames(expr),
                                           min_pct = MIN_PCT_DEG,
                                           log2fc_cutoff = LOG2FC_CUTOFF,
                                           fdr_cutoff = FDR_CUTOFF,
                                           max_cells_per_group = Inf) {
  cells1 <- intersect(cells1, colnames(expr))
  cells2 <- intersect(cells2, colnames(expr))
  if (length(cells1) < MIN_CELLS_PER_GROUP || length(cells2) < MIN_CELLS_PER_GROUP) {
    warning("Skipping comparison with too few cells: ", group1_name, " vs ", group2_name)
    return(data.frame())
  }
  
  if (is.finite(max_cells_per_group)) {
    if (length(cells1) > max_cells_per_group) cells1 <- sample(cells1, max_cells_per_group)
    if (length(cells2) > max_cells_per_group) cells2 <- sample(cells2, max_cells_per_group)
  }
  
  genes <- intersect(genes, rownames(expr))
  mat1 <- expr[genes, cells1, drop = FALSE]
  mat2 <- expr[genes, cells2, drop = FALSE]
  
  mean1 <- Matrix::rowMeans(mat1)
  mean2 <- Matrix::rowMeans(mat2)
  pct1 <- Matrix::rowMeans(mat1 > 0) * 100
  pct2 <- Matrix::rowMeans(mat2 > 0) * 100
  keep <- (pct1 >= min_pct * 100) | (pct2 >= min_pct * 100)
  genes <- genes[keep]
  
  pvals <- vapply(genes, function(g) {
    x <- as.numeric(mat1[g, ])
    y <- as.numeric(mat2[g, ])
    tryCatch(
      suppressWarnings(wilcox.test(x, y)$p.value),
      error = function(e) NA_real_
    )
  }, numeric(1))
  
  out <- data.frame(
    gene = genes,
    group1_name = group1_name,
    group2_name = group2_name,
    n_group1 = length(cells1),
    n_group2 = length(cells2),
    mean_group1 = mean1[genes],
    mean_group2 = mean2[genes],
    pct_group1 = pct1[genes],
    pct_group2 = pct2[genes],
    log2fc = log2((mean2[genes] + 0.01) / (mean1[genes] + 0.01)),
    p_value = pvals,
    stringsAsFactors = FALSE
  )
  out$padj <- p.adjust(out$p_value, method = "BH")
  out$significant <- is.finite(out$padj) & out$padj < fdr_cutoff & abs(out$log2fc) >= log2fc_cutoff
  out[order(out$padj, -abs(out$log2fc)), , drop = FALSE]
}

summarise_de_sample_support <- function(expr, meta, genes) {
  empty_out <- data.frame(
    gene = character(),
    n_samples_with_both_states = integer(),
    median_sample_log2fc_mature_vs_immature = numeric(),
    fraction_samples_higher_in_immature_like = numeric(),
    fraction_samples_higher_in_mature_like = numeric(),
    stringsAsFactors = FALSE
  )
  genes <- intersect(genes, rownames(expr))
  meta <- meta[colnames(expr), , drop = FALSE]
  meta$analysis_cell <- rownames(meta)
  keep <- meta$disease_group == "iNOA" &
    meta$analysis_group %in% c("iNOA immature-like Sertoli", "iNOA mature-like Sertoli") &
    !is.na(meta$sample_id)
  meta <- meta[keep, , drop = FALSE]
  if (nrow(meta) == 0) return(empty_out)
  
  rows <- lapply(split(meta, meta$sample_id), function(mm) {
    imm <- mm$analysis_cell[mm$analysis_group == "iNOA immature-like Sertoli"]
    mat <- mm$analysis_cell[mm$analysis_group == "iNOA mature-like Sertoli"]
    if (length(imm) < MIN_CELLS_PER_GROUP || length(mat) < MIN_CELLS_PER_GROUP) return(NULL)
    mean_imm <- Matrix::rowMeans(expr[genes, imm, drop = FALSE])
    mean_mat <- Matrix::rowMeans(expr[genes, mat, drop = FALSE])
    data.frame(
      sample_id = unique(mm$sample_id)[1],
      gene = genes,
      n_immature_like_cells = length(imm),
      n_mature_like_cells = length(mat),
      sample_log2fc_mature_vs_immature = log2((mean_mat + 0.01) / (mean_imm + 0.01)),
      stringsAsFactors = FALSE
    )
  })
  rows <- bind_rows(rows)
  if (nrow(rows) == 0) return(empty_out)
  
  rows %>%
    group_by(gene) %>%
    summarise(
      n_samples_with_both_states = n_distinct(sample_id),
      median_sample_log2fc_mature_vs_immature = median(sample_log2fc_mature_vs_immature, na.rm = TRUE),
      fraction_samples_higher_in_immature_like =
        mean(sample_log2fc_mature_vs_immature < 0, na.rm = TRUE),
      fraction_samples_higher_in_mature_like =
        mean(sample_log2fc_mature_vs_immature > 0, na.rm = TRUE),
      .groups = "drop"
    )
}

prepare_representative_gene_data <- function(expr, meta, genes) {
  genes <- intersect(genes, rownames(expr))
  rows <- lapply(genes, function(g) {
    data.frame(
      analysis_cell = colnames(expr),
      gene = g,
      expression = as.numeric(expr[g, ]),
      stringsAsFactors = FALSE
    )
  })
  bind_rows(rows) %>%
    left_join(meta, by = "analysis_cell") %>%
    filter(!is.na(analysis_group))
}

annotate_gene_category <- function(gene_symbol) {
  case_when(
    gene_symbol %in% unlist(SENESCENCE_MARKERS) ~ "Senescence markers",
    gene_symbol %in% unlist(STRESS_MARKERS) ~ "Stress response",
    gene_symbol %in% unlist(FUNCTIONAL_GENES) ~ "Sertoli function",
    TRUE ~ "Other DEGs"
  )
}

run_go_enrichment <- function(genes, label) {
  if (!HAS_CLUSTERPROFILER || !HAS_ORGDB || !HAS_ENRICHPLOT) {
    message("Skipping enrichment for ", label, ": clusterProfiler/org.Hs.eg.db/enrichplot not all available.")
    return(NULL)
  }
  genes <- unique(genes[!is.na(genes)])
  if (length(genes) < 10) {
    message("Skipping enrichment for ", label, ": fewer than 10 genes.")
    return(NULL)
  }
  
  tryCatch({
    gene_map <- suppressMessages(clusterProfiler::bitr(
      genes,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db::org.Hs.eg.db
    ))
    entrez <- unique(gene_map$ENTREZID)
    if (length(entrez) < 10) return(NULL)
    
    ego <- suppressMessages(clusterProfiler::enrichGO(
      gene = entrez,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db,
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.10,
      readable = TRUE
    ))
    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(NULL)
    
    out_table <- as.data.frame(ego)
    fwrite(out_table, file.path(TABLE_DIR, paste0(label, "_GO_BP_enrichment.csv")))
    p <- enrichplot::dotplot(ego, showCategory = 15) +
      ggtitle(paste0(label, " GO biological process enrichment")) +
      theme_pub(9)
    save_supp(paste0("FigS6_GO_BP_", label, ".pdf"), p, width = 8.5, height = 6.5)
    ego
  }, error = function(e) {
    message("Enrichment failed for ", label, ": ", e$message)
    NULL
  })
}

# ------------------------------------------------------------------------------
# 2. Load data and define Sertoli analysis groups
# ------------------------------------------------------------------------------

message("\n========== Part 5 Sertoli molecular characterization ==========\n")

priority_table <- run_multicelltype_priority_extension()

sertoli_ref <- load_sertoli_reference()
sertoli_ref <- assign_normal_reference_groups(sertoli_ref)

projected_sertoli <- load_part4_sertoli_projection()
disease_obj <- combine_disease_objects(load_disease_objects())
disease_sertoli <- attach_projection_metadata(disease_obj, projected_sertoli)
disease_sertoli <- assign_disease_sertoli_groups(disease_sertoli)

cell_count_summary <- bind_rows(
  data.frame(
    source = "Normal reference",
    analysis_group = as.character(sertoli_ref$analysis_group),
    stringsAsFactors = FALSE
  ) %>%
    count(source, analysis_group, name = "n_cells"),
  data.frame(
    source = "Disease",
    analysis_group = as.character(disease_sertoli$analysis_group),
    disease_group = as.character(disease_sertoli$disease_group),
    stringsAsFactors = FALSE
  ) %>%
    count(source, disease_group, analysis_group, name = "n_cells")
)
fwrite(cell_count_summary, file.path(TABLE_DIR, "part5_sertoli_analysis_group_cell_counts.csv"))

sample_group_composition <- data.frame(
  disease_group = as.character(disease_sertoli$disease_group),
  sample_id = as.character(disease_sertoli$sample_id),
  analysis_group = as.character(disease_sertoli$analysis_group),
  projection_quality = as.character(disease_sertoli$part4_projection_quality),
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(analysis_group)) %>%
  count(disease_group, sample_id, analysis_group, projection_quality, name = "n_cells") %>%
  group_by(disease_group, sample_id) %>%
  mutate(fraction_within_sample = n_cells / sum(n_cells)) %>%
  ungroup()
fwrite(sample_group_composition, file.path(TABLE_DIR, "part5_sertoli_subgroup_sample_composition.csv"))

all_marker_genes <- unique(c(
  unlist(MAIN_GENE_SETS),
  unlist(EXTENDED_GENE_SETS),
  unlist(SUBCATEGORY_GENE_SETS),
  REPRESENTATIVE_GENES
))

combined_marker <- build_combined_marker_expression(
  sertoli_ref,
  disease_sertoli,
  all_marker_genes
)
marker_expr <- combined_marker$expr
marker_meta <- combined_marker$meta
fwrite(marker_meta, file.path(TABLE_DIR, "part5_sertoli_cell_metadata_for_marker_analysis.csv"))

# ------------------------------------------------------------------------------
# 3. Gene-set score validation
# ------------------------------------------------------------------------------

main_score_long <- compute_gene_set_scores(marker_expr, marker_meta, MAIN_GENE_SETS)
subcategory_score_long <- compute_gene_set_scores(marker_expr, marker_meta, SUBCATEGORY_GENE_SETS)

fwrite(main_score_long, file.path(TABLE_DIR, "part5_sertoli_main_gene_set_scores_by_cell.csv"))
fwrite(subcategory_score_long, file.path(TABLE_DIR, "part5_sertoli_subcategory_gene_set_scores_by_cell.csv"))

main_score_summary <- summarise_scores(main_score_long)
subcategory_score_summary <- summarise_scores(subcategory_score_long)
main_score_sample_summary <- summarise_scores_by_sample(main_score_long)
subcategory_score_sample_summary <- summarise_scores_by_sample(subcategory_score_long)
main_score_sample_effects <- summarise_sample_score_effects(main_score_sample_summary, "Normal mature reference")
score_tests_vs_mature <- compare_score_groups(main_score_long, "Normal mature reference")

fwrite(main_score_summary, file.path(TABLE_DIR, "part5_sertoli_main_gene_set_score_summary.csv"))
fwrite(subcategory_score_summary, file.path(TABLE_DIR, "part5_sertoli_subcategory_gene_set_score_summary.csv"))
fwrite(main_score_sample_summary, file.path(TABLE_DIR, "part5_sertoli_main_gene_set_score_summary_by_sample.csv"))
fwrite(subcategory_score_sample_summary, file.path(TABLE_DIR, "part5_sertoli_subcategory_gene_set_score_summary_by_sample.csv"))
fwrite(main_score_sample_effects, file.path(TABLE_DIR, "part5_sertoli_main_gene_set_sample_level_effects.csv"))
fwrite(score_tests_vs_mature, file.path(TABLE_DIR, "part5_sertoli_main_gene_set_score_tests_vs_mature_reference.csv"))

main_groups <- c(
  "Normal immature reference",
  "Normal mature reference",
  "AZFa_Del Sertoli",
  "iNOA immature-like Sertoli",
  "iNOA mature-like Sertoli"
)

program_summary_cell <- main_score_sample_effects %>%
  filter(analysis_group %in% setdiff(main_groups, "Normal mature reference")) %>%
  transmute(
    score_set,
    analysis_group,
    n_cells,
    n_samples,
    median_score = median_of_sample_medians,
    q25_score = q25_sample_median,
    q75_score = q75_sample_median,
    delta_vs_mature_reference = median_delta_vs_mature_reference,
    fraction_samples_above_mature_reference
  ) %>%
  bind_rows(
    main_score_sample_summary %>%
      filter(analysis_group == "Normal mature reference") %>%
      group_by(score_set, analysis_group) %>%
      summarise(
        n_cells = sum(n_cells),
        n_samples = n_distinct(sample_id),
        median_score = median(.data$median_score, na.rm = TRUE),
        q25_score = quantile(.data$median_score, 0.25, na.rm = TRUE),
        q75_score = quantile(.data$median_score, 0.75, na.rm = TRUE),
        delta_vs_mature_reference = 0,
        fraction_samples_above_mature_reference = NA_real_,
        .groups = "drop"
      )
  ) %>%
  group_by(score_set) %>%
  mutate(
    mature_reference_median = median_score[analysis_group == "Normal mature reference"][1]
  ) %>%
  ungroup() %>%
  mutate(
    analysis_group = factor(as.character(analysis_group), levels = main_groups),
    analysis_group_label = factor(short_group_labels[as.character(analysis_group)], levels = short_group_labels[main_groups]),
    score_set_label = factor(short_program_labels[as.character(score_set)], levels = rev(unname(short_program_labels)))
  )

p_scores <- ggplot(
  program_summary_cell,
  aes(x = analysis_group_label, y = score_set_label)
) +
  geom_tile(fill = "grey97", color = "white", linewidth = 0.55) +
  geom_point(
    aes(fill = median_score, size = abs(delta_vs_mature_reference)),
    shape = 21,
    color = "grey18",
    stroke = 0.28
  ) +
  scale_fill_gradient2(
    low = "#4C78A8",
    mid = "white",
    high = "#C44E52",
    midpoint = 0,
    name = "Median\nprogram score"
  ) +
  scale_size_area(max_size = 9, name = "|Delta vs\nmature ref.|") +
  labs(
    x = NULL,
    y = NULL,
    title = "Sample-level Sertoli programs support stress-biased rewiring"
  ) +
  theme_cell(9) +
  theme(
    axis.text.x = element_text(size = 8.5),
    axis.text.y = element_text(face = "bold"),
    legend.position = "right",
    panel.border = element_rect(color = "grey85", fill = NA, linewidth = 0.35)
  )

save_main("Fig6A_sertoli_gene_set_program_scores_sample_level.pdf", p_scores, width = 7.3, height = 3.2)

p_scores_violin <- main_score_long %>%
  filter(analysis_group %in% main_groups) %>%
  mutate(analysis_group = factor(as.character(analysis_group), levels = main_groups)) %>%
  ggplot(aes(x = analysis_group, y = score, fill = analysis_group)) +
  geom_violin(width = 0.78, alpha = 0.82, color = NA, trim = TRUE) +
  geom_boxplot(width = 0.16, outlier.size = 0.2, color = "grey20", fill = "white", linewidth = 0.25) +
  facet_wrap(~score_set, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = analysis_group_palette, guide = "none") +
  labs(
    x = NULL,
    y = "Relative gene-set score",
    title = "Sertoli molecular programs after disease projection"
  ) +
  theme_pub(9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

save_supp("FigS6D_sertoli_gene_set_program_scores_violin.pdf", p_scores_violin, width = 11, height = 4.6)

sample_program_delta_data <- main_score_sample_summary %>%
  group_by(score_set) %>%
  mutate(
    mature_reference_sample_median =
      median(median_score[analysis_group == "Normal mature reference"], na.rm = TRUE),
    delta_vs_mature_reference = median_score - mature_reference_sample_median
  ) %>%
  ungroup() %>%
  filter(
    disease_group == "iNOA",
    analysis_group %in% c("iNOA immature-like Sertoli", "iNOA mature-like Sertoli")
  ) %>%
  mutate(
    analysis_group = factor(as.character(analysis_group),
                            levels = c("iNOA immature-like Sertoli", "iNOA mature-like Sertoli")),
    analysis_group_label = factor(short_group_labels[as.character(analysis_group)],
                                  levels = short_group_labels[c("iNOA immature-like Sertoli", "iNOA mature-like Sertoli")])
  )
fwrite(sample_program_delta_data, file.path(TABLE_DIR, "part5_iNOA_sertoli_sample_program_deltas.csv"))

if (nrow(sample_program_delta_data) > 0) {
  p_sample_delta <- ggplot(
    sample_program_delta_data,
    aes(x = analysis_group_label, y = delta_vs_mature_reference, group = sample_id, color = sample_id)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.35) +
    geom_line(linewidth = 0.55, alpha = 0.82) +
    geom_point(aes(size = n_cells), alpha = 0.9) +
    facet_wrap(~score_set, nrow = 1, scales = "free_y") +
    scale_size_area(max_size = 5.2, name = "Cells") +
    labs(
      x = NULL,
      y = "Sample median score delta vs mature normal reference",
      color = "iNOA sample",
      title = "iNOA Sertoli molecular shifts are evaluated at sample level"
    ) +
    theme_cell(9) +
    theme(
      legend.position = "right",
      axis.text.x = element_text(size = 8.5)
    )
  save_main("Fig6D_iNOA_sertoli_sample_level_program_deltas.pdf", p_sample_delta, width = 7.6, height = 3.7)
}

sub_heatmap <- subcategory_score_summary %>%
  filter(analysis_group %in% main_groups) %>%
  group_by(score_set) %>%
  mutate(scaled_median_score = as.numeric(scale(median_score))) %>%
  ungroup() %>%
  mutate(
    analysis_group = factor(as.character(analysis_group), levels = main_groups),
    score_set = factor(score_set, levels = rev(unique(score_set)))
  )

p_sub_heatmap <- ggplot(sub_heatmap, aes(x = analysis_group, y = score_set, fill = scaled_median_score)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(low = "#4C78A8", mid = "white", high = "#C44E52", midpoint = 0, name = "Scaled\nmedian") +
  labs(
    x = NULL,
    y = NULL,
    title = "Sertoli marker subprogram score summary"
  ) +
  theme_pub(8) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid = element_blank()
  )

save_supp("FigS6A_sertoli_marker_subprogram_score_heatmap.pdf", p_sub_heatmap, width = 8.5, height = 6.8)

# ------------------------------------------------------------------------------
# 4. Folded Part 7 extension: expanded mechanism programs
# ------------------------------------------------------------------------------

extended_gene_set_availability <- bind_rows(lapply(names(EXTENDED_GENE_SETS), function(gs) {
  available <- intersect(EXTENDED_GENE_SETS[[gs]], rownames(marker_expr))
  data.frame(
    score_set = gs,
    n_defined = length(unique(EXTENDED_GENE_SETS[[gs]])),
    n_available = length(available),
    available_genes = paste(available, collapse = ";"),
    stringsAsFactors = FALSE
  )
}))
fwrite(extended_gene_set_availability, file.path(TABLE_DIR, "part5_extended_gene_set_availability.csv"))

extended_score_long <- compute_gene_set_scores(marker_expr, marker_meta, EXTENDED_GENE_SETS)
extended_score_summary <- summarise_scores(extended_score_long)
extended_score_sample_summary <- summarise_scores_by_sample(extended_score_long)
extended_score_sample_effects <- summarise_sample_score_effects(extended_score_sample_summary, "Normal mature reference")

fwrite(extended_score_long, file.path(TABLE_DIR, "part5_extended_gene_set_scores_by_cell.csv"))
fwrite(extended_score_summary, file.path(TABLE_DIR, "part5_extended_gene_set_score_summary.csv"))
fwrite(extended_score_sample_summary, file.path(TABLE_DIR, "part5_extended_gene_set_score_summary_by_sample.csv"))
fwrite(extended_score_sample_effects, file.path(TABLE_DIR, "part5_extended_gene_set_sample_level_effects.csv"))

extended_groups <- c(
  "Normal immature reference",
  "Normal mature reference",
  "AZFa_Del Sertoli",
  "iNOA immature-like Sertoli",
  "iNOA mature-like Sertoli"
)

extended_score_plot <- extended_score_sample_summary %>%
  filter(analysis_group %in% extended_groups) %>%
  mutate(
    analysis_group = factor(as.character(analysis_group), levels = extended_groups),
    analysis_group_label = factor(short_group_labels[as.character(analysis_group)], levels = short_group_labels[extended_groups]),
    score_set_label = factor(short_program_labels[as.character(score_set)], levels = rev(short_program_labels[names(EXTENDED_GENE_SETS)]))
  )

if (nrow(extended_score_plot) > 0) {
  p_extended_scores <- ggplot(extended_score_plot, aes(x = analysis_group_label, y = median_score, fill = analysis_group)) +
    geom_hline(yintercept = 0, color = "grey65", linewidth = 0.28, linetype = "dashed") +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.55, color = "grey35") +
    geom_point(position = position_jitter(width = 0.12, height = 0), shape = 21, size = 2.1, color = "black", stroke = 0.25) +
    facet_wrap(~ score_set_label, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = analysis_group_palette, guide = "none") +
    labs(
      x = NULL,
      y = "Sample-median program score",
      title = "Extended Sertoli mechanism programs support the projected-state pattern",
      subtitle = "Folded-in Part 7 extension; sample medians reduce cell-level pseudoreplication"
    ) +
    theme_pub(8.5) +
    theme(axis.text.x = element_text(size = 7.8), strip.text = element_text(size = 8.2))
  
  save_supp("FigS6H_sertoli_extended_mechanism_program_scores.pdf", p_extended_scores, width = 9.4, height = 5.8)
}

extended_delta <- extended_score_sample_summary %>%
  filter(analysis_group %in% extended_groups) %>%
  group_by(score_set, analysis_group) %>%
  summarise(
    n_samples = n_distinct(sample_id),
    n_cells = sum(n_cells, na.rm = TRUE),
    median_of_sample_medians = median(median_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(score_set) %>%
  mutate(
    mature_reference_median = median_of_sample_medians[analysis_group == "Normal mature reference"][1],
    delta_vs_mature_reference = median_of_sample_medians - mature_reference_median
  ) %>%
  ungroup()
fwrite(extended_delta, file.path(TABLE_DIR, "part5_extended_gene_set_delta_vs_mature_reference.csv"))

extended_delta_plot <- extended_delta %>%
  mutate(
    analysis_group = factor(as.character(analysis_group), levels = extended_groups),
    analysis_group_label = factor(short_group_labels[as.character(analysis_group)], levels = short_group_labels[extended_groups]),
    score_set_label = factor(short_program_labels[as.character(score_set)], levels = rev(short_program_labels[names(EXTENDED_GENE_SETS)]))
  )

if (nrow(extended_delta_plot) > 0) {
  p_extended_delta <- ggplot(extended_delta_plot, aes(x = analysis_group_label, y = score_set_label, fill = delta_vs_mature_reference)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.2f", delta_vs_mature_reference)), size = 2.5, color = "black") +
    scale_fill_gradient2(low = "#4C78A8", mid = "white", high = "#C44E52", midpoint = 0, name = "Delta vs\nmature ref.") +
    labs(
      x = NULL,
      y = NULL,
      title = "Expanded mechanisms are summarized relative to mature normal Sertoli cells"
    ) +
    theme_pub(8.5) +
    theme(axis.text.x = element_text(size = 7.8), axis.text.y = element_text(face = "bold"))
  
  save_supp("FigS6I_sertoli_extended_program_delta_heatmap.pdf", p_extended_delta, width = 8.6, height = 4.4)
}

# ------------------------------------------------------------------------------
# 5. Representative marker expression and marker-gene comparisons
# ------------------------------------------------------------------------------

representative_data <- prepare_representative_gene_data(marker_expr, marker_meta, REPRESENTATIVE_GENES)
fwrite(representative_data, file.path(TABLE_DIR, "part5_sertoli_representative_gene_expression_by_cell.csv"))

representative_gene_order <- c(
  "CDKN1A", "CDKN2A", "HSPA5", "HMOX1", "CXCL8", "IL6",
  "AMH", "CLDN11", "KITLG", "FSHR", "SOX9", "GDNF"
)
representative_gene_order <- intersect(representative_gene_order, unique(representative_data$gene))

gene_dot_data <- representative_data %>%
  filter(analysis_group %in% main_groups, gene %in% representative_gene_order) %>%
  group_by(gene, analysis_group) %>%
  summarise(
    mean_expression = mean(expression, na.rm = TRUE),
    pct_expressing = mean(expression > 0, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  group_by(gene) %>%
  mutate(scaled_mean_expression = as.numeric(scale(mean_expression))) %>%
  ungroup() %>%
  mutate(
    analysis_group = factor(as.character(analysis_group), levels = main_groups),
    analysis_group_label = factor(short_group_labels[as.character(analysis_group)], levels = short_group_labels[main_groups]),
    gene = factor(gene, levels = rev(representative_gene_order))
  )

p_genes <- ggplot(gene_dot_data, aes(x = analysis_group_label, y = gene)) +
  geom_point(
    aes(size = pct_expressing, fill = scaled_mean_expression),
    shape = 21,
    color = "grey18",
    stroke = 0.25
  ) +
  scale_fill_gradient2(
    low = "#4C78A8",
    mid = "white",
    high = "#C44E52",
    midpoint = 0,
    name = "Scaled mean\nexpression"
  ) +
  scale_size_area(max_size = 7.5, limits = c(0, 100), name = "Expressing\ncells (%)") +
  labs(
    x = NULL,
    y = NULL,
    title = "Representative Sertoli markers support an immature stress state"
  ) +
  theme_cell(9) +
  theme(
    axis.text.x = element_text(size = 8.5),
    axis.text.y = element_text(face = "bold.italic"),
    legend.position = "right",
    panel.border = element_rect(color = "grey85", fill = NA, linewidth = 0.35)
  )

save_main("Fig6B_sertoli_representative_gene_expression.pdf", p_genes, width = 7.8, height = 5.2)

p_genes_violin <- representative_data %>%
  filter(analysis_group %in% main_groups) %>%
  mutate(analysis_group = factor(as.character(analysis_group), levels = main_groups)) %>%
  ggplot(aes(x = analysis_group, y = expression, fill = analysis_group)) +
  geom_violin(width = 0.78, alpha = 0.78, color = NA, trim = TRUE) +
  geom_boxplot(width = 0.14, outlier.size = 0.15, color = "grey20", fill = "white", linewidth = 0.22) +
  facet_wrap(~gene, ncol = 4, scales = "free_y") +
  scale_fill_manual(values = analysis_group_palette, guide = "none") +
  labs(
    x = NULL,
    y = "Log-normalized expression",
    title = "Representative Sertoli stress and function genes"
  ) +
  theme_pub(8) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

save_supp("FigS6E_sertoli_representative_gene_expression_violin.pdf", p_genes_violin, width = 11.5, height = 7.2)

marker_cells <- split(marker_meta$analysis_cell, as.character(marker_meta$analysis_group))
marker_comparisons <- bind_rows(
  compare_gene_expression_matrix(
    marker_expr,
    marker_cells[["Normal mature reference"]],
    marker_cells[["AZFa_Del Sertoli"]],
    "Normal mature reference",
    "AZFa_Del Sertoli",
    genes = combined_marker$genes,
    max_cells_per_group = Inf
  ) %>% mutate(comparison = "AZFa_Del_vs_Normal_mature"),
  compare_gene_expression_matrix(
    marker_expr,
    marker_cells[["Normal mature reference"]],
    marker_cells[["iNOA immature-like Sertoli"]],
    "Normal mature reference",
    "iNOA immature-like Sertoli",
    genes = combined_marker$genes,
    max_cells_per_group = Inf
  ) %>% mutate(comparison = "iNOA_immature_like_vs_Normal_mature"),
  compare_gene_expression_matrix(
    marker_expr,
    marker_cells[["Normal mature reference"]],
    marker_cells[["iNOA mature-like Sertoli"]],
    "Normal mature reference",
    "iNOA mature-like Sertoli",
    genes = combined_marker$genes,
    max_cells_per_group = Inf
  ) %>% mutate(comparison = "iNOA_mature_like_vs_Normal_mature")
) %>%
  mutate(biological_category = annotate_gene_category(gene)) %>%
  arrange(comparison, padj)

fwrite(marker_comparisons, file.path(TABLE_DIR, "part5_sertoli_marker_gene_comparisons.csv"))

# ------------------------------------------------------------------------------
# 6. iNOA immature-like vs mature-like differential expression
# ------------------------------------------------------------------------------

disease_expr <- get_assay_data_safe(disease_sertoli, ASSAY_USE, SLOT_USE)
disease_meta <- disease_sertoli@meta.data
inoa_immature_cells <- rownames(disease_meta)[disease_meta$analysis_group == "iNOA immature-like Sertoli"]
inoa_mature_cells <- rownames(disease_meta)[disease_meta$analysis_group == "iNOA mature-like Sertoli"]

deg_inoa <- compare_gene_expression_matrix(
  disease_expr,
  inoa_immature_cells,
  inoa_mature_cells,
  "iNOA immature-like Sertoli",
  "iNOA mature-like Sertoli",
  genes = rownames(disease_expr),
  max_cells_per_group = MAX_CELLS_PER_GROUP_FOR_DEG
) %>%
  mutate(
    biological_category = annotate_gene_category(gene),
    cell_level_direction = case_when(
      significant & log2fc < 0 ~ "Higher in iNOA immature-like Sertoli",
      significant & log2fc > 0 ~ "Higher in iNOA mature-like Sertoli",
      TRUE ~ "Not significant"
    )
  )

deg_sample_support <- summarise_de_sample_support(
  disease_expr,
  disease_meta,
  genes = deg_inoa$gene
)
fwrite(deg_sample_support, file.path(TABLE_DIR, "part5_iNOA_DEG_sample_level_support.csv"))

deg_inoa <- deg_inoa %>%
  left_join(deg_sample_support, by = "gene") %>%
  mutate(
    sample_supported = case_when(
      cell_level_direction == "Higher in iNOA immature-like Sertoli" &
        fraction_samples_higher_in_immature_like >= SAMPLE_CONSISTENCY_MIN_FRACTION ~ TRUE,
      cell_level_direction == "Higher in iNOA mature-like Sertoli" &
        fraction_samples_higher_in_mature_like >= SAMPLE_CONSISTENCY_MIN_FRACTION ~ TRUE,
      TRUE ~ FALSE
    ),
    direction = case_when(
      significant & sample_supported & log2fc < 0 ~ "Higher in iNOA immature-like Sertoli",
      significant & sample_supported & log2fc > 0 ~ "Higher in iNOA mature-like Sertoli",
      significant & !sample_supported ~ "Cell-level only",
      TRUE ~ "Not significant"
    )
  )

fwrite(deg_inoa, file.path(TABLE_DIR, "part5_iNOA_immature_like_vs_mature_like_DEG.csv"))

top_deg <- bind_rows(
  deg_inoa %>%
    filter(significant, sample_supported, log2fc < 0) %>%
    arrange(log2fc) %>%
    head(20) %>%
    mutate(direction = "Higher in iNOA immature-like Sertoli"),
  deg_inoa %>%
    filter(significant, sample_supported, log2fc > 0) %>%
    arrange(desc(log2fc)) %>%
    head(20) %>%
    mutate(direction = "Higher in iNOA mature-like Sertoli")
) %>%
  arrange(log2fc) %>%
  mutate(gene = factor(gene, levels = gene))

fwrite(top_deg, file.path(TABLE_DIR, "part5_iNOA_subgroup_top_DEGs.csv"))

volcano_label_genes <- unique(c(
  intersect(VOLCANO_IMMATURE_LABEL_GENES, deg_inoa$gene),
  intersect(VOLCANO_MATURE_LABEL_GENES, deg_inoa$gene)
))

volcano_data <- deg_inoa %>%
  mutate(
    neg_log10_padj = -log10(padj + 1e-300),
    plot_neg_log10_padj = pmin(neg_log10_padj, VOLCANO_Y_CAP),
    volcano_label_group = case_when(
      gene %in% VOLCANO_IMMATURE_LABEL_GENES & significant & sample_supported & log2fc < -LOG2FC_CUTOFF ~ "Immature-like key genes",
      gene %in% VOLCANO_MATURE_LABEL_GENES & significant & sample_supported & log2fc > LOG2FC_CUTOFF ~ "Mature-like key genes",
      TRUE ~ NA_character_
    ),
    label = ifelse(!is.na(volcano_label_group) & gene %in% volcano_label_genes, gene, NA_character_)
  )

volcano_colors <- c(
  "Higher in iNOA immature-like Sertoli" = "#B53B6F",
  "Higher in iNOA mature-like Sertoli" = "#4E79A7",
  "Cell-level only" = "grey55",
  "Not significant" = "grey82"
)

label_x <- volcano_data$log2fc[!is.na(volcano_data$label)]
volcano_x_limits <- c(
  min(quantile(volcano_data$log2fc, 0.002, na.rm = TRUE), label_x, na.rm = TRUE) - 0.25,
  max(quantile(volcano_data$log2fc, 0.998, na.rm = TRUE), label_x, na.rm = TRUE) + 0.25
)

p_volcano <- ggplot(volcano_data, aes(x = log2fc, y = plot_neg_log10_padj, color = direction)) +
  geom_point(data = volcano_data %>% filter(direction == "Not significant"), size = 0.45, alpha = 0.28) +
  geom_point(data = volcano_data %>% filter(direction != "Not significant"), size = 1.15, alpha = 0.78) +
  geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = "dashed", color = "grey45", linewidth = 0.35) +
  geom_vline(xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF), linetype = "dashed", color = "grey45", linewidth = 0.35) +
  scale_color_manual(values = volcano_colors, name = NULL) +
  coord_cartesian(xlim = volcano_x_limits, ylim = c(0, VOLCANO_Y_CAP), clip = "off") +
  labs(
    x = "log2 fold change (mature-like / immature-like)",
    y = expression(-log[10](FDR)),
    title = "iNOA Sertoli immature-like cells are enriched for stress-response genes"
  ) +
  theme_cell(9) +
  theme(
    legend.position = "bottom",
    legend.justification = "left",
    panel.grid.major.y = element_line(color = "grey93", linewidth = 0.25)
  )

if (HAS_GGREPEL && sum(!is.na(volcano_data$label)) > 0) {
  p_volcano <- p_volcano +
    ggrepel::geom_label_repel(
      data = volcano_data %>% filter(!is.na(label), volcano_label_group == "Immature-like key genes"),
      aes(label = label),
      size = 2.3,
      max.overlaps = Inf,
      box.padding = 0.3,
      point.padding = 0.2,
      segment.color = "grey55",
      segment.size = 0.2,
      force = 2.5,
      min.segment.length = 0,
      na.rm = TRUE,
      show.legend = FALSE
    ) +
    ggrepel::geom_label_repel(
      data = volcano_data %>% filter(!is.na(label), volcano_label_group == "Mature-like key genes"),
      aes(label = label),
      size = 2.3,
      max.overlaps = Inf,
      box.padding = 0.3,
      point.padding = 0.2,
      segment.color = "grey55",
      segment.size = 0.2,
      force = 2.5,
      min.segment.length = 0,
      na.rm = TRUE,
      show.legend = FALSE
    )
}

save_main("Fig6C_iNOA_immature_like_vs_mature_like_volcano.pdf", p_volcano, width = 7.2, height = 5.6)

if (nrow(top_deg) > 0) {
  p_top_deg <- ggplot(top_deg, aes(x = gene, y = log2fc, fill = direction)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.2) +
    geom_hline(yintercept = 0, color = "grey30", linewidth = 0.35) +
    coord_flip() +
    scale_fill_manual(values = volcano_colors[names(volcano_colors) != "Not significant"], name = NULL) +
    labs(
      x = NULL,
      y = "log2 fold change (mature-like / immature-like)",
      title = "Top DEGs separating iNOA Sertoli projected states"
    ) +
    theme_pub(9) +
    theme(legend.position = "bottom")
  
  save_supp("FigS6B_iNOA_subgroup_top_DEG_barplot.pdf", p_top_deg, width = 8.2, height = 8.5)
}

category_bar_data <- deg_inoa %>%
  filter(significant, sample_supported) %>%
  mutate(
    biological_category = factor(
      biological_category,
      levels = c("Senescence markers", "Stress response", "Sertoli function", "Other DEGs")
    ),
    direction = factor(
      direction,
      levels = c("Higher in iNOA immature-like Sertoli", "Higher in iNOA mature-like Sertoli")
    )
  ) %>%
  count(biological_category, direction, name = "n_genes") %>%
  mutate(n_plot = ifelse(direction == "Higher in iNOA immature-like Sertoli", -n_genes, n_genes))

fwrite(category_bar_data, file.path(TABLE_DIR, "part5_iNOA_DEG_biological_category_counts.csv"))

if (nrow(category_bar_data) > 0) {
  p_cat <- ggplot(category_bar_data, aes(x = biological_category, y = n_plot, fill = direction)) +
    geom_col(width = 0.65, alpha = 0.9) +
    geom_hline(yintercept = 0, color = "grey25", linewidth = 0.4) +
    geom_text(aes(label = abs(n_genes), vjust = ifelse(n_plot >= 0, -0.35, 1.15)), size = 3, fontface = "bold") +
    coord_flip() +
    scale_fill_manual(
      values = volcano_colors[c("Higher in iNOA immature-like Sertoli", "Higher in iNOA mature-like Sertoli")],
      name = NULL
    ) +
    scale_y_continuous(labels = function(x) abs(x)) +
    labs(
      x = NULL,
      y = "Number of DEGs",
      title = "Biological categories of iNOA Sertoli state-associated DEGs"
    ) +
    theme_pub(10) +
    theme(legend.position = "bottom")
  
  save_supp("FigS6C_iNOA_DEG_biological_category_barplot.pdf", p_cat, width = 8, height = 5.2)
}

# ------------------------------------------------------------------------------
# 7. Optional enrichment analysis
# ------------------------------------------------------------------------------

genes_higher_immature <- deg_inoa %>%
  filter(significant, sample_supported, log2fc < -LOG2FC_CUTOFF) %>%
  pull(gene)

genes_higher_mature <- deg_inoa %>%
  filter(significant, sample_supported, log2fc > LOG2FC_CUTOFF) %>%
  pull(gene)

ego_immature <- run_go_enrichment(genes_higher_immature, "iNOA_immature_like_up")
ego_mature <- run_go_enrichment(genes_higher_mature, "iNOA_mature_like_up")

if (!is.null(ego_immature)) {
  go_main_data <- as.data.frame(ego_immature) %>%
    mutate(
      gene_ratio_num = parse_gene_ratio(GeneRatio),
      neg_log10_padj = -log10(p.adjust + 1e-300)
    ) %>%
    arrange(p.adjust) %>%
    head(10) %>%
    arrange(gene_ratio_num) %>%
    mutate(Description_wrapped = factor(wrap_text(Description, width = 32), levels = wrap_text(Description, width = 32)))
  
  p_go_main <- ggplot(go_main_data, aes(x = gene_ratio_num, y = Description_wrapped)) +
    geom_segment(aes(x = 0, xend = gene_ratio_num, yend = Description_wrapped), color = "grey82", linewidth = 0.42) +
    geom_point(
      aes(size = Count, fill = neg_log10_padj),
      shape = 21,
      color = "grey18",
      stroke = 0.28
    ) +
    scale_x_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.08))) +
    scale_size_area(max_size = 8, name = "Gene count") +
    scale_fill_gradient(low = "#4E79A7", high = "#C44E52", name = "-log10(FDR)") +
    labs(
      x = "Gene ratio",
      y = NULL,
      title = "Immature-like iNOA Sertoli cells activate ER stress and glycoprotein programs"
    ) +
    theme_cell(9) +
    theme(
      axis.text.y = element_text(size = 8.2),
      legend.position = "right",
      panel.grid.major.x = element_line(color = "grey93", linewidth = 0.25)
    )
  
  save_supp("FigS6F_iNOA_immature_like_GO_BP_enrichment.pdf", p_go_main, width = 7.4, height = 4.9)
}

# ------------------------------------------------------------------------------
# 8. Figure/table placement plan
# ------------------------------------------------------------------------------

placement <- data.frame(
  placement = c(
    "Main Figure 6A",
    "Main Figure 6B",
    "Main Figure 6C",
    "Main Figure 6D",
    "Supplementary Figure S6A",
    "Supplementary Figure S6B",
    "Supplementary Figure S6C",
    "Supplementary Figure S6D",
    "Supplementary Figure S6E",
    "Supplementary Figure S6F",
    "Supplementary Figure S6G",
    "Supplementary Figure S6H",
    "Supplementary Figure S6I",
    "Supplementary Table S25",
    "Supplementary Table S26",
    "Supplementary Table S27",
    "Supplementary Table S28",
    "Supplementary Table S29",
    "Supplementary Table S30",
    "Supplementary Table S31",
    "Supplementary Table S32",
    "Supplementary Table S33",
    "Supplementary Table S34",
    "Supplementary Table S35",
    "Supplementary Table S36",
    "Supplementary Table S37",
    "Supplementary Table S38",
    "Supplementary Table S39"
  ),
  file = c(
    "figures/main/Fig6A_sertoli_gene_set_program_scores_sample_level.pdf",
    "figures/main/Fig6B_sertoli_representative_gene_expression.pdf",
    "figures/main/Fig6C_iNOA_immature_like_vs_mature_like_volcano.pdf",
    "figures/main/Fig6D_iNOA_sertoli_sample_level_program_deltas.pdf",
    "figures/supplementary/FigS6A_sertoli_marker_subprogram_score_heatmap.pdf",
    "figures/supplementary/FigS6B_iNOA_subgroup_top_DEG_barplot.pdf",
    "figures/supplementary/FigS6C_iNOA_DEG_biological_category_barplot.pdf",
    "figures/supplementary/FigS6D_sertoli_gene_set_program_scores_violin.pdf",
    "figures/supplementary/FigS6E_sertoli_representative_gene_expression_violin.pdf",
    "figures/supplementary/FigS6F_iNOA_immature_like_GO_BP_enrichment.pdf",
    "figures/supplementary/FigS6G_multicelltype_iNOA_developmental_deficit_priority.pdf",
    "figures/supplementary/FigS6H_sertoli_extended_mechanism_program_scores.pdf",
    "figures/supplementary/FigS6I_sertoli_extended_program_delta_heatmap.pdf",
    "tables/part5_sertoli_analysis_group_cell_counts.csv",
    "tables/part5_sertoli_projection_quality_filter_summary.csv",
    "tables/part5_sertoli_subgroup_sample_composition.csv",
    "tables/part5_sertoli_main_gene_set_score_summary.csv",
    "tables/part5_sertoli_main_gene_set_score_summary_by_sample.csv",
    "tables/part5_sertoli_main_gene_set_sample_level_effects.csv",
    "tables/part5_sertoli_main_gene_set_score_tests_vs_mature_reference.csv",
    "tables/part5_sertoli_marker_gene_comparisons.csv",
    "tables/part5_iNOA_immature_like_vs_mature_like_DEG.csv",
    "tables/part5_iNOA_DEG_sample_level_support.csv",
    "tables/*_GO_BP_enrichment.csv",
    "tables/part5_multicelltype_iNOA_developmental_deficit_priority.csv",
    "tables/part5_extended_gene_set_availability.csv",
    "tables/part5_extended_gene_set_score_summary_by_sample.csv",
    "tables/part5_extended_gene_set_delta_vs_mature_reference.csv"
  ),
  rationale = c(
    "Show sample-level stress/senescence/function program changes after projected Sertoli state assignment.",
    "Show representative marker genes supporting the gene-set-level pattern.",
    "Identify genes distinguishing iNOA immature-like and mature-like projected Sertoli states, highlighting sample-supported signals.",
    "Show whether iNOA sample-level program deltas are directionally consistent across projected Sertoli states.",
    "Show finer subprogram-level marker score patterns.",
    "Show top DEGs driving iNOA projected-state separation.",
    "Classify iNOA state-associated DEGs into senescence, stress, Sertoli function, and other groups.",
    "Show cell-level distributions supporting the compact main gene-set summary.",
    "Show cell-level representative gene distributions supporting the compact main dot plot.",
    "Optional GO enrichment among sample-supported genes higher in iNOA immature-like Sertoli.",
    "Folded Part 7 cell-type prioritization showing why Sertoli is a reasonable mechanism entry point.",
    "Folded Part 7 expanded mechanism programs, summarized by sample medians.",
    "Folded Part 7 expanded mechanism-program deltas relative to mature normal Sertoli reference.",
    "Cell counts for all Part 5 comparison groups.",
    "Document how Part 4 projection-quality filtering affected Sertoli cells retained for Part 5.",
    "Document sample composition of iNOA projected-state subgroups to control for sample imbalance.",
    "Gene-set score group summaries.",
    "Sample-level gene-set score summaries.",
    "Sample-level gene-set effect summaries relative to mature normal reference.",
    "Statistical tests for gene-set scores versus mature normal reference.",
    "Marker-gene-level disease versus normal comparisons.",
    "All-gene DE analysis for iNOA immature-like versus mature-like Sertoli states with sample-support annotations.",
    "Sample-level support statistics for iNOA immature-like versus mature-like DEGs.",
    "Optional GO enrichment result tables.",
    "Cell-type prioritization table from Part 4 sample-level projection summaries.",
    "Availability check for genes used in the folded Part 7 expanded mechanism programs.",
    "Sample-level expanded mechanism program score summaries.",
    "Expanded mechanism program deltas versus mature normal Sertoli reference."
  ),
  stringsAsFactors = FALSE
)

fwrite(placement, file.path(TABLE_DIR, "part5_sertoli_mechanism_figure_table_plan.csv"))

message("Part 5 Sertoli mechanism characterization finished: ", normalizePath(OUTDIR, winslash = "/"))

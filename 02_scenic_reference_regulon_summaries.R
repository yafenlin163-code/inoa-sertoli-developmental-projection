#!/usr/bin/env Rscript

# ==============================================================================
# 02_scenic_reference_regulon_summaries.R
#
# Purpose:
#   Summarize completed SCENIC/regulon outputs for normal postnatal testicular
#   cell-type references and generate manuscript source tables and figures.
#
# Scope:
#   This script reads precomputed SCENIC outputs. It does not run pySCENIC/
#   SCENIC regulon inference from expression matrices. The outputs are used as
#   supportive regulatory-activity evidence and should not be interpreted as
#   causal proof of Sertoli-cell developmental arrest.
#
# Expected input directory structure:
#   <SCENIC_RESULTS_DIR>/
#     scenic_sertoli_results/
#     scenic_leydig_results/
#     scenic_pmc_results/
#     scenic_endothelial_results/
#     scenic_macrophage_results/
#     scenic_t_cell_results/
#     scenic_spermatogonia_results/
#     scenic_spermatocyte_results/
#     scenic_spermatid_results/
#
# Each directory may contain:
#   age_correlation_regulons.csv
#   cluster_specificity_regulons.csv
#   GO_enrichment_summary_top5.csv
#   intermediate/age_mean_auc.rds
#   intermediate/cluster_auc.rds
#   intermediate/regulons_clean.rds
#   regulon_source_detailed.csv
#
# Outputs:
#   results/02_scenic_reference_regulon_summaries/
#     Fig2A_SCENIC_QC_summary.pdf
#     Fig2B_age_regulon_landscape.pdf
#     Fig2C_age_dynamic_regulon_heatmap.pdf
#     Fig2D_GO_bubble_age_dynamic_regulons.pdf
#     Fig2E_key_TF_age_trajectories.pdf
#     FigS2_GO_bubble_age_dynamic_regulons_full.pdf
#     FigS2_cluster_marker_regulon_bubble.pdf
#     tables/*.csv
#
# Default paths can be overridden with environment variables:
#   PROJECT_DIR, SCENIC_RESULTS_DIR, OUTPUT_DIR
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(patchwork)
  library(stringr)
  library(forcats)
  library(scales)
  library(RColorBrewer)
})

# ------------------------------------------------------------------------------
# 0. User configuration
# ------------------------------------------------------------------------------

set.seed(42)

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = getwd())
BASE_DIR <- Sys.getenv(
  "SCENIC_RESULTS_DIR",
  unset = file.path(PROJECT_DIR, "results", "scenic")
)
OUT_DIR <- Sys.getenv(
  "OUTPUT_DIR",
  unset = file.path(PROJECT_DIR, "results", "02_scenic_reference_regulon_summaries")
)
TABLE_DIR <- file.path(OUT_DIR, "tables")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)

celltype_order <- c(
  "Sertoli", "Leydig", "PMC", "Endothelial", "Macrophage",
  "T_cell", "Spermatogonia", "Spermatocyte", "Spermatid"
)

celltype_dirs <- c(
  Sertoli       = "scenic_sertoli_results",
  Leydig        = "scenic_leydig_results",
  PMC           = "scenic_pmc_results",
  Endothelial   = "scenic_endothelial_results",
  Macrophage    = "scenic_macrophage_results",
  T_cell        = "scenic_t_cell_results",
  Spermatogonia = "scenic_spermatogonia_results",
  Spermatocyte  = "scenic_spermatocyte_results",
  Spermatid     = "scenic_spermatid_results"
)

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

# Direction palette inspired by the user's reference colors.
direction_palette <- c(
  "Up with age" = "#BB6A8A",
  "Down with age" = "#91B0CB",
  "Stable" = "grey78"
)

direction_annotation_palette <- c(
  "Up with age" = "#8F3F63",
  "Down with age" = "#4F7FA3",
  "Stable" = "#C9C9C9"
)

key_tf_map <- list(
  Sertoli       = c("SOX9", "WT1", "GATA4", "DMRT1", "AR", "NR5A1", "FOXO1"),
  Leydig        = c("NR5A1", "GATA4", "DLX5", "TCF21", "AR"),
  PMC           = c("MYOCD", "SRF", "TEAD1", "TEAD3", "KLF4"),
  Endothelial   = c("ERG", "FLI1", "KLF2", "KLF4", "SOX17"),
  Macrophage    = c("SPI1", "CEBPB", "MAFB", "IRF8", "NFKB1"),
  T_cell        = c("TCF7", "LEF1", "RUNX3", "TBX21", "GATA3"),
  Spermatogonia = c("DMRT1", "SOHLH1", "SOHLH2", "ID4", "ZBTB16"),
  Spermatocyte  = c("MEIOSIN", "STRA8", "DMRT1", "MYBL1", "TCFL5"),
  Spermatid     = c("CREM", "RFX2", "RFX3", "YBX2", "SPZ1")
)

# Biological terms used for GO display.
# The full enrichment table is still exported; these rules only control the
# concise main-figure version of the GO bubble plot.
go_main_include_patterns <- c(
  "differentiation", "maturation", "spermatogenesis", "meiosis", "germ cell",
  "male sex", "sex differentiation", "steroid", "androgen", "cholesterol",
  "cell adhesion", "cell-cell adhesion", "cell-substrate", "junction",
  "blood-testis", "barrier", "extracellular matrix", "endothelium",
  "angiogenesis", "immune", "leukocyte", "lymphocyte", "T cell",
  "interferon", "antiviral", "inflammatory", "cell cycle", "mitotic",
  "checkpoint", "apoptotic"
)

go_main_exclude_patterns <- c(
  "forebrain", "brain", "neuron", "neural", "kidney", "renal",
  "roof of mouth", "mouth", "respiratory", "lung", "skeletal muscle",
  "striated muscle", "muscle tissue", "osteoblast", "bone"
)

# ------------------------------------------------------------------------------
# 1. Helper functions
# ------------------------------------------------------------------------------

theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 2),
      plot.subtitle = element_text(color = "grey35", hjust = 0, size = base_size),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.35),
      legend.title = element_text(face = "bold"),
      legend.key.size = unit(0.45, "cm"),
      strip.background = element_rect(fill = "grey94", color = NA),
      strip.text = element_text(face = "bold", color = "black")
    )
}

safe_read_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  readr::read_csv(path, show_col_types = FALSE)
}

safe_read_rds <- function(path) {
  if (!file.exists(path)) return(NULL)
  readRDS(path)
}

bind_or_empty <- function(rows, columns) {
  if (length(rows) == 0) {
    out <- as.data.frame(setNames(replicate(length(columns), logical(0), simplify = FALSE), columns))
    return(out)
  }
  bind_rows(rows)
}

count_source <- function(source_df, pattern, fixed = FALSE) {
  if (is.null(source_df) || !"source" %in% colnames(source_df)) return(NA_integer_)
  if (fixed) {
    sum(source_df$source == pattern, na.rm = TRUE)
  } else {
    sum(grepl(pattern, source_df$source), na.rm = TRUE)
  }
}

standardize_age_corr <- function(df, celltype) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  colnames(df)[colnames(df) %in% c("rho.rho", "rho.S")] <- "rho"
  needed <- c("TF", "rho", "padj")
  if (!all(needed %in% colnames(df))) return(NULL)
  df %>%
    mutate(
      cell_type = celltype,
      rho = as.numeric(rho),
      padj = as.numeric(padj),
      neg_log10_fdr = -log10(padj + 1e-300),
      direction = case_when(
        !is.na(padj) & padj < 0.05 & rho > 0 ~ "Up with age",
        !is.na(padj) & padj < 0.05 & rho < 0 ~ "Down with age",
        TRUE ~ "Stable"
      ),
      abs_rho = abs(rho)
    )
}

get_age_mean_long <- function(celltype, result_dir) {
  mat <- safe_read_rds(file.path(result_dir, "intermediate", "age_mean_auc.rds"))
  if (is.null(mat)) return(NULL)
  mat <- as.matrix(mat)
  data.frame(
    TF = rep(rownames(mat), times = ncol(mat)),
    age = rep(colnames(mat), each = nrow(mat)),
    AUC = as.numeric(mat),
    cell_type = celltype,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      age_numeric = suppressWarnings(as.numeric(as.character(age))),
      cell_type = factor(cell_type, levels = celltype_order)
    )
}

cap_values <- function(x, lower = -2, upper = 2) {
  pmin(pmax(x, lower), upper)
}

format_go_label <- function(x, width = 42) {
  x <- str_replace(
    x,
    "positive regulation of adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains",
    "positive regulation of adaptive immune response"
  )
  str_wrap(x, width = width)
}

make_go_bubble <- function(df, title, label_width = 38) {
  term_order <- df %>%
    dplyr::group_by(Description) %>%
    dplyr::summarise(score = max(neg_log10_p, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(score) %>%
    dplyr::pull(Description)
  
  df <- df %>%
    dplyr::mutate(
      cell_type = factor(cell_type, levels = celltype_order),
      Description = factor(Description, levels = unique(term_order)),
      Description_label = format_go_label(as.character(Description), width = label_width),
      Description_label = factor(
        Description_label,
        levels = unique(format_go_label(term_order, width = label_width))
      )
    )
  
  ggplot(df, aes(x = cell_type, y = Description_label)) +
    geom_point(aes(size = Count, color = neg_log10_p), alpha = 0.88) +
    scale_color_gradient(low = "#E8EDF2", high = "#9A3D3A", name = "-log10(P)") +
    scale_size_continuous(name = "Gene count", range = c(1.5, 7)) +
    labs(
      x = NULL,
      y = "GO biological process",
      title = title
    ) +
    theme_pub(10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 7.5, lineheight = 0.9)
    )
}

save_plot <- function(filename, plot, width, height) {
  ggsave(
    file.path(OUT_DIR, filename),
    plot = plot,
    width = width,
    height = height,
    device = cairo_pdf
  )
}

# ------------------------------------------------------------------------------
# 2. Load all SCENIC result tables
# ------------------------------------------------------------------------------

message("Loading SCENIC results...")

age_corr_all <- list()
cluster_spec_all <- list()
go_all <- list()
age_auc_all <- list()
qc_rows <- list()
regulon_size_rows <- list()

for (ct in names(celltype_dirs)) {
  result_dir <- file.path(BASE_DIR, celltype_dirs[[ct]])
  if (!dir.exists(result_dir)) {
    message("  Skip missing directory: ", result_dir)
    next
  }
  
  age_corr <- standardize_age_corr(
    safe_read_csv(file.path(result_dir, "age_correlation_regulons.csv")),
    ct
  )
  if (!is.null(age_corr)) age_corr_all[[ct]] <- age_corr
  
  cluster_spec <- safe_read_csv(file.path(result_dir, "cluster_specificity_regulons.csv"))
  if (!is.null(cluster_spec) && all(c("TF", "cluster", "log2FC", "padj") %in% colnames(cluster_spec))) {
    cluster_spec_all[[ct]] <- cluster_spec %>%
      mutate(
        cell_type = ct,
        log2FC = as.numeric(log2FC),
        padj = as.numeric(padj),
        neg_log10_fdr = -log10(padj + 1e-300)
      )
  }
  
  go_summary <- safe_read_csv(file.path(result_dir, "GO_enrichment_summary_top5.csv"))
  if (!is.null(go_summary) && all(c("TF", "Description", "pvalue", "Count") %in% colnames(go_summary))) {
    go_all[[ct]] <- go_summary %>%
      mutate(
        cell_type = ct,
        pvalue = as.numeric(pvalue),
        Count = as.numeric(Count),
        neg_log10_p = -log10(pvalue + 1e-300)
      )
  }
  
  age_auc <- get_age_mean_long(ct, result_dir)
  if (!is.null(age_auc)) age_auc_all[[ct]] <- age_auc
  
  regulons <- safe_read_rds(file.path(result_dir, "intermediate", "regulons_clean.rds"))
  if (is.null(regulons)) {
    regulons <- safe_read_rds(file.path(result_dir, "intermediate", "regulons.rds"))
  }
  
  source_df <- safe_read_csv(file.path(result_dir, "regulon_source_detailed.csv"))
  if (!is.null(regulons)) {
    regulon_size_df <- data.frame(
      cell_type = ct,
      TF = names(regulons),
      n_targets = as.integer(lengths(regulons)),
      stringsAsFactors = FALSE
    )
    
    if (!is.null(source_df) && all(c("TF", "n_raw", "source") %in% colnames(source_df))) {
      source_meta <- source_df %>%
        dplyr::mutate(
          TF = as.character(TF),
          n_raw = as.integer(n_raw),
          source = as.character(source)
        ) %>%
        dplyr::select(TF, n_raw, source)
      
      regulon_size_df <- regulon_size_df %>%
        dplyr::left_join(source_meta, by = "TF")
    } else {
      regulon_size_df$n_raw <- NA_integer_
      regulon_size_df$source <- NA_character_
    }
    
    regulon_size_df <- regulon_size_df %>%
      dplyr::mutate(
        regulon_class = dplyr::case_when(
          !is.na(n_raw) & n_raw <= 50 ~ "top50-derived",
          !is.na(n_raw) & n_raw >= 500 ~ "top500-derived",
          is.na(n_raw) & n_targets < 100 ~ "compact",
          is.na(n_raw) & n_targets >= 100 ~ "broad",
          TRUE ~ "other"
        )
      )
    
    regulon_size_rows[[ct]] <- regulon_size_df
  }
  
  qc_rows[[ct]] <- data.frame(
    cell_type = ct,
    n_regulons = if (is.null(regulons)) NA_integer_ else length(regulons),
    median_targets = if (is.null(regulons)) NA_real_ else median(lengths(regulons)),
    motif_pruned = count_source(source_df, "motif_pruned", fixed = TRUE),
    fallback = count_source(source_df, "fallback"),
    stringsAsFactors = FALSE
  )
}

age_corr_all <- bind_or_empty(
  age_corr_all,
  c("TF", "rho", "padj", "cell_type", "neg_log10_fdr", "direction", "abs_rho")
) %>%
  mutate(cell_type = factor(cell_type, levels = celltype_order))

cluster_spec_all <- bind_or_empty(
  cluster_spec_all,
  c("TF", "cluster", "log2FC", "padj", "cell_type", "neg_log10_fdr")
) %>%
  mutate(cell_type = factor(cell_type, levels = celltype_order))

go_all <- bind_or_empty(
  go_all,
  c("TF", "Description", "pvalue", "Count", "cell_type", "neg_log10_p")
) %>%
  mutate(cell_type = factor(cell_type, levels = celltype_order))

age_auc_all <- bind_or_empty(
  age_auc_all,
  c("TF", "age", "AUC", "cell_type", "age_numeric")
)
regulon_size_all <- bind_or_empty(
  regulon_size_rows,
  c("cell_type", "TF", "n_targets", "n_raw", "source", "regulon_class")
) %>%
  mutate(
    cell_type = factor(cell_type, levels = celltype_order),
    regulon_class = factor(
      regulon_class,
      levels = c("top50-derived", "compact", "top500-derived", "broad", "other")
    )
  )

qc_summary <- bind_or_empty(
  qc_rows,
  c("cell_type", "n_regulons", "median_targets", "motif_pruned", "fallback")
) %>%
  mutate(cell_type = factor(cell_type, levels = celltype_order)) %>%
  arrange(cell_type)

write_csv(qc_summary, file.path(TABLE_DIR, "part2_scenic_qc_summary.csv"))
write_csv(age_corr_all, file.path(TABLE_DIR, "part2_age_correlation_all_celltypes.csv"))
write_csv(cluster_spec_all, file.path(TABLE_DIR, "part2_cluster_specificity_all_celltypes.csv"))
write_csv(go_all, file.path(TABLE_DIR, "part2_go_summary_all_celltypes.csv"))
write_csv(regulon_size_all, file.path(TABLE_DIR, "part2_regulon_size_distribution.csv"))

if (nrow(age_corr_all) == 0) {
  stop("No age_correlation_regulons.csv files were found. Please check BASE_DIR and result folders.")
}

# ------------------------------------------------------------------------------
# 3. Fig2A: SCENIC QC summary
# ------------------------------------------------------------------------------

message("Plotting Fig2A...")

p_qc_regulons <- ggplot(qc_summary, aes(x = cell_type, y = n_regulons, fill = cell_type)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.25) +
  scale_fill_manual(values = celltype_palette, guide = "none") +
  labs(x = NULL, y = "No. regulons", title = "Regulon recovery") +
  theme_pub(10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_qc_targets <- ggplot(regulon_size_all, aes(x = cell_type, y = n_targets)) +
  geom_boxplot(
    aes(fill = cell_type),
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.45,
    color = "grey35",
    linewidth = 0.35
  ) +
  geom_jitter(
    aes(color = regulon_class),
    width = 0.18,
    height = 0,
    size = 1.35,
    alpha = 0.82
  ) +
  scale_fill_manual(values = celltype_palette, guide = "none") +
  scale_color_manual(
    values = c(
      "top50-derived" = "#4D83B8",
      "compact" = "#4D83B8",
      "top500-derived" = "#C05A4B",
      "broad" = "#C05A4B",
      "other" = "grey55"
    ),
    breaks = c("top50-derived", "compact", "top500-derived", "broad"),
    labels = c("top50-derived", "compact (<100 targets)", "top500-derived", "broad (>=100 targets)"),
    name = "Regulon class"
  ) +
  scale_y_log10(
    breaks = c(10, 30, 50, 100, 300, 500),
    labels = c("10", "30", "50", "100", "300", "500")
  ) +
  labs(x = NULL, y = "Targets per regulon (log10)", title = "Regulon size distribution") +
  theme_pub(10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

fig2a <- p_qc_regulons | p_qc_targets
save_plot("Fig2A_SCENIC_QC_summary.pdf", fig2a, width = 11, height = 4.6)

# ------------------------------------------------------------------------------
# 4. Fig2B: Age-associated regulon landscape
# ------------------------------------------------------------------------------

message("Plotting Fig2B...")

age_count <- age_corr_all %>%
  filter(direction != "Stable") %>%
  count(cell_type, direction, name = "n") %>%
  complete(cell_type, direction = c("Up with age", "Down with age"), fill = list(n = 0)) %>%
  mutate(
    signed_n = ifelse(direction == "Down with age", -n, n),
    cell_type = factor(cell_type, levels = celltype_order)
  )

celltype_by_age_signal <- age_count %>%
  group_by(cell_type) %>%
  summarise(total_dynamic = sum(n), .groups = "drop") %>%
  arrange(total_dynamic) %>%
  pull(cell_type)

age_count <- age_count %>%
  mutate(cell_type_plot = factor(cell_type, levels = celltype_by_age_signal))

p_age_count <- ggplot(age_count, aes(x = signed_n, y = cell_type_plot, fill = direction)) +
  geom_vline(xintercept = 0, color = "grey30", linewidth = 0.35) +
  geom_col(width = 0.68, color = "white", linewidth = 0.25) +
  scale_fill_manual(values = direction_palette[c("Up with age", "Down with age")]) +
  scale_x_continuous(labels = abs, expand = expansion(mult = c(0.08, 0.08))) +
  labs(
    x = "No. significant regulons",
    y = NULL,
    title = "Age-associated regulon dynamics",
    fill = NULL
  ) +
  theme_pub(10) +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.25)
  )

top_age <- age_corr_all %>%
  filter(!is.na(padj), !is.na(rho)) %>%
  group_by(cell_type) %>%
  slice_max(order_by = abs_rho, n = 4, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    cell_type = factor(cell_type, levels = celltype_order),
    label = paste(cell_type, TF, sep = " | "),
    label = factor(label, levels = rev(unique(label[order(cell_type, abs_rho)]))),
    direction = ifelse(rho >= 0, "Up with age", "Down with age")
  )

p_age_dot <- ggplot(top_age, aes(x = rho, y = label)) +
  geom_vline(xintercept = 0, color = "grey80", linewidth = 0.35) +
  geom_point(aes(size = neg_log10_fdr, color = rho), alpha = 0.9) +
  scale_color_gradient2(
    low = direction_palette[["Down with age"]], mid = "grey96", high = direction_palette[["Up with age"]],
    midpoint = 0, name = "Spearman rho"
  ) +
  scale_size_continuous(name = "-log10(FDR)", range = c(1.4, 5.4)) +
  scale_x_continuous(limits = c(-0.8, 0.8), breaks = c(-0.6, -0.3, 0, 0.3, 0.6)) +
  labs(
    x = "Spearman rho with age",
    y = "Top age-associated regulons",
    title = "Leading age-associated TF regulons"
  ) +
  theme_pub(10) +
  theme(
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.25),
    axis.text.y = element_text(size = 8)
  )

fig2b <- p_age_count | p_age_dot
save_plot("Fig2B_age_regulon_landscape.pdf", fig2b, width = 12, height = 7.8)
write_csv(age_count, file.path(TABLE_DIR, "part2_age_regulon_counts.csv"))

# ------------------------------------------------------------------------------
# 5. Fig2C: Heatmap of top age-dynamic regulons across ages
# ------------------------------------------------------------------------------

message("Plotting Fig2C...")

top_heatmap_tfs <- age_corr_all %>%
  filter(!is.na(padj), !is.na(rho), padj < 0.05) %>%
  group_by(cell_type) %>%
  slice_max(order_by = abs_rho, n = 8, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(cell_type, TF, rho, padj, direction)

heatmap_long <- age_auc_all %>%
  inner_join(top_heatmap_tfs, by = c("cell_type", "TF")) %>%
  filter(!is.na(age_numeric))

if (nrow(heatmap_long) > 0) {
  heatmap_scaled <- heatmap_long %>%
    group_by(cell_type, TF) %>%
    mutate(AUC_z = as.numeric(scale(AUC))) %>%
    ungroup() %>%
    mutate(
      row_id = paste(cell_type, TF, sep = " | "),
      AUC_z = cap_values(AUC_z, -2, 2)
    ) %>%
    dplyr::select(row_id, age_numeric, AUC_z) %>%
    distinct() %>%
    pivot_wider(names_from = age_numeric, values_from = AUC_z)
  
  heatmap_mat <- as.data.frame(heatmap_scaled)
  rownames(heatmap_mat) <- heatmap_mat$row_id
  heatmap_mat$row_id <- NULL
  heatmap_mat <- as.matrix(heatmap_mat)
  heatmap_mat <- heatmap_mat[, order(as.numeric(colnames(heatmap_mat))), drop = FALSE]
  
  row_anno <- top_heatmap_tfs %>%
    mutate(row_id = paste(cell_type, TF, sep = " | ")) %>%
    distinct(row_id, cell_type, direction) %>%
    as.data.frame()
  rownames(row_anno) <- row_anno$row_id
  row_anno$row_id <- NULL
  row_anno <- row_anno[rownames(heatmap_mat), , drop = FALSE]
  
  anno_colors <- list(
    cell_type = celltype_palette,
    direction = direction_annotation_palette
  )
  
  pdf(file.path(OUT_DIR, "Fig2C_age_dynamic_regulon_heatmap.pdf"),
      width = 8.5, height = max(7, nrow(heatmap_mat) * 0.18))
  pheatmap(
    heatmap_mat,
    color = colorRampPalette(c(direction_palette[["Down with age"]], "white", direction_palette[["Up with age"]]))(100),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_row = row_anno,
    annotation_colors = anno_colors,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 7,
    fontsize_col = 10,
    border_color = NA,
    main = "Age-dynamic regulon activity across normal development"
  )
  dev.off()
  
  write_csv(top_heatmap_tfs, file.path(TABLE_DIR, "part2_top_age_dynamic_regulons_for_heatmap.csv"))
} else {
  message("  No significant age-dynamic regulons available for Fig2C.")
}

# ------------------------------------------------------------------------------
# 6. Fig2D: GO bubble plot for age-dynamic regulons
# ------------------------------------------------------------------------------

message("Plotting Fig2D...")

if (nrow(go_all) > 0) {
  age_dynamic_tf <- age_corr_all %>%
    filter(!is.na(padj), padj < 0.05) %>%
    dplyr::select(cell_type, TF, rho, padj, direction)
  
  go_age_dynamic <- go_all %>%
    inner_join(age_dynamic_tf, by = c("cell_type", "TF")) %>%
    mutate(
      main_include = str_detect(
        str_to_lower(Description),
        paste(str_to_lower(go_main_include_patterns), collapse = "|")
      ),
      main_exclude = str_detect(
        str_to_lower(Description),
        paste(str_to_lower(go_main_exclude_patterns), collapse = "|")
      )
    ) %>%
    group_by(cell_type, Description) %>%
    summarise(
      Count = max(Count, na.rm = TRUE),
      neg_log10_p = max(neg_log10_p, na.rm = TRUE),
      n_TF = n_distinct(TF),
      main_include = any(main_include),
      main_exclude = any(main_exclude),
      .groups = "drop"
    )
  
  go_main <- go_age_dynamic %>%
    filter(main_include, !main_exclude) %>%
    arrange(desc(n_TF), desc(neg_log10_p)) %>%
    group_by(cell_type) %>%
    slice_head(n = 4) %>%
    ungroup()
  
  if (nrow(go_main) == 0) {
    message("  No GO terms passed the main-figure biological filter; falling back to top enriched terms.")
    go_main <- go_age_dynamic %>%
      filter(!main_exclude) %>%
      arrange(desc(n_TF), desc(neg_log10_p)) %>%
      group_by(cell_type) %>%
      slice_head(n = 4) %>%
      ungroup()
  }
  
  # Keep the main figure compact by retaining the strongest terms globally.
  main_terms <- go_main %>%
    group_by(Description) %>%
    summarise(score = max(neg_log10_p, na.rm = TRUE), n_celltypes = n_distinct(cell_type), .groups = "drop") %>%
    arrange(desc(n_celltypes), desc(score)) %>%
    slice_head(n = 24) %>%
    pull(Description)
  
  go_main <- go_main %>%
    filter(Description %in% main_terms)
  
  go_full <- go_age_dynamic %>%
    arrange(desc(main_include), main_exclude, desc(n_TF), desc(neg_log10_p)) %>%
    group_by(cell_type) %>%
    slice_head(n = 8) %>%
    ungroup()
  
  p_go_main <- make_go_bubble(
    go_main,
    "Functional programs linked to age-dynamic regulons",
    label_width = 36
  )
  save_plot("Fig2D_GO_bubble_age_dynamic_regulons.pdf", p_go_main, width = 9.5, height = 7.5)
  
  p_go_full <- make_go_bubble(
    go_full,
    "GO enrichment of age-dynamic regulons",
    label_width = 42
  )
  save_plot("FigS2_GO_bubble_age_dynamic_regulons_full.pdf", p_go_full, width = 10.5, height = 9.5)
  
  write_csv(go_main, file.path(TABLE_DIR, "part2_go_bubble_terms_main_figure.csv"))
  write_csv(go_full, file.path(TABLE_DIR, "part2_go_bubble_terms_full_supplement.csv"))
  write_csv(go_age_dynamic, file.path(TABLE_DIR, "part2_go_age_dynamic_all_terms.csv"))
  write_csv(go_main, file.path(TABLE_DIR, "part2_go_bubble_terms.csv"))
} else {
  message("  GO summary table not found. Run the GO enrichment step first.")
}

# ------------------------------------------------------------------------------
# 7. Fig2E: Key TF regulon trajectories across age
# ------------------------------------------------------------------------------

message("Plotting Fig2E...")

trajectory_tfs <- bind_rows(lapply(celltype_order, function(ct) {
  available_tfs <- age_auc_all %>%
    filter(cell_type == ct, !is.na(age_numeric)) %>%
    pull(TF) %>%
    unique()
  
  curated_tfs <- intersect(key_tf_map[[ct]], available_tfs)
  
  age_tfs <- age_corr_all %>%
    filter(cell_type == ct, TF %in% available_tfs, !is.na(rho), !is.na(padj)) %>%
    arrange(desc(abs_rho)) %>%
    pull(TF) %>%
    unique()
  
  keep_tfs <- unique(c(curated_tfs, age_tfs))
  keep_tfs <- head(keep_tfs, 4)
  
  if (length(keep_tfs) == 0) return(NULL)
  
  data.frame(cell_type = ct, TF = keep_tfs, stringsAsFactors = FALSE)
}))

trajectory_df <- age_auc_all %>%
  inner_join(trajectory_tfs, by = c("cell_type", "TF")) %>%
  filter(!is.na(age_numeric))

if (nrow(trajectory_df) > 0) {
  trajectory_df <- trajectory_df %>%
    group_by(cell_type, TF) %>%
    filter(n_distinct(age_numeric) >= 3) %>%
    mutate(AUC_z = as.numeric(scale(AUC))) %>%
    ungroup() %>%
    mutate(
      cell_type = factor(cell_type, levels = celltype_order),
      trend = ifelse(
        age_numeric == max(age_numeric, na.rm = TRUE),
        AUC_z - AUC_z[which.min(age_numeric)],
        NA_real_
      )
    )
  
  label_df <- trajectory_df %>%
    group_by(cell_type, TF) %>%
    filter(age_numeric == max(age_numeric, na.rm = TRUE)) %>%
    slice_tail(n = 1) %>%
    ungroup()
  
  p_traj <- ggplot(trajectory_df, aes(x = age_numeric, y = AUC_z, group = TF)) +
    geom_hline(yintercept = 0, color = "grey88", linewidth = 0.3) +
    geom_line(aes(color = TF), linewidth = 0.9, alpha = 0.92) +
    geom_point(aes(fill = TF), shape = 21, color = "white", stroke = 0.25, size = 2.1) +
    geom_text_repel(
      data = label_df,
      aes(label = TF, color = TF),
      size = 3,
      fontface = "italic",
      direction = "y",
      nudge_x = 0.75,
      segment.color = "grey75",
      segment.linewidth = 0.25,
      min.segment.length = 0,
      box.padding = 0.25,
      show.legend = FALSE
    ) +
    facet_wrap(~cell_type, scales = "free_y", ncol = 3) +
    scale_x_continuous(breaks = sort(unique(trajectory_df$age_numeric)), expand = expansion(mult = c(0.03, 0.18))) +
    scale_color_manual(values = rep(c("#BB6A8A", "#91B0CB", "#6FA39A", "#C7A15A", "#8F7AAE", "#D98B72"), 20)) +
    scale_fill_manual(values = rep(c("#BB6A8A", "#91B0CB", "#6FA39A", "#C7A15A", "#8F7AAE", "#D98B72"), 20)) +
    labs(
      x = "Age",
      y = "Regulon activity z-score",
      title = "Developmental trajectories of selected TF regulons"
    ) +
    theme_pub(9) +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_line(color = "grey92", linewidth = 0.25),
      strip.text = element_text(size = 9, face = "bold")
    )
  
  save_plot("Fig2E_key_TF_age_trajectories.pdf", p_traj, width = 11.5, height = 8.2)
  write_csv(trajectory_tfs, file.path(TABLE_DIR, "part2_key_tf_trajectories_tfs.csv"))
} else {
  message("  No key TF trajectories available for Fig2E.")
}

# ------------------------------------------------------------------------------
# 8. Supplementary: cluster marker regulon bubble plot
# ------------------------------------------------------------------------------

message("Plotting supplementary cluster regulon bubble...")

if (nrow(cluster_spec_all) > 0) {
  cluster_plot <- cluster_spec_all %>%
    filter(!is.na(padj), padj < 0.05, log2FC > 0.5) %>%
    group_by(cell_type, cluster) %>%
    slice_max(order_by = log2FC, n = 4, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      cluster_id = paste(cell_type, cluster, sep = " C"),
      TF = fct_reorder(TF, log2FC, .fun = max)
    )
  
  if (nrow(cluster_plot) > 0) {
    p_cluster <- ggplot(cluster_plot, aes(x = cluster_id, y = TF)) +
      geom_point(aes(size = neg_log10_fdr, color = log2FC), alpha = 0.9) +
      scale_color_gradient(low = "#F3E7C6", high = "#B9473D", name = "log2FC") +
      scale_size_continuous(name = "-log10(FDR)", range = c(1.2, 6)) +
      labs(
        x = "Cell type-specific cluster",
        y = "Marker regulon",
        title = "Cluster-specific marker regulons"
      ) +
      theme_pub(9) +
      theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7))
    
    save_plot("FigS2_cluster_marker_regulon_bubble.pdf", p_cluster, width = 12, height = 9)
    write_csv(cluster_plot, file.path(TABLE_DIR, "part2_cluster_marker_regulons_for_bubble.csv"))
  }
}

message("Done. Figures saved to: ", normalizePath(OUT_DIR, winslash = "/"))

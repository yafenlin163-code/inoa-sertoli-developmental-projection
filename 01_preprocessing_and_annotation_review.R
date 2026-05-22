

# -------------------------------------------------------------------------
# 01_preprocessing_and_annotation_review.R
#
# Purpose:
#   Review quality-control metrics and major testicular cell-type annotations
#   in public processed human testicular single-cell objects, then generate
#   marker-validation and overview plots used for the manuscript.
#
# Scope:
#   This script starts from processed/annotated Seurat objects. It does not
#   perform FASTQ alignment, raw count generation, or full de novo preprocessing.
#
# Expected inputs:
#   - merged_normal Seurat object containing normal postnatal testicular cells
#   - optional merged_raw and merged_qc Seurat objects for QC comparison plots
#
# Default input paths can be overridden with environment variables:
#   MERGED_NORMAL_RDS, MERGED_RAW_RDS, MERGED_QC_RDS, PROJECT_DIR, OUTPUT_DIR
# -------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(ggplotify)
  library(scales)
  library(viridis)
})

set.seed(42)

project_dir <- Sys.getenv("PROJECT_DIR", unset = getwd())
output_dir <- Sys.getenv(
  "OUTPUT_DIR",
  unset = file.path(project_dir, "results", "01_preprocessing_and_annotation_review")
)
plot_dir <- file.path(output_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

merged_normal_rds <- Sys.getenv(
  "MERGED_NORMAL_RDS",
  unset = file.path(project_dir, "data", "processed", "merged_normal.rds")
)
merged_raw_rds <- Sys.getenv(
  "MERGED_RAW_RDS",
  unset = file.path(project_dir, "data", "processed", "merged_raw.rds")
)
merged_qc_rds <- Sys.getenv(
  "MERGED_QC_RDS",
  unset = file.path(project_dir, "data", "processed", "merged_qc.rds")
)

load_object <- function(object_name, path, required = TRUE) {
  if (exists(object_name, envir = .GlobalEnv)) {
    return(get(object_name, envir = .GlobalEnv))
  }
  if (file.exists(path)) {
    return(readRDS(path))
  }
  if (required) {
    stop(
      object_name,
      " was not found in the current R session and the file does not exist: ",
      path,
      call. = FALSE
    )
  }
  NULL
}

merged_normal <- load_object("merged_normal", merged_normal_rds, required = TRUE)
merged_raw <- load_object("merged_raw", merged_raw_rds, required = FALSE)
merged_qc <- load_object("merged_qc", merged_qc_rds, required = FALSE)
has_qc_objects <- !is.null(merged_raw) && !is.null(merged_qc)

required_meta <- c("sample", "cell_type", "age_numeric")
missing_meta <- setdiff(required_meta, colnames(merged_normal@meta.data))
if (length(missing_meta) > 0) {
  stop(
    "merged_normal is missing required metadata columns: ",
    paste(missing_meta, collapse = ", "),
    call. = FALSE
  )
}

if (!"seurat_clusters" %in% colnames(merged_normal@meta.data)) {
  merged_normal$seurat_clusters <- as.character(Idents(merged_normal))
}

merged_normal$cell_type <- as.character(merged_normal$cell_type)
present_types <- sort(unique(merged_normal$cell_type))

age_levels <- unique(merged_normal$age_numeric)
if (all(!is.na(suppressWarnings(as.numeric(as.character(age_levels)))))) {
  age_levels <- as.character(sort(unique(as.numeric(as.character(age_levels)))))
} else {
  age_levels <- as.character(age_levels)
}
merged_normal$age_group_for_plot <- factor(as.character(merged_normal$age_numeric), levels = age_levels)

base_celltype_colors <- c(
  Sertoli = "#1B9E77",
  Leydig = "#D95F02",
  PMC = "#7570B3",
  Endothelial = "#E7298A",
  Macrophage = "#66A61E",
  T_cell = "#E6AB02",
  Spermatogonia = "#A6761D",
  Spermatocyte = "#666666",
  Spermatid = "#1F78B4"
)

missing_color_types <- setdiff(present_types, names(base_celltype_colors))
if (length(missing_color_types) > 0) {
  extra_colors <- hue_pal()(length(missing_color_types))
  names(extra_colors) <- missing_color_types
  base_celltype_colors <- c(base_celltype_colors, extra_colors)
}
celltype_colors_use <- base_celltype_colors[present_types]

# -------------------------------------------------------------------------
# Marker genes
# -------------------------------------------------------------------------
testis_markers <- list(
  Sertoli = c("SOX9", "AMH", "CLDN11", "WT1", "FSHR", "GATA4", "VIM"),
  Leydig = c("CYP11A1", "CYP17A1", "STAR", "HSD3B1", "INSL3", "NR5A1"),
  PMC = c("ACTA2", "MYH11", "MYOCD", "DES", "CNN1"),
  Endothelial = c("PECAM1", "CD34", "VWF", "ENG", "CDH5"),
  Macrophage = c("CD68", "CD163", "MRC1", "CSF1R", "ADGRE1"),
  T_cell = c("CD3D", "CD3E", "CD3G", "CD8A", "CD4"),
  Spermatogonia = c("MAGEA4", "DAZL", "STRA8", "ID4", "FGFR3", "GFRA1", "UTF1"),
  Spermatocyte = c("SYCP1", "SYCP3", "MLH3", "SPO11", "MEIOC", "HORMAD1"),
  Spermatid = c("TNP1", "TNP2", "PRM1", "PRM2", "ACR", "ACRV1", "SPEM1")
)

top_markers <- lapply(testis_markers, function(genes) {
  intersect(genes, rownames(merged_normal))
})
dot_genes <- unique(unlist(lapply(top_markers, head, 3), use.names = FALSE))

if (length(dot_genes) == 0) {
  stop("No marker genes were found in merged_normal.", call. = FALSE)
}

# -------------------------------------------------------------------------
# Compute t-SNE before figure generation
# -------------------------------------------------------------------------
if (!"tsne" %in% names(merged_normal@reductions)) {
  merged_normal <- RunTSNE(
    merged_normal,
    dims = 1:30,
    perplexity = 30,
    seed.use = 42,
    check_duplicates = FALSE
  )
}

# -------------------------------------------------------------------------
# Main Figure 1: marker validation
# -------------------------------------------------------------------------
p_dot <- DotPlot(
  merged_normal,
  features = dot_genes,
  group.by = "seurat_clusters"
) +
  RotatedAxis() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  ggtitle("Marker gene expression by cluster") +
  theme(axis.text.x = element_text(size = 7))

avg_expr <- tryCatch(
  {
    AverageExpression(
      merged_normal,
      features = dot_genes,
      group.by = "seurat_clusters",
      assays = "RNA",
      layer = "data"
    )$RNA
  },
  error = function(e) {
    AverageExpression(
      merged_normal,
      features = dot_genes,
      group.by = "seurat_clusters",
      assays = "RNA",
      slot = "data"
    )$RNA
  }
)

avg_scaled <- t(scale(t(as.matrix(avg_expr))))
avg_scaled[is.na(avg_scaled)] <- 0
avg_scaled <- pmin(pmax(avg_scaled, -2), 2)

ph <- pheatmap(
  avg_scaled,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  fontsize = 8,
  main = "Average marker expression across clusters (z-score)",
  silent = TRUE
)
p_heatmap <- as.ggplot(ph)

p_marker <- p_dot / p_heatmap +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(
    title = "Marker-based validation of cell type annotation",
    subtitle = "Top: DotPlot (dot size = detection rate; color = average expression)\nBottom: Heatmap (cluster-level z-score)"
  )

ggsave(
  file.path(plot_dir, "Fig1_marker_validation.pdf"),
  p_marker,
  width = 18,
  height = 18
)

# -------------------------------------------------------------------------
# Main Figure 2: annotated t-SNE
# -------------------------------------------------------------------------
p_ct <- DimPlot(
  merged_normal,
  reduction = "tsne",
  group.by = "cell_type",
  cols = celltype_colors_use,
  pt.size = 0.4,
  label = TRUE,
  label.size = 4,
  repel = TRUE
) +
  ggtitle("Annotated t-SNE of major cell types") +
  theme(legend.position = "right")

ggsave(
  file.path(plot_dir, "Fig2_annotated_tsne.pdf"),
  p_ct,
  width = 10,
  height = 8
)

# -------------------------------------------------------------------------
# Main Figure 3: age-stratified t-SNE
# -------------------------------------------------------------------------
p_ct_age <- DimPlot(
  merged_normal,
  reduction = "tsne",
  group.by = "cell_type",
  cols = celltype_colors_use,
  split.by = "age_group_for_plot",
  pt.size = 0.3,
  ncol = 3
) +
  ggtitle("Cell type distribution across age groups (t-SNE)")

ggsave(
  file.path(plot_dir, "Fig3_age_stratified_tsne.pdf"),
  p_ct_age,
  width = 18,
  height = 10
)

# -------------------------------------------------------------------------
# Main Figure 4: cell-type composition
# -------------------------------------------------------------------------
prop_data <- merged_normal@meta.data %>%
  group_by(sample, cell_type, age_numeric) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(sample = factor(sample, levels = unique(sample[order(age_numeric)])))

p_prop <- ggplot(prop_data, aes(x = sample, y = prop, fill = cell_type)) +
  geom_col(position = "stack", width = 0.8) +
  scale_fill_manual(values = celltype_colors_use) +
  scale_y_continuous(labels = percent) +
  labs(
    x = "Sample (ordered by age)",
    y = "Cell proportion",
    fill = "Cell type",
    title = "Cell type composition across samples"
  ) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(plot_dir, "Fig4_celltype_composition.pdf"),
  p_prop,
  width = 12,
  height = 6
)

# -------------------------------------------------------------------------
# Main Figure 5: selected target cell types on global t-SNE
# -------------------------------------------------------------------------
target_celltypes <- c("Sertoli", "Leydig", "PMC", "Endothelial")
target_celltypes <- target_celltypes[target_celltypes %in% present_types]

highlight_plots <- lapply(target_celltypes, function(cell_type_name) {
  merged_normal$highlight_tmp <- factor(
    ifelse(merged_normal$cell_type == cell_type_name, cell_type_name, "Other"),
    levels = c("Other", cell_type_name)
  )

  color_map <- c("Other" = "grey85", "#4DAF4A")
  names(color_map)[2] <- cell_type_name

  DimPlot(
    merged_normal,
    reduction = "tsne",
    group.by = "highlight_tmp",
    pt.size = 0.3,
    order = cell_type_name
  ) +
    scale_color_manual(values = color_map) +
    ggtitle(paste0(cell_type_name, " in global t-SNE")) +
    theme_classic(base_size = 12) +
    theme(legend.position = "none")
})

if (length(highlight_plots) > 0) {
  p_highlight <- wrap_plots(highlight_plots, ncol = 2)

  ggsave(
    file.path(plot_dir, "Fig5_selected_celltypes_global_tsne.pdf"),
    p_highlight,
    width = 12,
    height = 10
  )
}

# -------------------------------------------------------------------------
# Supplementary Figure S1: QC metrics
# -------------------------------------------------------------------------
if (has_qc_objects) {
  sample_order <- c(
    "2_years", "5_years", "8_years", "11_years",
    "Young1", "Young2", "Young3", "Young4",
    "KS_1", "KS_2", "KS_3", "iNOA_1", "iNOA_2", "iNOA_3", "AZFa_Del_1"
  )

  sample_order_present_raw <- intersect(sample_order, unique(merged_raw$sample))
  merged_raw$sample <- factor(merged_raw$sample, levels = sample_order_present_raw)

  sample_order_present_qc <- intersect(sample_order, unique(merged_qc$sample))
  merged_qc$sample <- factor(merged_qc$sample, levels = sample_order_present_qc)

  qc_features <- c("nFeature_RNA", "percent.mt", "percent.ribo", "ISG_score")
  qc_features_raw <- intersect(qc_features, colnames(merged_raw@meta.data))
  qc_features_after <- intersect(qc_features, colnames(merged_qc@meta.data))

  if (length(qc_features_raw) > 0 && length(qc_features_after) > 0) {
    p_qc_raw <- VlnPlot(
      merged_raw,
      features = qc_features_raw,
      group.by = "sample",
      pt.size = 0,
      ncol = 2
    ) +
      plot_annotation(title = "QC metrics before filtering")

    p_qc_after <- VlnPlot(
      merged_qc,
      features = qc_features_after,
      group.by = "sample",
      pt.size = 0,
      ncol = 2
    ) +
      plot_annotation(title = "QC metrics after filtering")

    ggsave(
      file.path(plot_dir, "FigS1_qc_metrics.pdf"),
      p_qc_raw / p_qc_after,
      width = 16,
      height = 18
    )
  }
}

# -------------------------------------------------------------------------
# Supplementary Figure S2: t-SNE overview
# -------------------------------------------------------------------------
p1 <- DimPlot(merged_normal, reduction = "tsne", group.by = "sample", pt.size = 0.3) +
  ggtitle("Sample")
p2 <- DimPlot(merged_normal, reduction = "tsne", group.by = "age_group_for_plot", pt.size = 0.3) +
  scale_color_viridis_d(option = "C") +
  ggtitle("Age")
p3 <- DimPlot(merged_normal, reduction = "tsne", label = TRUE, pt.size = 0.3) +
  ggtitle("Clusters")

if ("percent.mt" %in% colnames(merged_normal@meta.data)) {
  p4 <- FeaturePlot(merged_normal, reduction = "tsne", features = "percent.mt", pt.size = 0.3) +
    ggtitle("Mitochondrial percentage")
} else {
  p4 <- ggplot() +
    theme_void() +
    ggtitle("Mitochondrial percentage unavailable")
}

ggsave(
  file.path(plot_dir, "FigS2_tsne_overview.pdf"),
  (p1 | p2) / (p3 | p4),
  width = 14,
  height = 12
)

# -------------------------------------------------------------------------
# Supplementary Figure S3: module-score t-SNE
# -------------------------------------------------------------------------
score_cols <- grep("_score$", colnames(merged_normal@meta.data), value = TRUE)
if (length(score_cols) > 0) {
  p_scores <- FeaturePlot(
    merged_normal,
    reduction = "tsne",
    features = score_cols,
    ncol = 4,
    pt.size = 0.2,
    order = TRUE
  ) &
    scale_color_gradient2(low = "grey90", mid = "#C0D6EA", high = "#D3635F", midpoint = 0)

  ggsave(
    file.path(plot_dir, "FigS3_module_score_tsne.pdf"),
    p_scores,
    width = 20,
    height = 15
  )
}

# -------------------------------------------------------------------------
# Supplementary Figure S4: sample distributions of selected target cell types
# -------------------------------------------------------------------------
sample_dist_plots <- list()
for (cell_type_name in target_celltypes) {
  csv_file <- file.path(output_dir, paste0(cell_type_name, "_sample_distribution.csv"))
  if (!file.exists(csv_file)) next

  sample_dist <- read.csv(csv_file)
  sample_dist$sample <- factor(sample_dist$sample, levels = unique(sample_dist$sample))

  p_tmp <- ggplot(sample_dist, aes(x = sample, y = Percent_in_sample, fill = category)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(
      title = paste0(cell_type_name, ": proportion within each sample"),
      x = "Sample",
      y = "Proportion (%)",
      fill = "Group"
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top",
      plot.title = element_text(face = "bold", size = 11)
    )

  sample_dist_plots[[cell_type_name]] <- p_tmp
}

if (length(sample_dist_plots) > 0) {
  ggsave(
    file.path(plot_dir, "FigS4_selected_celltype_sample_distribution.pdf"),
    wrap_plots(sample_dist_plots, ncol = 2),
    width = 14,
    height = 10
  )
}

# -------------------------------------------------------------------------
# Supplementary Figure S5: integrated t-SNE overview
# -------------------------------------------------------------------------
p_tsne1 <- DimPlot(
  merged_normal,
  reduction = "tsne",
  group.by = "cell_type",
  cols = celltype_colors_use,
  pt.size = 0.3,
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Cell type annotation (t-SNE)")
p_tsne2 <- DimPlot(merged_normal, reduction = "tsne", group.by = "sample", pt.size = 0.3) +
  ggtitle("Sample distribution (t-SNE)")
p_tsne3 <- DimPlot(merged_normal, reduction = "tsne", group.by = "age_group_for_plot", pt.size = 0.3) +
  scale_color_viridis_d(option = "C") +
  ggtitle("Age distribution (t-SNE)")
p_tsne4 <- p4

ggsave(
  file.path(plot_dir, "FigS5_tsne_overview.pdf"),
  (p_tsne1 | p_tsne2) / (p_tsne3 | p_tsne4),
  width = 16,
  height = 14
)

message("Annotation review and visualization script finished.")

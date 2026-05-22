#!/usr/bin/env Rscript

# ==============================================================================
# 06_sertoli_evidence_synthesis.R
#
# Purpose:
#   Integrate the normal Sertoli developmental reference, disease projection
#   summaries, and optional Part 5 molecular characterization outputs into a
#   focused evidence chain for an immature-like, stress-biased Sertoli-cell
#   phenotype in idiopathic non-obstructive azoospermia.
#
# Important interpretation:
#   This script strengthens a developmental-reference association model. It does
#   not by itself prove causal developmental arrest.
#
# Expected inputs, discovered automatically when present:
#   results/03_developmental_reference_construction/tables/*.csv
#   results/04_disease_cell_projection/tables/*.csv
#   results/05_sertoli_molecular_characterization/tables/*.csv      optional
#   Legacy Part 3-5 output directory names are also supported.
#
# Outputs:
#   results/06_sertoli_evidence_synthesis/tables/
#   results/06_sertoli_evidence_synthesis/figures/main/
#   results/06_sertoli_evidence_synthesis/figures/supplementary/
#   results/06_sertoli_evidence_synthesis/reports/
#
# Default paths can be overridden with environment variables:
#   PROJECT_DIR, PART6_BASE_DIR, PART6_OUTPUT_DIR
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

HAS_GGREPEL <- requireNamespace("ggrepel", quietly = TRUE)
HAS_PATCHWORK <- requireNamespace("patchwork", quietly = TRUE)

# ------------------------------------------------------------------------------
# 0. Configuration and path discovery
# ------------------------------------------------------------------------------

USER_BASE_DIR <- Sys.getenv("PART6_BASE_DIR", unset = Sys.getenv("PROJECT_DIR", unset = ""))

rel_exists <- function(base, rel_paths) {
  any(file.exists(file.path(base, rel_paths)))
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
    "results/03_developmental_reference_construction/tables/all_celltypes_pseudotime_age_summary.csv",
    "part3/tables/all_celltypes_pseudotime_age_summary.csv",
    "part3_multicelltype_pseudotime_results/tables/all_celltypes_pseudotime_age_summary.csv"
  )
  part4_rel <- c(
    "results/04_disease_cell_projection/tables/disease_projection_summary_by_sample.csv",
    "part4/disease_projection_summary_by_sample.csv",
    "part4/tables/disease_projection_summary_by_sample.csv",
    "part4_disease_projection_results/tables/disease_projection_summary_by_sample.csv"
  )

  seed_dirs <- unique(c(
    USER_BASE_DIR,
    getwd(),
    get_script_dir()
  ))
  seed_dirs <- seed_dirs[!is.na(seed_dirs) & nzchar(seed_dirs)]
  seed_dirs <- seed_dirs[dir.exists(seed_dirs)]
  candidate_dirs <- unique(unlist(lapply(seed_dirs, parent_dirs), use.names = FALSE))

  for (base in candidate_dirs) {
    if (rel_exists(base, part3_rel) && rel_exists(base, part4_rel)) return(base)
  }

  stop(
    "Could not detect the project base directory.\n",
    "Please either set your working directory to the project root before source(), e.g.\n",
    "  setwd('/path/to/Regulatory rewiring of Sertoli cell maturation')\n",
    "or set:\n",
    "  Sys.setenv(PART6_BASE_DIR = '/path/to/Regulatory rewiring of Sertoli cell maturation')\n",
    "The project root should contain Part 3 and Part 4 result files."
  )
}

BASE_DIR <- normalizePath(detect_base_dir(), winslash = "/", mustWork = TRUE)

is_absolute_path <- function(path) {
  grepl("^[A-Za-z]:[/\\\\]", path) || startsWith(path, "/") || startsWith(path, "\\\\")
}

OUTDIR_SETTING <- Sys.getenv("PART6_OUTPUT_DIR", unset = file.path("results", "06_sertoli_evidence_synthesis"))
OUTDIR <- if (is_absolute_path(OUTDIR_SETTING)) OUTDIR_SETTING else file.path(BASE_DIR, OUTDIR_SETTING)
TABLE_DIR <- file.path(OUTDIR, "tables")
FIG_DIR <- file.path(OUTDIR, "figures")
MAIN_FIG_DIR <- file.path(FIG_DIR, "main")
SUPP_FIG_DIR <- file.path(FIG_DIR, "supplementary")
REPORT_DIR <- file.path(OUTDIR, "reports")

dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(MAIN_FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(SUPP_FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(REPORT_DIR, showWarnings = FALSE, recursive = TRUE)

SERTOLI_CELLTYPE <- "Sertoli"
MATURE_REFERENCE_LABEL <- "Normal mature reference"
IMMATURE_REFERENCE_LABEL <- "Normal immature reference"
set.seed(20260501)

# Color palette adapted from the user-provided reference.
# Pink tones are used for iNOA / stress-shifted states; blue tones for the
# normal maturation axis and AZFa_Del comparator; grey-lavender tones for
# supporting/reference elements.
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
  lavender = "#A497BA"
)

palette_table <- data.frame(
  color_name = names(PAL),
  hex = unlist(PAL, use.names = FALSE),
  recommended_use = c(
    "iNOA main color and high end of molecular heatmap",
    "secondary iNOA / soft pink emphasis",
    "light iNOA fill and low-opacity annotations",
    "diverging heatmap midpoint / pale background",
    "normal mature reference and low end of molecular heatmap",
    "AZFa_Del comparator",
    "normal immature reference and normal trajectory ribbon",
    "optional secondary blue-lavender accent",
    "axes, dashed reference lines, and neutral labels",
    "light neutral outlines",
    "panel strips and pale neutral fills",
    "optional cluster / support accent"
  ),
  stringsAsFactors = FALSE
)
fwrite(palette_table, file.path(TABLE_DIR, "part6_visual_palette.csv"))

find_file <- function(rel_paths, required = TRUE, label = "file") {
  for (rel in rel_paths) {
    f <- file.path(BASE_DIR, rel)
    if (file.exists(f)) return(normalizePath(f, winslash = "/", mustWork = TRUE))
  }
  if (required) {
    stop("Missing ", label, ". Tried: ", paste(file.path(BASE_DIR, rel_paths), collapse = " ; "))
  }
  NA_character_
}

PART3_AGE_SUMMARY <- find_file(
  c(
    "results/03_developmental_reference_construction/tables/all_celltypes_pseudotime_age_summary.csv",
    "part3/tables/all_celltypes_pseudotime_age_summary.csv",
    "part3_multicelltype_pseudotime_results/tables/all_celltypes_pseudotime_age_summary.csv"
  ),
  label = "Part 3 age summary"
)

PART3_VALUES <- find_file(
  c(
    "results/03_developmental_reference_construction/tables/all_celltypes_pseudotime_values.csv",
    "part3/tables/all_celltypes_pseudotime_values.csv",
    "part3_multicelltype_pseudotime_results/tables/all_celltypes_pseudotime_values.csv"
  ),
  required = FALSE,
  label = "Part 3 cell-level pseudotime values"
)

PART3_QUALITY <- find_file(
  c(
    "results/03_developmental_reference_construction/tables/all_celltypes_pseudotime_quality_metrics.csv",
    "part3/tables/all_celltypes_pseudotime_quality_metrics.csv",
    "part3_multicelltype_pseudotime_results/tables/all_celltypes_pseudotime_quality_metrics.csv"
  ),
  required = FALSE,
  label = "Part 3 pseudotime quality metrics"
)

PART4_SAMPLE_SUMMARY <- find_file(
  c(
    "results/04_disease_cell_projection/tables/disease_projection_summary_by_sample.csv",
    "part4/disease_projection_summary_by_sample.csv",
    "part4/tables/disease_projection_summary_by_sample.csv",
    "part4_disease_projection_results/tables/disease_projection_summary_by_sample.csv"
  ),
  label = "Part 4 sample-level projection summary"
)

PART4_GROUP_SUMMARY <- find_file(
  c(
    "results/04_disease_cell_projection/tables/disease_projection_summary_by_group.csv",
    "part4/disease_projection_summary_by_group.csv",
    "part4/tables/disease_projection_summary_by_group.csv",
    "part4_disease_projection_results/tables/disease_projection_summary_by_group.csv"
  ),
  required = FALSE,
  label = "Part 4 group-level projection summary"
)

PART4_PROJECTED_VALUES <- find_file(
  c(
    "results/04_disease_cell_projection/tables/disease_projected_pseudotime_values.csv",
    "part4/disease_projected_pseudotime_values.csv",
    "part4/tables/disease_projected_pseudotime_values.csv",
    "part4_disease_projection_results/tables/disease_projected_pseudotime_values.csv"
  ),
  required = FALSE,
  label = "Part 4 cell-level projected pseudotime values"
)

PART5_SCORE_SUMMARY <- find_file(
  c(
    "results/05_sertoli_molecular_characterization/tables/part5_sertoli_main_gene_set_score_summary.csv",
    "part5/tables/part5_sertoli_main_gene_set_score_summary.csv",
    "part5_sertoli_mechanism_characterization_results/tables/part5_sertoli_main_gene_set_score_summary.csv",
    "part5_sertoli_mechanism_validation_results/tables/part5_sertoli_main_gene_set_score_summary.csv"
  ),
  required = FALSE,
  label = "Part 5 main gene-set score summary"
)

PART5_SCORE_TESTS <- find_file(
  c(
    "results/05_sertoli_molecular_characterization/tables/part5_sertoli_main_gene_set_score_tests_vs_mature_reference.csv",
    "part5/tables/part5_sertoli_main_gene_set_score_tests_vs_mature_reference.csv",
    "part5_sertoli_mechanism_characterization_results/tables/part5_sertoli_main_gene_set_score_tests_vs_mature_reference.csv",
    "part5_sertoli_mechanism_validation_results/tables/part5_sertoli_main_gene_set_score_tests_vs_mature_reference.csv"
  ),
  required = FALSE,
  label = "Part 5 score tests"
)

message("Part 6 base directory: ", BASE_DIR)

# ------------------------------------------------------------------------------
# 1. Helper functions
# ------------------------------------------------------------------------------

to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 2),
      plot.subtitle = element_text(color = PAL$grey_dark, size = base_size),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.35),
      legend.title = element_text(face = "bold"),
      legend.key.size = grid::unit(0.42, "cm"),
      strip.background = element_rect(fill = PAL$grey_light, color = NA),
      strip.text = element_text(face = "bold", color = "black"),
      plot.margin = margin(5.5, 7, 5.5, 5.5)
    )
}

save_fig <- function(name, plot, width, height, section = c("main", "supplementary")) {
  section <- match.arg(section)
  out_dir <- if (section == "main") MAIN_FIG_DIR else SUPP_FIG_DIR
  out_file <- file.path(out_dir, name)
  ggsave(out_file, plot = plot, width = width, height = height, device = cairo_pdf)
  message("Saved figure: ", normalizePath(out_file, winslash = "/", mustWork = FALSE))
  invisible(out_file)
}

save_legacy_fig <- function(name, plot, width, height) {
  out_file <- file.path(FIG_DIR, name)
  ggsave(out_file, plot = plot, width = width, height = height, device = cairo_pdf)
  message("Saved legacy-compatible figure: ", normalizePath(out_file, winslash = "/", mustWork = FALSE))
  invisible(out_file)
}

add_panel_label <- function(plot, label) {
  plot +
    labs(tag = label) +
    theme(
      plot.tag = element_text(face = "bold", size = 12),
      plot.tag.position = c(0.01, 0.99)
    )
}

bootstrap_median_ci <- function(x, n_boot = 1000, conf = 0.95) {
  x <- x[is.finite(x)]
  if (length(x) < 2) {
    return(data.frame(
      median = ifelse(length(x) == 1, x[[1]], NA_real_),
      ci_low = NA_real_,
      ci_high = NA_real_,
      n = length(x)
    ))
  }
  boots <- replicate(n_boot, median(sample(x, length(x), replace = TRUE), na.rm = TRUE))
  alpha <- (1 - conf) / 2
  data.frame(
    median = median(x, na.rm = TRUE),
    ci_low = unname(quantile(boots, alpha, na.rm = TRUE)),
    ci_high = unname(quantile(boots, 1 - alpha, na.rm = TRUE)),
    n = length(x)
  )
}

infer_reference_age <- function(pt, ref_age_summary) {
  ref <- ref_age_summary %>%
    arrange(median_pseudotime) %>%
    distinct(median_pseudotime, .keep_all = TRUE)
  x <- ref$median_pseudotime
  y <- ref$age_numeric
  if (length(x) < 2) return(rep(NA_real_, length(pt)))
  pt_clip <- pmin(pmax(pt, min(x, na.rm = TRUE)), max(x, na.rm = TRUE))
  approx(x = x, y = y, xout = pt_clip, rule = 2, ties = "ordered")$y
}

nearest_reference_age <- function(pt, ref_age_summary) {
  vapply(pt, function(z) {
    if (!is.finite(z)) return(NA_real_)
    ref_age_summary$age_numeric[which.min(abs(ref_age_summary$median_pseudotime - z))]
  }, numeric(1))
}

signed_support <- function(condition) {
  if (isTRUE(condition)) "supports"
  else if (identical(condition, FALSE)) "does_not_support"
  else "not_tested"
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.finite(x), formatC(x, digits = digits, format = "f"), "NA")
}

ensure_columns <- function(df, cols, value = NA_real_) {
  for (cc in cols) {
    if (!cc %in% colnames(df)) df[[cc]] <- value
  }
  df
}

# ------------------------------------------------------------------------------
# 2. Load normal reference and disease projection summaries
# ------------------------------------------------------------------------------

age_summary <- fread(PART3_AGE_SUMMARY) %>%
  filter(.data$cell_type == SERTOLI_CELLTYPE) %>%
  mutate(
    age_numeric = to_num(age_numeric),
    n_cells = to_num(n_cells),
    median_pseudotime = to_num(median_pseudotime),
    q25 = to_num(q25),
    q75 = to_num(q75),
    median_maturation_score = to_num(median_maturation_score)
  ) %>%
  arrange(age_numeric)

if (nrow(age_summary) < 2) {
  stop("Part 3 Sertoli age summary has fewer than two age points.")
}

immature_age <- min(age_summary$age_numeric, na.rm = TRUE)
mature_age <- max(age_summary$age_numeric, na.rm = TRUE)
mature_ref <- age_summary %>% filter(age_numeric == mature_age) %>% slice(1)
immature_ref <- age_summary %>% filter(age_numeric == immature_age) %>% slice(1)
mature_ref_pt <- mature_ref$median_pseudotime[[1]]
immature_ref_pt <- immature_ref$median_pseudotime[[1]]
reference_pt_span <- mature_ref_pt - immature_ref_pt

sample_projection_raw <- fread(PART4_SAMPLE_SUMMARY) %>%
  ensure_columns(c(
    "median_projected_maturation_score",
    "median_maturation_deficit_vs_mature",
    "wasserstein_to_mature_reference"
  ))

sample_projection <- sample_projection_raw %>%
  filter(.data$cell_type == SERTOLI_CELLTYPE) %>%
  mutate(
    n_cells = to_num(n_cells),
    median_projected_pseudotime = to_num(median_projected_pseudotime),
    median_pseudotime_deficit_vs_mature = to_num(median_pseudotime_deficit_vs_mature),
    median_projected_maturation_score = to_num(median_projected_maturation_score),
    median_maturation_deficit_vs_mature = to_num(median_maturation_deficit_vs_mature),
    wasserstein_to_mature_reference = to_num(wasserstein_to_mature_reference),
    disease_group = as.character(disease_group),
    sample_id = as.character(sample_id),
    inferred_normal_equiv_age = infer_reference_age(median_projected_pseudotime, age_summary),
    nearest_normal_age = nearest_reference_age(median_projected_pseudotime, age_summary),
    developmental_age_deficit_vs_mature = mature_age - inferred_normal_equiv_age,
    normalized_maturity_position = ifelse(
      is.finite(reference_pt_span) && abs(reference_pt_span) > .Machine$double.eps,
      (median_projected_pseudotime - immature_ref_pt) / reference_pt_span,
      NA_real_
    ),
    normalized_maturity_position = pmin(pmax(normalized_maturity_position, 0), 1)
  ) %>%
  arrange(disease_group, sample_id)

if (nrow(sample_projection) == 0) {
  stop("No Sertoli rows were found in Part 4 sample-level projection summary.")
}

fwrite(sample_projection, file.path(TABLE_DIR, "part6_sertoli_sample_level_developmental_projection.csv"))

if (!is.na(PART4_GROUP_SUMMARY)) {
  group_projection_raw <- fread(PART4_GROUP_SUMMARY) %>%
    ensure_columns(c(
      "median_projected_maturation_score",
      "median_maturation_deficit_vs_mature",
      "mean_wasserstein_to_mature_reference"
    ))

  group_projection <- group_projection_raw %>%
    filter(.data$cell_type == SERTOLI_CELLTYPE) %>%
    mutate(
      n_samples = to_num(n_samples),
      n_cells = to_num(n_cells),
      median_projected_pseudotime = to_num(median_projected_pseudotime),
      median_pseudotime_deficit_vs_mature = to_num(median_pseudotime_deficit_vs_mature),
      median_projected_maturation_score = to_num(median_projected_maturation_score),
      median_maturation_deficit_vs_mature = to_num(median_maturation_deficit_vs_mature),
      mean_wasserstein_to_mature_reference = to_num(mean_wasserstein_to_mature_reference),
      inferred_normal_equiv_age = infer_reference_age(median_projected_pseudotime, age_summary),
      developmental_age_deficit_vs_mature = mature_age - inferred_normal_equiv_age
    )
  fwrite(group_projection, file.path(TABLE_DIR, "part6_sertoli_group_level_developmental_projection.csv"))
} else {
  group_projection <- data.frame()
}

# ------------------------------------------------------------------------------
# 3. Optional cell-level projection robustness, if Part 4 per-cell table exists
# ------------------------------------------------------------------------------

cell_projection <- NULL
cell_projection_bootstrap <- data.frame()

if (!is.na(PART4_PROJECTED_VALUES)) {
  projected <- fread(PART4_PROJECTED_VALUES) %>%
    ensure_columns("projected_maturation_score")
  needed_cols <- c("cell", "cell_type", "disease_group", "projected_pseudotime")
  if (all(needed_cols %in% colnames(projected))) {
    if (!"sample_id" %in% colnames(projected)) projected$sample_id <- projected$disease_group
    cell_projection <- projected %>%
      filter(.data$cell_type == SERTOLI_CELLTYPE) %>%
      mutate(
        disease_group = as.character(disease_group),
        sample_id = as.character(sample_id),
        projected_pseudotime = to_num(projected_pseudotime),
        projected_maturation_score = to_num(projected_maturation_score)
      ) %>%
      filter(is.finite(projected_pseudotime))

    if (nrow(cell_projection) > 0) {
      cell_projection_bootstrap <- bind_rows(lapply(split(cell_projection, cell_projection$sample_id), function(df) {
        ci <- bootstrap_median_ci(df$projected_pseudotime)
        ci$disease_group <- unique(df$disease_group)[[1]]
        ci$sample_id <- unique(df$sample_id)[[1]]
        ci
      })) %>%
        mutate(
          inferred_normal_equiv_age = infer_reference_age(median, age_summary),
          developmental_age_deficit_vs_mature = mature_age - inferred_normal_equiv_age
        ) %>%
        select(disease_group, sample_id, n, median, ci_low, ci_high,
               inferred_normal_equiv_age, developmental_age_deficit_vs_mature)

      fwrite(cell_projection, file.path(TABLE_DIR, "part6_sertoli_cell_level_projected_pseudotime_values.csv"))
      fwrite(cell_projection_bootstrap, file.path(TABLE_DIR, "part6_sertoli_cell_level_projection_bootstrap_ci.csv"))
    }
  } else {
    warning("Part 4 cell-level table was found but lacks required columns: ",
            paste(setdiff(needed_cols, colnames(projected)), collapse = ", "))
  }
}

# ------------------------------------------------------------------------------
# 4. Figures: normal reference, disease projection, and molecular characterization
# ------------------------------------------------------------------------------

disease_palette <- c(
  "iNOA" = PAL$pink_dark,
  "AZFa_Del" = PAL$blue_mid
)

all_disease_groups <- unique(sample_projection$disease_group)
missing_colors <- setdiff(all_disease_groups, names(disease_palette))
if (length(missing_colors) > 0) {
  disease_palette <- c(disease_palette, setNames(rep(PAL$grey_dark, length(missing_colors)), missing_colors))
}

disease_shape <- c("iNOA" = 24, "AZFa_Del" = 22)
missing_shapes <- setdiff(all_disease_groups, names(disease_shape))
if (length(missing_shapes) > 0) {
  disease_shape <- c(disease_shape, setNames(rep(21, length(missing_shapes)), missing_shapes))
}

p_maturation <- NULL
p_part5 <- NULL
p_density <- NULL
p_bootstrap <- NULL

p_ref_projection <- ggplot(age_summary, aes(x = age_numeric, y = median_pseudotime)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = PAL$blue_light, alpha = 0.22) +
  geom_line(color = PAL$blue_dark, linewidth = 0.95) +
  geom_point(aes(size = n_cells), shape = 21, fill = PAL$grey_light, color = PAL$blue_dark, stroke = 0.55) +
  geom_hline(yintercept = mature_ref_pt, color = PAL$grey_dark, linetype = "dashed", linewidth = 0.32) +
  geom_point(
    data = sample_projection,
    aes(
      x = inferred_normal_equiv_age,
      y = median_projected_pseudotime,
      fill = disease_group,
      shape = disease_group
    ),
    color = "black",
    size = 3.1,
    stroke = 0.35,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = disease_palette, name = "Disease group") +
  scale_shape_manual(values = disease_shape, name = "Disease group") +
  scale_size_area(max_size = 5.2, name = "Normal cells") +
  scale_x_continuous(breaks = sort(unique(age_summary$age_numeric))) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(
    x = "Normal age / disease equivalent age",
    y = "Sertoli developmental pseudotime",
    title = "Disease Sertoli cells project to early normal developmental states",
    subtitle = "Normal median and IQR are shown as the line and band; dashed line marks mature normal reference"
  ) +
  theme_pub(10) +
  theme(legend.position = "right")

if (HAS_GGREPEL) {
  p_ref_projection <- p_ref_projection +
    ggrepel::geom_text_repel(
      data = sample_projection,
      aes(
        x = inferred_normal_equiv_age,
        y = median_projected_pseudotime,
        label = sample_id
      ),
      inherit.aes = FALSE,
      size = 3,
      box.padding = 0.25,
      point.padding = 0.2,
      segment.color = PAL$grey_dark,
      segment.size = 0.2,
      max.overlaps = Inf,
      show.legend = FALSE
    )
}

save_fig("Fig7A_Sertoli_projection_on_normal_development_axis.pdf",
         p_ref_projection, width = 7.4, height = 4.7, section = "main")
save_legacy_fig("Fig7A_sertoli_projection_on_normal_reference.pdf",
                p_ref_projection, width = 7.4, height = 4.7)

p_deficit_data <- sample_projection %>%
  mutate(
    disease_group = factor(disease_group, levels = unique(c("AZFa_Del", "iNOA", disease_group))),
    sample_id = factor(sample_id, levels = sample_id[order(disease_group, median_pseudotime_deficit_vs_mature)])
  )

p_deficit <- p_deficit_data %>%
  ggplot(aes(x = sample_id, y = median_pseudotime_deficit_vs_mature, fill = disease_group)) +
  geom_hline(yintercept = 0, color = PAL$grey_dark, linewidth = 0.35) +
  geom_col(width = 0.68, color = "white", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.2f", median_pseudotime_deficit_vs_mature)),
            vjust = -0.35, size = 3, color = "black") +
  scale_fill_manual(values = disease_palette, name = "Disease group") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.16))) +
  labs(
    x = NULL,
    y = "Pseudotime deficit vs mature reference",
    title = "iNOA Sertoli cells show a consistent maturation deficit",
    subtitle = "Positive values indicate projection earlier than mature normal Sertoli cells"
  ) +
  theme_pub(10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right")

save_fig("Fig7B_iNOA_Sertoli_pseudotime_deficit_by_sample.pdf",
         p_deficit, width = 7.0, height = 4.5, section = "main")
save_legacy_fig("Fig7B_sertoli_pseudotime_deficit_by_sample.pdf",
                p_deficit, width = 7.0, height = 4.5)

if ("median_projected_maturation_score" %in% colnames(sample_projection) &&
    any(is.finite(sample_projection$median_projected_maturation_score))) {
  p_maturation <- ggplot() +
    geom_path(
      data = age_summary,
      aes(x = median_pseudotime, y = median_maturation_score),
      color = PAL$grey_dark,
      linewidth = 0.65,
      alpha = 0.75
    ) +
    geom_point(
      data = age_summary,
      aes(x = median_pseudotime, y = median_maturation_score),
      shape = 21,
      fill = PAL$grey_light,
      color = PAL$grey_dark,
      size = 2.4,
      stroke = 0.4,
      alpha = 0.9
    ) +
    geom_text(
      data = age_summary,
      aes(x = median_pseudotime, y = median_maturation_score, label = paste0(age_numeric, "y")),
      size = 2.6,
      nudge_y = 0.04,
      color = PAL$grey_dark
    ) +
    geom_vline(xintercept = mature_ref_pt, linetype = "dashed", color = PAL$grey_dark, linewidth = 0.35) +
    geom_hline(yintercept = mature_ref$median_maturation_score[[1]], linetype = "dashed", color = PAL$grey_dark, linewidth = 0.35) +
    geom_point(
      data = sample_projection,
    aes(
      x = median_projected_pseudotime,
      y = median_projected_maturation_score,
      fill = disease_group,
      shape = disease_group,
      size = n_cells
      ),
      color = "black",
      stroke = 0.35,
      alpha = 0.95
    ) +
    scale_fill_manual(values = disease_palette, name = "Disease group") +
    scale_shape_manual(values = disease_shape, name = "Disease group") +
    scale_size_area(max_size = 5.2, name = "Disease cells") +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    labs(
      x = "Median projected pseudotime",
      y = "Median projected maturation score",
      title = "Projected pseudotime is coupled to marker-based maturation",
      subtitle = "Grey line shows normal age medians; disease symbols show sample medians"
    ) +
    theme_pub(10) +
    theme(legend.position = "right")

  if (HAS_GGREPEL) {
    p_maturation <- p_maturation +
      ggrepel::geom_text_repel(
        data = sample_projection,
        aes(
          x = median_projected_pseudotime,
          y = median_projected_maturation_score,
          label = sample_id
        ),
        inherit.aes = FALSE,
        size = 3,
        box.padding = 0.25,
        point.padding = 0.2,
        segment.color = PAL$grey_dark,
        segment.size = 0.2,
        max.overlaps = Inf,
        show.legend = FALSE
      )
  }

  save_fig("Fig7C_projected_pseudotime_maturation_score_coupling.pdf",
           p_maturation, width = 7.0, height = 4.8, section = "main")
  save_legacy_fig("Fig7C_sertoli_projected_pseudotime_vs_maturation_score.pdf",
                  p_maturation, width = 7.0, height = 4.8)
}

if (!is.null(cell_projection) && nrow(cell_projection) > 0 && !is.na(PART3_VALUES)) {
  normal_values <- fread(PART3_VALUES) %>%
    filter(.data$cell_type == SERTOLI_CELLTYPE) %>%
    mutate(
      age_numeric = to_num(age_numeric),
      pseudotime_final = to_num(pseudotime_final)
    ) %>%
    filter(age_numeric %in% c(immature_age, mature_age), is.finite(pseudotime_final)) %>%
    transmute(
      group = paste0("Normal age ", age_numeric),
      pseudotime = pseudotime_final,
      source = "Normal reference"
    )

  disease_values <- cell_projection %>%
    transmute(
      group = disease_group,
      pseudotime = projected_pseudotime,
      source = "Disease projection"
    )

  density_df <- bind_rows(normal_values, disease_values) %>%
    mutate(
      group = factor(
        group,
        levels = unique(c(paste0("Normal age ", immature_age), paste0("Normal age ", mature_age), "AZFa_Del", "iNOA", group))
      )
    )

  density_colors <- c(
    setNames(PAL$blue_light, paste0("Normal age ", immature_age)),
    setNames(PAL$blue_dark, paste0("Normal age ", mature_age)),
    "AZFa_Del" = PAL$blue_mid,
    "iNOA" = PAL$pink_dark
  )

  p_density <- ggplot(density_df, aes(x = pseudotime, color = group, fill = group)) +
    geom_density(alpha = 0.14, linewidth = 0.9, adjust = 1.1) +
    geom_vline(xintercept = mature_ref_pt, linetype = "dashed", color = PAL$grey_dark, linewidth = 0.35) +
    scale_color_manual(values = density_colors, name = NULL, na.value = PAL$grey_dark) +
    scale_fill_manual(values = density_colors, name = NULL, na.value = PAL$grey_light) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    labs(
      x = "Projected / reference pseudotime",
      y = "Density",
      title = "Cell-level projected pseudotime distributions",
      subtitle = "Shown as supporting evidence; sample-level medians remain the primary statistical unit"
    ) +
    theme_pub(10) +
    theme(legend.position = "top")

  save_fig("FigS7A_cell_level_projected_pseudotime_density.pdf",
           p_density, width = 7.6, height = 4.6, section = "supplementary")
  save_legacy_fig("Fig7D_sertoli_cell_level_projected_pseudotime_density.pdf",
                  p_density, width = 7.6, height = 4.6)

  if (nrow(cell_projection_bootstrap) > 0) {
    p_bootstrap <- cell_projection_bootstrap %>%
      mutate(
        disease_group = factor(disease_group, levels = unique(c("AZFa_Del", "iNOA", disease_group))),
        sample_id = factor(sample_id, levels = sample_id[order(disease_group, median)])
      ) %>%
      ggplot(aes(x = sample_id, y = median, ymin = ci_low, ymax = ci_high, color = disease_group)) +
      geom_hline(yintercept = mature_ref_pt, color = PAL$grey_dark, linetype = "dashed", linewidth = 0.35) +
      geom_pointrange(linewidth = 0.55, fatten = 2.4) +
      scale_color_manual(values = disease_palette, name = "Disease group") +
      scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
      labs(
        x = NULL,
        y = "Median projected pseudotime with bootstrap CI",
        title = "Cell-level bootstrap intervals support early iNOA projections",
        subtitle = "Dashed line marks mature normal Sertoli median pseudotime"
      ) +
      theme_pub(10) +
      theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right")

    save_fig("FigS7B_cell_level_projected_pseudotime_bootstrap_CI.pdf",
             p_bootstrap, width = 7.0, height = 4.4, section = "supplementary")
  }
}

# Optional Part 5 molecular characterization heatmap.
part5_molecular_summary <- data.frame()

if (!is.na(PART5_SCORE_SUMMARY)) {
  score_summary <- fread(PART5_SCORE_SUMMARY) %>%
    mutate(
      analysis_group = as.character(analysis_group),
      score_set = as.character(score_set),
      median_score = to_num(median_score),
      n_cells = to_num(n_cells)
    )

  molecular_groups <- c(
    IMMATURE_REFERENCE_LABEL,
    MATURE_REFERENCE_LABEL,
    "AZFa_Del Sertoli",
    "iNOA immature-like Sertoli",
    "iNOA mature-like Sertoli"
  )

  part5_molecular_summary <- score_summary %>%
    filter(.data$analysis_group %in% molecular_groups) %>%
    group_by(score_set) %>%
    mutate(
      mature_reference_median = median_score[analysis_group == MATURE_REFERENCE_LABEL][1],
      delta_vs_mature_reference = median_score - mature_reference_median
    ) %>%
    ungroup() %>%
    mutate(
      analysis_group = factor(analysis_group, levels = molecular_groups),
      analysis_group_label = dplyr::recode(
        as.character(analysis_group),
        "Normal immature reference" = "Normal\nimmature",
        "Normal mature reference" = "Normal\nmature",
        "AZFa_Del Sertoli" = "AZFa_Del\nSertoli",
        "iNOA immature-like Sertoli" = "iNOA\nimmature-like",
        "iNOA mature-like Sertoli" = "iNOA\nmature-like"
      ),
      analysis_group_label = factor(
        analysis_group_label,
        levels = c("Normal\nimmature", "Normal\nmature", "AZFa_Del\nSertoli", "iNOA\nimmature-like", "iNOA\nmature-like")
      ),
      score_set_label = dplyr::recode(
        score_set,
        "Stress response" = "Stress\nresponse",
        "Sertoli function" = "Sertoli\nfunction",
        "Senescence" = "Senescence"
      ),
      score_set_label = factor(score_set_label, levels = rev(c("Senescence", "Stress\nresponse", "Sertoli\nfunction")))
    )

  fwrite(part5_molecular_summary, file.path(TABLE_DIR, "part6_part5_molecular_score_delta_summary.csv"))

  p_part5 <- ggplot(part5_molecular_summary, aes(x = analysis_group_label, y = score_set_label, fill = delta_vs_mature_reference)) +
    geom_tile(color = "white", linewidth = 0.35) +
    geom_text(aes(label = sprintf("%.2f", delta_vs_mature_reference)), size = 2.7, color = "black") +
    scale_fill_gradient2(
      low = PAL$blue_dark,
      mid = PAL$pink_pale,
      high = PAL$pink_dark,
      midpoint = 0,
      name = "Delta vs\nmature ref."
    ) +
    labs(
      x = NULL,
      y = NULL,
      title = "Projected states show coordinated molecular rewiring",
      subtitle = "Scores are summarized relative to mature normal Sertoli reference"
    ) +
    theme_pub(9) +
    theme(
      axis.text.x = element_text(size = 8.2, lineheight = 0.92),
      axis.text.y = element_text(face = "bold"),
      legend.position = "right"
    )

  save_fig("Fig7D_molecular_program_characterization_heatmap.pdf",
           p_part5, width = 7.2, height = 4.5, section = "main")
  save_legacy_fig("Fig7E_part5_molecular_characterization_heatmap.pdf",
                  p_part5, width = 7.2, height = 4.5)
}

if (HAS_PATCHWORK && !is.null(p_maturation) && !is.null(p_part5)) {
  p_main_composite <- (
    add_panel_label(p_ref_projection, "A") | add_panel_label(p_deficit, "B")
  ) / (
    add_panel_label(p_maturation, "C") | add_panel_label(p_part5, "D")
  ) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = "Immature-like Sertoli-cell phenotype in iNOA",
      theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0))
    )

  save_fig("Fig7_Sertoli_immature_like_phenotype_main.pdf",
           p_main_composite, width = 14.2, height = 9.2, section = "main")
} else if (HAS_PATCHWORK && !is.null(p_maturation)) {
  p_main_composite <- (
    add_panel_label(p_ref_projection, "A") | add_panel_label(p_deficit, "B")
  ) / add_panel_label(p_maturation, "C") +
    patchwork::plot_annotation(
      title = "Immature-like Sertoli-cell phenotype in iNOA",
      theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0))
    )

  save_fig("Fig7_Sertoli_immature_like_phenotype_main.pdf",
           p_main_composite, width = 13.4, height = 8.8, section = "main")
}

# ------------------------------------------------------------------------------
# 5. Evidence synthesis tables and report
# ------------------------------------------------------------------------------

inoa_samples <- sample_projection %>% filter(disease_group == "iNOA")
azfa_samples <- sample_projection %>% filter(disease_group == "AZFa_Del")

inoa_group <- if (nrow(group_projection) > 0) {
  group_projection %>% filter(disease_group == "iNOA") %>% slice(1)
} else {
  data.frame()
}

azfa_group <- if (nrow(group_projection) > 0) {
  group_projection %>% filter(disease_group == "AZFa_Del") %>% slice(1)
} else {
  data.frame()
}

quality_rho <- NA_real_
raw_quality_rho <- NA_real_
age_ordered_quality_rho <- NA_real_
if (!is.na(PART3_QUALITY)) {
  quality <- fread(PART3_QUALITY) %>%
    filter(.data$cell_type == SERTOLI_CELLTYPE)
  quality_rho <- quality %>%
    filter(metric == "spearman_final_age") %>%
    pull(value) %>%
    to_num() %>%
    .[1]
  raw_quality_rho <- quality %>%
    filter(metric == "spearman_raw_age") %>%
    pull(value) %>%
    to_num() %>%
    .[1]
  age_ordered_quality_rho <- quality %>%
    filter(metric == "spearman_age_ordered_age") %>%
    pull(value) %>%
    to_num() %>%
    .[1]
}

n_normal_ages <- dplyr::n_distinct(age_summary$age_numeric[is.finite(age_summary$age_numeric)])
min_normal_age_cells <- if (nrow(age_summary) > 0) min(age_summary$n_cells, na.rm = TRUE) else NA_real_

all_inoa_lower_than_mature <- nrow(inoa_samples) > 0 &&
  all(inoa_samples$median_projected_pseudotime < mature_ref_pt, na.rm = TRUE)

all_inoa_low_deficit <- nrow(inoa_samples) > 0 &&
  all(inoa_samples$median_pseudotime_deficit_vs_mature > 0, na.rm = TRUE)

inoa_lower_than_azfa <- nrow(inoa_samples) > 0 && nrow(azfa_samples) > 0 &&
  max(inoa_samples$median_projected_pseudotime, na.rm = TRUE) <
  min(azfa_samples$median_projected_pseudotime, na.rm = TRUE)

median_inoa_pt <- if (nrow(inoa_samples) > 0) median(inoa_samples$median_projected_pseudotime, na.rm = TRUE) else NA_real_
median_azfa_pt <- if (nrow(azfa_samples) > 0) median(azfa_samples$median_projected_pseudotime, na.rm = TRUE) else NA_real_
median_inoa_age <- if (nrow(inoa_samples) > 0) median(inoa_samples$inferred_normal_equiv_age, na.rm = TRUE) else NA_real_
median_azfa_age <- if (nrow(azfa_samples) > 0) median(azfa_samples$inferred_normal_equiv_age, na.rm = TRUE) else NA_real_

molecular_available <- nrow(part5_molecular_summary) > 0
molecular_support_note <- if (molecular_available) {
  "Part 5 molecular characterization summary was integrated; inspect Fig7D and part6_part5_molecular_score_delta_summary.csv."
} else {
  "Part 5 molecular characterization summary was not found, so molecular-program support was not integrated in Part 6."
}

evidence_table <- bind_rows(
  data.frame(
    evidence_step = "Raw normal Sertoli trajectory has developmental signal",
    metric = "Spearman(raw DPT pseudotime, age)",
    result = paste0("raw DPT rho=", fmt_num(raw_quality_rho)),
    support = signed_support(is.finite(raw_quality_rho) && raw_quality_rho > 0),
    caveat = "Raw DPT-age association is the non-age-calibrated trajectory evidence and should be reported with the final axis.",
    stringsAsFactors = FALSE
  ),
  data.frame(
    evidence_step = "Final normal Sertoli reference is age-ordered",
    metric = "Spearman(final pseudotime, age)",
    result = paste0("final rho=", fmt_num(quality_rho), "; age-ordered rho=", fmt_num(age_ordered_quality_rho)),
    support = signed_support(is.finite(quality_rho) && quality_rho >= 0.7),
    caveat = "Final pseudotime is age-calibrated in Part 3; use it as the reference axis, not as independent proof of age association.",
    stringsAsFactors = FALSE
  ),
  data.frame(
    evidence_step = "Normal Sertoli reference has usable age coverage",
    metric = "Number of normal ages and minimum cells per age",
    result = paste0(n_normal_ages, " ages; minimum cells per age=", fmt_num(min_normal_age_cells, digits = 0)),
    support = signed_support(n_normal_ages >= 2 && is.finite(min_normal_age_cells) && min_normal_age_cells >= 50),
    caveat = "Coverage supports projection interpretability but does not replace independent disease-sample replication.",
    stringsAsFactors = FALSE
  ),
  data.frame(
    evidence_step = "iNOA Sertoli projects earlier than mature normal reference",
    metric = "Sample-level median projected pseudotime",
    result = paste0(
      "median iNOA=", fmt_num(median_inoa_pt),
      "; mature normal=", fmt_num(mature_ref_pt),
      "; inferred iNOA age=", fmt_num(median_inoa_age), "y"
    ),
    support = signed_support(all_inoa_lower_than_mature && all_inoa_low_deficit),
    caveat = "Use sample-level summaries as primary evidence to avoid cell-level pseudoreplication.",
    stringsAsFactors = FALSE
  ),
  data.frame(
    evidence_step = "iNOA pattern is consistent across available iNOA samples",
    metric = "All iNOA samples have positive pseudotime deficit vs mature reference",
    result = paste0(sum(inoa_samples$median_pseudotime_deficit_vs_mature > 0, na.rm = TRUE),
                    "/", nrow(inoa_samples), " iNOA samples positive"),
    support = signed_support(all_inoa_low_deficit),
    caveat = paste0(
      "There are ", nrow(inoa_samples),
      " iNOA samples; this supports within-dataset consistency but is not definitive population-level proof."
    ),
    stringsAsFactors = FALSE
  ),
  data.frame(
    evidence_step = "Disease-control comparison",
    metric = "iNOA projected pseudotime lower than AZFa_Del",
    result = paste0(
      "median iNOA=", fmt_num(median_inoa_pt),
      "; AZFa_Del=", fmt_num(median_azfa_pt),
      "; inferred AZFa_Del age=", fmt_num(median_azfa_age), "y"
    ),
    support = signed_support(inoa_lower_than_azfa),
    caveat = paste0(
      "AZFa_Del has ", nrow(azfa_samples),
      " sample(s), so this should be framed as a comparator, not a formal disease-general conclusion."
    ),
    stringsAsFactors = FALSE
  ),
  data.frame(
    evidence_step = "Molecular characterization supports phenotype interpretation",
    metric = "Part 5 stress/senescence/function scores",
    result = molecular_support_note,
    support = ifelse(molecular_available, "supports_if_direction_matches", "not_tested"),
    caveat = "Molecular-program support should preferably emphasize genes not used to orient Part 3 pseudotime.",
    stringsAsFactors = FALSE
  ),
  data.frame(
    evidence_step = "Causal interpretation boundary",
    metric = "Observational single-cell association",
    result = "Supports an immature-like Sertoli-cell phenotype model; does not prove that Sertoli arrest causes iNOA.",
    support = "interpretation_boundary",
    caveat = "Causality would require perturbation, longitudinal, spatial, or independent validation evidence.",
    stringsAsFactors = FALSE
  )
)

fwrite(evidence_table, file.path(TABLE_DIR, "part6_developmental_arrest_evidence_chain.csv"))

manifest <- data.frame(
  item = c(
    "Main Figure 7 composite",
    "Main Figure 7A",
    "Main Figure 7B",
    "Main Figure 7C",
    "Main Figure 7D",
    "Supplementary Figure S7A",
    "Supplementary Figure S7B",
    "Supplementary Table S31",
    "Supplementary Table S32",
    "Supplementary Table S33",
    "Report"
  ),
  file = c(
    "figures/main/Fig7_Sertoli_immature_like_phenotype_main.pdf",
    "figures/main/Fig7A_Sertoli_projection_on_normal_development_axis.pdf",
    "figures/main/Fig7B_iNOA_Sertoli_pseudotime_deficit_by_sample.pdf",
    "figures/main/Fig7C_projected_pseudotime_maturation_score_coupling.pdf",
    "figures/main/Fig7D_molecular_program_characterization_heatmap.pdf",
    "figures/supplementary/FigS7A_cell_level_projected_pseudotime_density.pdf",
    "figures/supplementary/FigS7B_cell_level_projected_pseudotime_bootstrap_CI.pdf",
    "tables/part6_sertoli_sample_level_developmental_projection.csv",
    "tables/part6_sertoli_group_level_developmental_projection.csv",
    "tables/part6_developmental_arrest_evidence_chain.csv",
    "reports/part6_developmental_arrest_logic_chain_report.txt"
  ),
  role = c(
    "Recommended manuscript-ready main figure integrating projection, sample consistency, maturation coupling, and molecular characterization.",
    "Places disease Sertoli sample medians onto the normal Sertoli maturation axis.",
    "Shows sample-level pseudotime deficit relative to mature normal Sertoli cells.",
    "Shows whether projected pseudotime and marker-based maturation score move together.",
    "Links projected states to molecular programs from Part 5.",
    "Cell-level density support; use as supplementary because cells are not independent samples.",
    "Cell-level bootstrap interval support; use as supplementary because bootstrap precision reflects cell count.",
    "Sample-level developmental projection table.",
    "Group-level developmental projection table.",
    "Evidence-chain summary with caveats.",
    "Plain-language interpretation for manuscript writing."
  ),
  stringsAsFactors = FALSE
)
fwrite(manifest, file.path(TABLE_DIR, "part6_figure_table_manifest.csv"))

placement <- data.frame(
  recommended_placement = c(
    "Main Figure 7",
    "Main Figure 7A",
    "Main Figure 7B",
    "Main Figure 7C",
    "Main Figure 7D",
    "Supplementary Figure S7A",
    "Supplementary Figure S7B",
    "Supplementary Table S31",
    "Supplementary Table S32",
    "Supplementary Table S33"
  ),
  file = c(
    "figures/main/Fig7_Sertoli_immature_like_phenotype_main.pdf",
    "figures/main/Fig7A_Sertoli_projection_on_normal_development_axis.pdf",
    "figures/main/Fig7B_iNOA_Sertoli_pseudotime_deficit_by_sample.pdf",
    "figures/main/Fig7C_projected_pseudotime_maturation_score_coupling.pdf",
    "figures/main/Fig7D_molecular_program_characterization_heatmap.pdf",
    "figures/supplementary/FigS7A_cell_level_projected_pseudotime_density.pdf",
    "figures/supplementary/FigS7B_cell_level_projected_pseudotime_bootstrap_CI.pdf",
    "tables/part6_sertoli_sample_level_developmental_projection.csv",
    "tables/part6_sertoli_group_level_developmental_projection.csv",
    "tables/part6_developmental_arrest_evidence_chain.csv"
  ),
  priority = c(
    "Use",
    "Use as source panel if not using composite",
    "Use as source panel if not using composite",
    "Use as source panel if not using composite",
    "Use as source panel if not using composite",
    "Supplement",
    "Supplement",
    "Supplement",
    "Supplement",
    "Supplement"
  ),
  reason = c(
    "Most concise main evidence chain for the Sertoli developmental-arrest model.",
    "Primary visual evidence that disease Sertoli states map onto early normal developmental age.",
    "Sample-level consistency; reduces concern about cell-level pseudoreplication.",
    "Links projected pseudotime with a marker-based maturation score.",
    "Adds molecular mechanism support through Part 5 stress/senescence/function programs.",
    "Useful visual distribution support, but cells should not be treated as independent biological replicates.",
    "Useful robustness support for cell-level projected medians.",
    "Gives the sample-level values underlying Fig7A-B.",
    "Gives group-level summaries for text reporting.",
    "Gives the final logic-chain claims and caveats."
  ),
  stringsAsFactors = FALSE
)
fwrite(placement, file.path(TABLE_DIR, "part6_figure_placement_recommendations.csv"))

filename_map <- data.frame(
  old_name = c(
    "Fig7A_sertoli_projection_on_normal_reference.pdf",
    "Fig7B_sertoli_pseudotime_deficit_by_sample.pdf",
    "Fig7C_sertoli_projected_pseudotime_vs_maturation_score.pdf",
    "Fig7D_sertoli_cell_level_projected_pseudotime_density.pdf",
    "Fig7E_part5_molecular_characterization_heatmap.pdf"
  ),
  optimized_name = c(
    "figures/main/Fig7A_Sertoli_projection_on_normal_development_axis.pdf",
    "figures/main/Fig7B_iNOA_Sertoli_pseudotime_deficit_by_sample.pdf",
    "figures/main/Fig7C_projected_pseudotime_maturation_score_coupling.pdf",
    "figures/supplementary/FigS7A_cell_level_projected_pseudotime_density.pdf",
    "figures/main/Fig7D_molecular_program_characterization_heatmap.pdf"
  ),
  placement_change = c(
    "Keep as Main Fig7A",
    "Keep as Main Fig7B",
    "Keep as Main Fig7C",
    "Move to Supplementary FigS7A",
    "Rename to Main Fig7D"
  ),
  stringsAsFactors = FALSE
)
fwrite(filename_map, file.path(TABLE_DIR, "part6_optimized_filename_map.csv"))

recommended_sentence <- paste(
  "Recommended wording:",
  "iNOA Sertoli cells are shifted toward an early projected developmental state",
  "relative to the mature normal Sertoli reference, consistent with a",
  "immature-like Sertoli-cell state.",
  "Because these data are observational, the result should be framed as a",
  "developmental-reference model rather than direct causal proof."
)

report_lines <- c(
  "Part 6 Sertoli evidence synthesis",
  "=================================",
  "",
  paste0("Base directory: ", BASE_DIR),
  paste0("Normal mature Sertoli age: ", mature_age),
  paste0("Normal mature Sertoli median pseudotime: ", fmt_num(mature_ref_pt)),
  paste0("iNOA sample-level median projected pseudotime: ", fmt_num(median_inoa_pt)),
  paste0("iNOA inferred normal-equivalent age: ", fmt_num(median_inoa_age), " years"),
  paste0("AZFa_Del sample-level median projected pseudotime: ", fmt_num(median_azfa_pt)),
  paste0("AZFa_Del inferred normal-equivalent age: ", fmt_num(median_azfa_age), " years"),
  "",
  "Evidence-chain summary:",
  paste0("- ", evidence_table$evidence_step, ": ", evidence_table$result, " [", evidence_table$support, "]"),
  "",
  "Interpretation:",
  recommended_sentence,
  "",
  "Recommended figure placement:",
  "- Main Figure 7: figures/main/Fig7_Sertoli_immature_like_phenotype_main.pdf",
  "- Main Figure 7A-D source panels are saved separately in figures/main/.",
  "- Supplementary Figure S7A-B: cell-level density and bootstrap CI, saved in figures/supplementary/.",
  "- See tables/part6_figure_placement_recommendations.csv and tables/part6_optimized_filename_map.csv.",
  "",
  "Main caveats:",
  "- Part 3 final pseudotime is age-calibrated; report raw DPT-age rho alongside the final rho.",
  "- Use sample-level disease summaries as the primary statistical unit.",
  "- With one AZFa_Del sample, the disease-control comparison is supportive but not a formal population-level test.",
  "- Causal wording requires perturbation, longitudinal, spatial, or independent cohort validation.",
  "",
  "Generated output directory:",
  normalizePath(OUTDIR, winslash = "/", mustWork = FALSE)
)

writeLines(report_lines, file.path(REPORT_DIR, "part6_developmental_arrest_logic_chain_report.txt"))

message("Part 6 evidence synthesis finished: ", normalizePath(OUTDIR, winslash = "/"))

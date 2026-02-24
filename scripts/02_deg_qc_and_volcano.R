#!/usr/bin/env Rscript
# ══════════════════════════════════════════════════════════════════════════════
# 02_deg_qc_and_volcano.R
# SBMA RNAseq Generalized Pipeline — DEG QC Diagnostics & Volcano Plot
#
# Generates quality-control histograms, an MA-like plot, a threshold
# sensitivity table, and a labeled volcano plot for a single contrast.
#
# Usage:
#   Rscript scripts/02_deg_qc_and_volcano.R \
#     --config config/config.yaml \
#     --dataset SBMA_iPSC_Q51 \
#     --contrast disease_vs_control_vehicle
# ══════════════════════════════════════════════════════════════════════════════

# ── Source pipeline libraries ────────────────────────────────────────────────
script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
project_root <- dirname(script_dir)
source(file.path(project_root, "lib", "config_utils.R"))
source(file.path(project_root, "lib", "data_utils.R"))
source(file.path(project_root, "lib", "plot_utils.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(cowplot)
  library(readr)
})

# ── Parse CLI arguments ─────────────────────────────────────────────────────
args <- parse_pipeline_args()

# ── Load config and contrast info ────────────────────────────────────────────
cfg      <- load_config(args$config)
contrast <- get_contrast(cfg, args$dataset, args$contrast)
out_dir  <- get_output_dir(cfg, args$dataset, args$contrast)

# ── Resolve colors and theme ─────────────────────────────────────────────────
colors    <- if (!is.null(cfg$colors)) cfg$colors else default_plot_colors()
base_size <- if (!is.null(cfg$plots$base_size)) cfg$plots$base_size else 14
pub_theme <- get_theme_publication(base_size)
plot_dpi  <- if (!is.null(cfg$plots$dpi)) cfg$plots$dpi else 300
plot_fmts <- if (!is.null(cfg$plots$format)) cfg$plots$format else c("pdf", "png")

# Thresholds from config
lfc_thresh  <- cfg$thresholds$lfc
padj_thresh <- cfg$thresholds$padj

cat("================================================================\n")
cat(sprintf("  02 DEG QC & Volcano — %s\n", contrast$label))
cat(sprintf("  Dataset : %s | Contrast: %s\n", args$dataset, args$contrast))
cat("================================================================\n\n")

# ── Load prepped DEG data ────────────────────────────────────────────────────
deg_rds <- file.path(out_dir, "data", "deg_all.rds")
if (!file.exists(deg_rds)) {
  stop(sprintf("Prepped DEG file not found: %s\nRun script 01 first.", deg_rds),
       call. = FALSE)
}

deg <- readRDS(deg_rds)
cat(sprintf("Loaded %d genes from %s\n\n", nrow(deg), deg_rds))

# ══════════════════════════════════════════════════════════════════════════════
# 1. QC HISTOGRAMS (combined 3-panel)
# ══════════════════════════════════════════════════════════════════════════════

# -- (a) P-value histogram ---------------------------------------------------
p_pval <- ggplot(deg, aes(x = pvalue)) +
  geom_histogram(bins = 50, fill = colors$down, colour = "white",
                 linewidth = 0.2, alpha = 0.85) +
  labs(
    title = "P-value Distribution",
    x     = "P-value",
    y     = "Count"
  ) +
  pub_theme

# -- (b) Adjusted p-value histogram ------------------------------------------
p_padj <- ggplot(deg, aes(x = padj)) +
  geom_histogram(bins = 50, fill = colors$up, colour = "white",
                 linewidth = 0.2, alpha = 0.85) +
  labs(
    title = "Adjusted P-value Distribution",
    x     = "Adjusted P-value",
    y     = "Count"
  ) +
  pub_theme

# -- (c) Log2FC distribution -------------------------------------------------
p_lfc <- ggplot(deg, aes(x = log2fc)) +
  geom_histogram(bins = 50, fill = "grey50", colour = "white",
                 linewidth = 0.2, alpha = 0.85) +
  geom_vline(xintercept = c(-lfc_thresh, lfc_thresh),
             linetype = "dashed", colour = colors$up, linewidth = 0.6) +
  labs(
    title = expression(log[2]~"FC Distribution"),
    x     = expression(log[2]~"Fold Change"),
    y     = "Count"
  ) +
  pub_theme

# -- Combine into a single 3-panel figure ------------------------------------
qc_hist_panel <- cowplot::plot_grid(
  p_pval, p_padj, p_lfc,
  ncol   = 3,
  labels = c("A", "B", "C"),
  align  = "h",
  axis   = "tb"
)

save_plot(
  plot    = qc_hist_panel,
  path    = file.path(out_dir, "plots", "qc", "qc_histograms"),
  width   = 15,
  height  = 5,
  dpi     = plot_dpi,
  formats = plot_fmts
)

# ══════════════════════════════════════════════════════════════════════════════
# 2. MA-LIKE PLOT (mean expression proxy vs log2FC)
# ══════════════════════════════════════════════════════════════════════════════

# Use baseMean / AveExpr if present; otherwise fall back to -log10(pvalue) as
# a proxy for expression level (genes with more counts yield smaller p-values).
mean_expr_col <- NULL
for (candidate in c("baseMean", "AveExpr", "logCPM", "mean_expression")) {
  if (candidate %in% colnames(deg)) {
    mean_expr_col <- candidate
    break
  }
}

if (!is.null(mean_expr_col)) {
  ma_df <- deg %>%
    mutate(mean_expr = .data[[mean_expr_col]])
  ma_xlab <- mean_expr_col
} else {
  # Proxy: -log10(pvalue) correlates with expression level
  ma_df <- deg %>%
    mutate(mean_expr = -log10(pmax(pvalue, 1e-300)))
  ma_xlab <- expression(-log[10]~"(P-value) [expression proxy]")
}

# Classify for coloring if not already present
if (!"diffexpressed" %in% colnames(ma_df)) {
  ma_df <- classify_deg(ma_df, lfc_thresh = lfc_thresh, padj_thresh = padj_thresh)
}

color_map_ma <- c(
  "Upregulated"   = colors$up,
  "Downregulated" = colors$down,
  "NS"            = colors$ns
)

p_ma <- ggplot(ma_df, aes(x = mean_expr, y = log2fc,
                           colour = diffexpressed)) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_hline(yintercept = c(-lfc_thresh, lfc_thresh),
             linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.3) +
  scale_colour_manual(values = color_map_ma, name = NULL) +
  labs(
    title = sprintf("MA-like Plot — %s", contrast$label),
    x     = ma_xlab,
    y     = expression(log[2]~"Fold Change")
  ) +
  pub_theme +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1)))

save_plot(
  plot    = p_ma,
  path    = file.path(out_dir, "plots", "qc", "ma_plot"),
  width   = 8,
  height  = 6,
  dpi     = plot_dpi,
  formats = plot_fmts
)

# ══════════════════════════════════════════════════════════════════════════════
# 3. THRESHOLD SENSITIVITY TABLE
# ══════════════════════════════════════════════════════════════════════════════

lfc_grid  <- c(0.5, 1.0, 1.5, 2.0)
padj_grid <- c(0.01, 0.05, 0.1)

threshold_summary <- expand.grid(
  lfc_threshold  = lfc_grid,
  padj_threshold = padj_grid,
  stringsAsFactors = FALSE
) %>%
  rowwise() %>%
  mutate(
    n_up   = sum(deg$log2fc >  lfc_threshold & deg$padj < padj_threshold,
                 na.rm = TRUE),
    n_down = sum(deg$log2fc < -lfc_threshold & deg$padj < padj_threshold,
                 na.rm = TRUE),
    n_total = n_up + n_down
  ) %>%
  ungroup() %>%
  arrange(padj_threshold, lfc_threshold)

# Save table
tables_dir <- file.path(out_dir, "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
threshold_csv <- file.path(tables_dir, "deg_threshold_summary.csv")
readr::write_csv(threshold_summary, threshold_csv)
cat(sprintf("Saved threshold summary: %s\n", threshold_csv))

# Print to console
cat("\n--- Threshold Sensitivity Summary ---\n")
print(as.data.frame(threshold_summary), row.names = FALSE)
cat("\n")

# ══════════════════════════════════════════════════════════════════════════════
# 4. VOLCANO PLOT
# ══════════════════════════════════════════════════════════════════════════════

# Ensure classification with config thresholds
deg_classified <- classify_deg(deg, lfc_thresh = lfc_thresh,
                               padj_thresh = padj_thresh)

p_volcano <- make_volcano(
  df          = deg_classified,
  title       = sprintf("Volcano — %s", contrast$label),
  lfc_thresh  = lfc_thresh,
  padj_thresh = padj_thresh,
  colors      = colors,
  top_n_label = 15
)

save_plot(
  plot    = p_volcano,
  path    = file.path(out_dir, "plots", "volcano", "volcano_plot"),
  width   = 10,
  height  = 8,
  dpi     = plot_dpi,
  formats = plot_fmts
)

# ══════════════════════════════════════════════════════════════════════════════
# 5. CONSOLE SUMMARY
# ══════════════════════════════════════════════════════════════════════════════

n_up   <- sum(deg_classified$diffexpressed == "Upregulated",   na.rm = TRUE)
n_down <- sum(deg_classified$diffexpressed == "Downregulated", na.rm = TRUE)
n_ns   <- sum(deg_classified$diffexpressed == "NS",            na.rm = TRUE)

cat("================================================================\n")
cat("  DEG Summary (config thresholds)\n")
cat(sprintf("  |log2FC| > %g  &  padj < %g\n", lfc_thresh, padj_thresh))
cat("================================================================\n")
cat(sprintf("  Total genes tested : %d\n", nrow(deg_classified)))
cat(sprintf("  Upregulated        : %d\n", n_up))
cat(sprintf("  Downregulated      : %d\n", n_down))
cat(sprintf("  Not Significant    : %d\n", n_ns))
cat(sprintf("  Genes with NA padj : %d\n", sum(is.na(deg$padj))))
cat(sprintf("  Genes with NA pval : %d\n", sum(is.na(deg$pvalue))))
cat("================================================================\n")
cat(sprintf("  Median |log2FC| (all)  : %.3f\n",
            median(abs(deg$log2fc), na.rm = TRUE)))
cat(sprintf("  Median |log2FC| (sig)  : %.3f\n",
            median(abs(deg_classified$log2fc[deg_classified$diffexpressed != "NS"]),
                   na.rm = TRUE)))
cat("================================================================\n\n")

cat("Script 02 completed successfully.\n")

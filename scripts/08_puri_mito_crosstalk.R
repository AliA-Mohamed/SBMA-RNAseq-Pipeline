#!/usr/bin/env Rscript
# ==============================================================================
# 08_puri_mito_crosstalk.R
# Mitochondrial-Purinergic Crosstalk Analysis
#
# Integrates mitochondrial and purinergic gene-level results to identify
# overlap genes, build a 28-gene crosstalk gene set across 7 functional
# groups, generate heatmap and forest plot visualizations, extract key
# gene LFCs, and produce a mechanistic narrative summary.
#
# Usage:
#   Rscript scripts/08_puri_mito_crosstalk.R --config config/config.yaml \
#       --dataset SBMA_iPSC_Q51 --contrast disease_vs_control_vehicle
# ==============================================================================

# ---- Resolve project root and source utilities -------------------------------

script_dir   <- dirname(sub("--file=", "",
                             grep("--file=", commandArgs(trailingOnly = FALSE),
                                  value = TRUE)))
project_root <- dirname(script_dir)
setwd(project_root)

source(file.path(project_root, "lib", "config_utils.R"))
source(file.path(project_root, "lib", "data_utils.R"))
source(file.path(project_root, "lib", "plot_utils.R"))
source(file.path(project_root, "lib", "gene_lists.R"))
source(file.path(project_root, "lib", "crosstalk_utils.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# ---- Parse CLI arguments and load config -------------------------------------

args <- parse_pipeline_args()

cfg <- load_config(args$config)

contrast_cfg <- get_contrast(cfg, args$dataset, args$contrast)
out_dir      <- get_output_dir(cfg, args$dataset, args$contrast)

ds <- args$dataset
ct <- args$contrast

message("=== 08_puri_mito_crosstalk.R ===")
message(sprintf("Dataset:  %s", ds))
message(sprintf("Contrast: %s (%s)", ct, contrast_cfg$label))

# ---- Derive paths ------------------------------------------------------------

data_dir   <- file.path(out_dir, "data")
tables_dir <- file.path(out_dir, "tables")
plots_dir  <- file.path(out_dir, "plots", "mitochondrial")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir,  recursive = TRUE, showWarnings = FALSE)

# ---- Load DEG data -----------------------------------------------------------

deg_all_path <- file.path(data_dir, "deg_all.rds")

if (!file.exists(deg_all_path)) {
  stop(sprintf("deg_all.rds not found: %s\nRun the prep step first.",
               deg_all_path), call. = FALSE)
}

deg_all <- readRDS(deg_all_path)
message(sprintf("Loaded DEG table: %d genes.", nrow(deg_all)))

# ---- Load purinergic and mitochondrial summaries -----------------------------

puri_summary_path <- file.path(tables_dir, "purinergic_gene_summary.csv")
mito_summary_path <- file.path(tables_dir, "mitochondrial_gene_summary.csv")

if (!file.exists(puri_summary_path)) {
  stop(sprintf("Purinergic summary not found: %s\nRun the purinergic analysis first.",
               puri_summary_path), call. = FALSE)
}
if (!file.exists(mito_summary_path)) {
  stop(sprintf("Mitochondrial summary not found: %s\nRun the mitochondrial analysis first.",
               mito_summary_path), call. = FALSE)
}

puri_df <- readr::read_csv(puri_summary_path, show_col_types = FALSE)
mito_df <- readr::read_csv(mito_summary_path, show_col_types = FALSE)

message(sprintf("Loaded purinergic summary: %d genes.", nrow(puri_df)))
message(sprintf("Loaded mitochondrial summary: %d genes.", nrow(mito_df)))

# ==============================================================================
# 1. Overlap Analysis
# ==============================================================================

message("\n--- Overlap Analysis ---")

overlap_df <- compute_overlap(puri_df, mito_df)

message(sprintf("Overlap genes found: %d", nrow(overlap_df)))

if (nrow(overlap_df) > 0) {
  message(sprintf("Overlap genes: %s",
                  paste(overlap_df$gene_symbol, collapse = ", ")))
}

overlap_out <- file.path(tables_dir, "mito_purinergic_overlap_genes.csv")
readr::write_csv(overlap_df, overlap_out)
message(sprintf("Saved: %s", overlap_out))

# ==============================================================================
# 2. Crosstalk Gene Set Analysis
# ==============================================================================

message("\n--- Crosstalk Gene Set Analysis ---")

crosstalk_genes <- get_crosstalk_genes()
message(sprintf("Crosstalk gene set: %d genes across %d functional groups.",
                nrow(crosstalk_genes),
                length(unique(crosstalk_genes$function_group))))

crosstalk_df <- build_crosstalk_df(deg_all, crosstalk_genes)

n_found   <- sum(!is.na(crosstalk_df$log2fc))
n_sig     <- sum(crosstalk_df$sig, na.rm = TRUE)
message(sprintf("Matched in DEG table: %d / %d genes.", n_found, nrow(crosstalk_df)))
message(sprintf("Significant (padj < 0.05): %d genes.", n_sig))

crosstalk_out <- file.path(tables_dir, "mito_purinergic_crosstalk_genes.csv")
readr::write_csv(crosstalk_df, crosstalk_out)
message(sprintf("Saved: %s", crosstalk_out))

# ==============================================================================
# 3. Crosstalk Summary by Functional Group
# ==============================================================================

message("\n--- Crosstalk Group Summary ---")

group_summary <- summarize_crosstalk_by_group(crosstalk_df)

group_summary_out <- file.path(tables_dir, "crosstalk_group_summary.csv")
readr::write_csv(group_summary, group_summary_out)
message(sprintf("Saved: %s", group_summary_out))

# ==============================================================================
# 4. Crosstalk Heatmap
# ==============================================================================

message("\n--- Crosstalk Heatmap ---")

# Prepare heatmap data: genes with available LFC values
heatmap_df <- crosstalk_df %>%
  filter(!is.na(log2fc)) %>%
  arrange(function_group, gene_symbol)

if (nrow(heatmap_df) > 0) {
  # Rename LFC column to the contrast label for display
  lfc_col_name <- contrast_cfg$label
  heatmap_df[[lfc_col_name]] <- heatmap_df$log2fc

  hm <- make_heatmap(
    gene_df      = heatmap_df,
    value_cols   = lfc_col_name,
    category_col = "function_group",
    title        = sprintf("Mito-Purinergic Crosstalk — %s", contrast_cfg$label),
    colors       = default_plot_colors()
  )

  hm_height <- max(6, 1.5 + nrow(heatmap_df) * 0.3)

  hm_path <- file.path(plots_dir, "crosstalk_heatmap")
  save_plot(hm, hm_path,
            width   = 8,
            height  = hm_height,
            dpi     = cfg$plots$dpi,
            formats = cfg$plots$format)
} else {
  message("No genes with LFC data available for heatmap — skipping.")
}

# ==============================================================================
# 5. Crosstalk Forest Plot
# ==============================================================================

message("\n--- Crosstalk Forest Plot ---")

forest_df <- crosstalk_df %>%
  filter(!is.na(log2fc)) %>%
  mutate(
    significance = case_when(
      sig & direction == "up"   ~ "Significant Up",
      sig & direction == "down" ~ "Significant Down",
      TRUE                      ~ "Not Significant"
    ),
    significance = factor(significance,
                          levels = c("Significant Down", "Not Significant",
                                     "Significant Up"))
  ) %>%
  arrange(function_group, desc(log2fc))

if (nrow(forest_df) > 0) {
  # Add contrast label as the condition column for the forest plot
  lfc_col_name <- contrast_cfg$label
  forest_df[[lfc_col_name]] <- forest_df$log2fc

  fp <- make_forest_plot(
    gene_df      = forest_df,
    value_cols   = lfc_col_name,
    category_col = "function_group",
    sig_col      = "significance",
    colors       = default_plot_colors()
  )

  fp <- fp +
    labs(title = sprintf("Mito-Purinergic Crosstalk — %s", contrast_cfg$label))

  fp_height <- max(6, 1.5 + nrow(forest_df) * 0.3)

  fp_path <- file.path(plots_dir, "crosstalk_forest_plot")
  save_plot(fp, fp_path,
            width   = 10,
            height  = fp_height,
            dpi     = cfg$plots$dpi,
            formats = cfg$plots$format)
} else {
  message("No genes with LFC data available for forest plot — skipping.")
}

# ==============================================================================
# 6. Mechanism Summary
# ==============================================================================

message("\n--- Mechanism Summary ---")

mechanism_text <- generate_mechanism_summary(crosstalk_df, contrast_cfg$label)

mechanism_out <- file.path(tables_dir, "crosstalk_mechanism_summary.txt")
writeLines(mechanism_text, mechanism_out)
message(sprintf("Saved: %s", mechanism_out))

# ==============================================================================
# 7. Key Gene LFC Extraction
# ==============================================================================

message("\n--- Key Gene LFC Extraction ---")

key_genes <- c(
  "P2RX7",   # Purinergic receptor
  "NLRP3",   # Inflammasome
  "CASP1",   # Caspase-1
  "IL1B",    # Interleukin-1 beta
  "SOD2",    # Mito ROS defense
  "VDAC1",   # Mito/purinergic overlap
  "PANX1",   # ATP release
  "ATP5F1A"  # ATP production
)

key_gene_df <- tibble::tibble(gene_symbol = key_genes) %>%
  left_join(
    deg_all %>%
      filter(!is.na(gene_symbol), gene_symbol != "") %>%
      arrange(padj) %>%
      distinct(gene_symbol, .keep_all = TRUE) %>%
      dplyr::select(gene_symbol, log2fc, pvalue, padj),
    by = "gene_symbol"
  ) %>%
  mutate(
    sig = if_else(!is.na(padj) & padj < 0.05, TRUE, FALSE),
    direction = case_when(
      is.na(log2fc)     ~ NA_character_,
      sig & log2fc > 0  ~ "up",
      sig & log2fc < 0  ~ "down",
      TRUE              ~ "ns"
    )
  )

key_gene_out <- file.path(tables_dir, "crosstalk_key_genes.csv")
readr::write_csv(key_gene_df, key_gene_out)
message(sprintf("Saved: %s", key_gene_out))

# ==============================================================================
# Summary
# ==============================================================================

message("\n========== Mito-Purinergic Crosstalk Summary ==========")
message(sprintf("Contrast: %s (%s)", ct, contrast_cfg$label))
message(sprintf("DEG table:          %d genes", nrow(deg_all)))
message(sprintf("Purinergic genes:   %d", nrow(puri_df)))
message(sprintf("Mitochondrial genes:%d", nrow(mito_df)))
message(paste(rep("-", 52), collapse = ""))

message(sprintf("Overlap genes:      %d", nrow(overlap_df)))
if (nrow(overlap_df) > 0) {
  message(sprintf("  Genes: %s", paste(overlap_df$gene_symbol, collapse = ", ")))
}

message(paste(rep("-", 52), collapse = ""))
message(sprintf("Crosstalk genes:    %d total, %d matched, %d significant",
                nrow(crosstalk_df), n_found, n_sig))

message("\nGroup-level breakdown:")
message(sprintf("  %-30s %5s %5s %5s %8s",
                "Function Group", "Total", "Sig", "Up", "Mean LFC"))
message(paste(rep("-", 60), collapse = ""))

for (i in seq_len(nrow(group_summary))) {
  row <- group_summary[i, ]
  lfc_str <- if (is.nan(row$mean_lfc)) "NA" else sprintf("%.3f", row$mean_lfc)
  message(sprintf("  %-30s %5d %5d %5d %8s",
                  row$function_group,
                  row$n_genes,
                  row$n_sig,
                  row$n_up,
                  lfc_str))
}

message(paste(rep("-", 60), collapse = ""))

message("\nKey gene LFCs:")
message(sprintf("  %-10s %8s %10s %s",
                "Gene", "LFC", "padj", "Direction"))
message(paste(rep("-", 40), collapse = ""))

for (i in seq_len(nrow(key_gene_df))) {
  row <- key_gene_df[i, ]
  lfc_str  <- if (is.na(row$log2fc)) "NA" else sprintf("%.3f", row$log2fc)
  padj_str <- if (is.na(row$padj))   "NA" else sprintf("%.2e", row$padj)
  dir_str  <- if (is.na(row$direction)) "not found" else row$direction
  message(sprintf("  %-10s %8s %10s %s",
                  row$gene_symbol, lfc_str, padj_str, dir_str))
}

message(paste(rep("-", 40), collapse = ""))

message("\nMechanism summary:")
message(mechanism_text)

message("\n====================================================")
message("08_puri_mito_crosstalk.R completed successfully.")

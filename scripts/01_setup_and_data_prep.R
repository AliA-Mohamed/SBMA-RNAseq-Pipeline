#!/usr/bin/env Rscript
# ══════════════════════════════════════════════════════════════════════════════
# 01_setup_and_data_prep.R — Load config, prep DEG data for one contrast
# ══════════════════════════════════════════════════════════════════════════════
#
# Usage:
#   Rscript scripts/01_setup_and_data_prep.R \
#     --config config/config.yaml \
#     --dataset SBMA_iPSC_Q51 \
#     --contrast disease_vs_control_vehicle
#
# Description:
#   Reads a pipeline config.yaml, loads the DEG results file for a single
#   dataset/contrast pair, classifies genes, builds ranked lists (if
#   full_ranked is enabled), and saves all prepared objects to the output
#   directory for use by downstream pipeline steps.
#
# Outputs (under results/<dataset>/<contrast>/data/):
#   deg_all.rds         — Full classified DEG tibble
#   deg_splits.rds      — Named list of gene-symbol vectors (up/down/all_sig/background)
#   ranked_symbol.rds   — Named numeric vector ranked by symbol (only if full_ranked)
#   ranked_entrez.rds   — Named numeric vector ranked by Entrez ID (only if full_ranked)
#   contrast_info.rds   — Metadata list for downstream scripts
# ══════════════════════════════════════════════════════════════════════════════

# ── Load libraries ───────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(tidyr)
})

# ── Resolve project root and source utility modules ──────────────────────────
script_dir   <- dirname(sub("--file=", "",
                             grep("--file=", commandArgs(trailingOnly = FALSE),
                                  value = TRUE)))
project_root <- dirname(script_dir)

source(file.path(project_root, "lib", "config_utils.R"))
source(file.path(project_root, "lib", "data_utils.R"))

# ── Parse CLI arguments ──────────────────────────────────────────────────────
args <- parse_pipeline_args()

cat("════════════════════════════════════════════════════════════════════\n")
cat("  SBMA RNAseq Pipeline — Step 01: Setup & Data Preparation\n")
cat("════════════════════════════════════════════════════════════════════\n")
cat(sprintf("  Config   : %s\n", args$config))
cat(sprintf("  Dataset  : %s\n", args$dataset))
cat(sprintf("  Contrast : %s\n", args$contrast))
cat("────────────────────────────────────────────────────────────────────\n\n")

# ── 1. Load and validate configuration ───────────────────────────────────────
cat("[1/9] Loading configuration...\n")
cfg <- load_config(args$config)
cat(sprintf("       Species        : %s\n", cfg$species))
cat(sprintf("       OrgDb          : %s\n", cfg$orgdb))
cat(sprintf("       LFC threshold  : %s\n", cfg$thresholds$lfc))
cat(sprintf("       padj threshold : %s\n", cfg$thresholds$padj))

# ── 2. Extract contrast configuration ────────────────────────────────────────
cat("\n[2/9] Extracting contrast configuration...\n")
contrast_cfg <- get_contrast(cfg, args$dataset, args$contrast)
cat(sprintf("       Label              : %s\n", contrast_cfg$label))
cat(sprintf("       DEG file           : %s\n", contrast_cfg$deg_file))
cat(sprintf("       Tissue type        : %s\n", contrast_cfg$tissue_type))
cat(sprintf("       Full ranked        : %s\n", contrast_cfg$full_ranked))
cat(sprintf("       Androgen treatment : %s\n", contrast_cfg$androgen_treatment))

# ── 3. Create output directory structure ─────────────────────────────────────
cat("\n[3/9] Creating output directories...\n")
out_dir <- get_output_dir(cfg, args$dataset, args$contrast)
data_dir <- file.path(out_dir, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("       Output base : %s\n", out_dir))

# ── 4. Load species-specific OrgDb ───────────────────────────────────────────
cat("\n[4/9] Loading OrgDb annotation package...\n")
orgdb_pkg <- cfg$orgdb

if (!requireNamespace(orgdb_pkg, quietly = TRUE)) {
  stop(sprintf("OrgDb package '%s' is not installed. Install with:\n  BiocManager::install('%s')",
               orgdb_pkg, orgdb_pkg),
       call. = FALSE)
}

suppressPackageStartupMessages(library(orgdb_pkg, character.only = TRUE))
cat(sprintf("       Loaded: %s\n", orgdb_pkg))

# ── 5. Load DEG table ────────────────────────────────────────────────────────
cat("\n[5/9] Loading DEG table...\n")
deg_df <- load_deg_table(contrast_cfg)
cat(sprintf("       Rows loaded           : %d\n", nrow(deg_df)))
cat(sprintf("       Unique gene symbols   : %d\n", n_distinct(deg_df$gene_symbol)))
cat(sprintf("       log2FC range          : [%.2f, %.2f]\n",
            min(deg_df$log2fc, na.rm = TRUE),
            max(deg_df$log2fc, na.rm = TRUE)))

# ── 5b. Data provenance checks ──────────────────────────────────────────────
cat("\n[5b/9] Data provenance checks...\n")

# Check if padj == pvalue (no FDR correction applied)
if (all(deg_df$padj == deg_df$pvalue, na.rm = TRUE)) {
  cat("       WARNING: padj values are identical to raw p-values.\n")
  cat("       No FDR correction was applied to this dataset.\n")
  cat("       ORA results using padj thresholds may have inflated false\n")
  cat("       discovery rates. GSEA (rank-based) is not affected.\n")
}

# Detect suspicious log2fc patterns
has_negative_lfc <- any(deg_df$log2fc < 0, na.rm = TRUE)
max_abs_lfc      <- max(abs(deg_df$log2fc), na.rm = TRUE)

if (!has_negative_lfc && max_abs_lfc > 5) {
  cat("       WARNING: log2fc values are all non-negative with large range.\n")
  cat("       This may indicate raw fold changes (not log2-transformed).\n")
}

# Cross-check gene count vs full_ranked setting
n_input_genes <- nrow(deg_df)
if (n_input_genes < 1000 && isTRUE(contrast_cfg$full_ranked)) {
  cat(sprintf("       WARNING: Only %d genes but full_ranked=TRUE.\n", n_input_genes))
  cat("       A full ranked list typically contains >10,000 genes.\n")
}
if (n_input_genes > 5000 && !isTRUE(contrast_cfg$full_ranked)) {
  cat(sprintf("       NOTE: %d genes with full_ranked=FALSE.\n", n_input_genes))
  cat("       Consider setting full_ranked=TRUE to enable GSEA.\n")
}

# ── 6. Classify DEGs ─────────────────────────────────────────────────────────
cat("\n[6/9] Classifying DEGs...\n")
lfc_thresh  <- cfg$thresholds$lfc
padj_thresh <- cfg$thresholds$padj

deg_df <- classify_deg(deg_df, lfc_thresh, padj_thresh)

n_up   <- sum(deg_df$diffexpressed == "Upregulated")
n_down <- sum(deg_df$diffexpressed == "Downregulated")
n_ns   <- sum(deg_df$diffexpressed == "NS")

cat(sprintf("       Thresholds: |log2FC| > %s, padj < %s\n", lfc_thresh, padj_thresh))
cat(sprintf("       Upregulated   : %d\n", n_up))
cat(sprintf("       Downregulated : %d\n", n_down))
cat(sprintf("       NS            : %d\n", n_ns))
cat(sprintf("       Total DEGs    : %d\n", n_up + n_down))

# ── 7. Split DEGs by direction ───────────────────────────────────────────────
cat("\n[7/9] Splitting DEGs by direction...\n")
deg_splits <- split_deg(deg_df)
cat(sprintf("       Up genes         : %d\n", length(deg_splits$up)))
cat(sprintf("       Down genes       : %d\n", length(deg_splits$down)))
cat(sprintf("       All significant  : %d\n", length(deg_splits$all_sig)))
cat(sprintf("       Background genes : %d\n", length(deg_splits$background)))

# ── 8. Build ranked lists (if full_ranked) ───────────────────────────────────
cat("\n[8/9] Building ranked lists...\n")
ranked_symbol <- NULL
ranked_entrez <- NULL

if (isTRUE(contrast_cfg$full_ranked)) {
  cat("       full_ranked = TRUE -> creating ranked lists\n")

  ranked_symbol <- make_ranked_list(deg_df)
  cat(sprintf("       Ranked (symbol) : %d genes\n", length(ranked_symbol)))

  ranked_entrez <- make_ranked_entrez(deg_df, orgdb_pkg)
  cat(sprintf("       Ranked (entrez) : %d genes\n", length(ranked_entrez)))
} else {
  cat("       full_ranked = FALSE -> skipping ranked list generation\n")
}

# ── 9. Save outputs ──────────────────────────────────────────────────────────
cat("\n[9/9] Saving outputs...\n")

saveRDS(deg_df, file.path(data_dir, "deg_all.rds"))
cat(sprintf("       Saved: %s\n", file.path(data_dir, "deg_all.rds")))

saveRDS(deg_splits, file.path(data_dir, "deg_splits.rds"))
cat(sprintf("       Saved: %s\n", file.path(data_dir, "deg_splits.rds")))

if (isTRUE(contrast_cfg$full_ranked)) {
  saveRDS(ranked_symbol, file.path(data_dir, "ranked_symbol.rds"))
  cat(sprintf("       Saved: %s\n", file.path(data_dir, "ranked_symbol.rds")))

  saveRDS(ranked_entrez, file.path(data_dir, "ranked_entrez.rds"))
  cat(sprintf("       Saved: %s\n", file.path(data_dir, "ranked_entrez.rds")))
}

# -- Contrast info metadata --
contrast_info <- list(
  dataset_id         = args$dataset,
  contrast_id        = args$contrast,
  label              = contrast_cfg$label,
  tissue_type        = contrast_cfg$tissue_type,
  full_ranked        = contrast_cfg$full_ranked,
  androgen_treatment = contrast_cfg$androgen_treatment,
  thresholds         = list(lfc = lfc_thresh, padj = padj_thresh),
  n_up               = n_up,
  n_down             = n_down,
  n_total            = n_up + n_down
)

saveRDS(contrast_info, file.path(data_dir, "contrast_info.rds"))
cat(sprintf("       Saved: %s\n", file.path(data_dir, "contrast_info.rds")))

# ── Summary ──────────────────────────────────────────────────────────────────
cat("\n════════════════════════════════════════════════════════════════════\n")
cat("  Step 01 complete.\n")
cat(sprintf("  Dataset  : %s\n", args$dataset))
cat(sprintf("  Contrast : %s (%s)\n", args$contrast, contrast_cfg$label))
cat(sprintf("  DEGs     : %d up, %d down, %d total\n", n_up, n_down, n_up + n_down))
cat(sprintf("  Output   : %s\n", data_dir))
cat("════════════════════════════════════════════════════════════════════\n")

#!/usr/bin/env Rscript
# ==============================================================================
# 07_mitochondrial_module.R
# Mitochondrial function deep-dive for SBMA RNAseq pipeline
#
# Performs comprehensive analysis of mitochondrial gene expression:
#   - Curated mito gene overlay on DEG results
#   - Per-category summary statistics
#   - OXPHOS complex-by-complex analysis with t-tests
#   - Violin and barplot visualizations
#   - Category heatmaps (full + focused)
#   - KEGG pathview for OXPHOS (optional)
#   - Custom ORA per complex
#   - Search global ORA results for mitochondrial terms
#
# Usage:
#   Rscript scripts/07_mitochondrial_module.R --config config/config.yaml \
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
source(file.path(project_root, "lib", "enrichment_utils.R"))
source(file.path(project_root, "lib", "plot_utils.R"))
source(file.path(project_root, "lib", "gene_lists.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(clusterProfiler)
})

# ---- Parse CLI arguments and load config -------------------------------------

args <- parse_pipeline_args()

cfg <- load_config(args$config)

contrast_cfg <- get_contrast(cfg, args$dataset, args$contrast)
out_dir      <- get_output_dir(cfg, args$dataset, args$contrast)

ds <- args$dataset
ct <- args$contrast

message("=== 07_mitochondrial_module.R ===")
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
message(sprintf("Loaded %d genes from deg_all.rds.", nrow(deg_all)))

# ---- Species parameters ------------------------------------------------------

sp_params   <- get_species_params(cfg)
orgdb       <- sp_params$orgdb
kegg_org    <- sp_params$kegg_organism

if (!requireNamespace(orgdb, quietly = TRUE)) {
  stop(sprintf("Annotation package '%s' is not installed. "
               , "Install with BiocManager::install('%s').",
               orgdb, orgdb), call. = FALSE)
}
library(orgdb, character.only = TRUE)

# ---- Thresholds --------------------------------------------------------------

lfc_thresh  <- cfg$thresholds$lfc
padj_thresh <- cfg$thresholds$padj

# ==============================================================================
# STEP 4-6: Mitochondrial Gene Analysis
# ==============================================================================

message("\n--- Mitochondrial Gene Analysis ---")

# Step 4: Get curated mitochondrial genes
mito_genes <- get_mitochondrial_genes()
message(sprintf("Curated mitochondrial gene list: %d genes across %d categories.",
                nrow(mito_genes), length(unique(mito_genes$category))))

# Step 5: Merge with DEG data
mito_deg <- merge_with_deg(
  gene_list   = mito_genes,
  deg_df      = deg_all,
  lfc_thresh  = lfc_thresh,
  padj_thresh = padj_thresh
)

n_measured <- sum(!is.na(mito_deg$log2fc))
n_sig      <- sum(mito_deg$sig, na.rm = TRUE)
message(sprintf("Measured: %d / %d | Significant: %d (Up: %d, Down: %d)",
                n_measured, nrow(mito_deg), n_sig,
                sum(mito_deg$direction == "Up", na.rm = TRUE),
                sum(mito_deg$direction == "Down", na.rm = TRUE)))

# Step 6: Save mitochondrial gene summary
mito_summary_path <- file.path(tables_dir, "mitochondrial_gene_summary.csv")
write.csv(mito_deg, mito_summary_path, row.names = FALSE)
message(sprintf("Saved: %s", mito_summary_path))

# ==============================================================================
# STEP 7-8: Per-Category Summary
# ==============================================================================

message("\n--- Per-Category Summary ---")

category_summary <- mito_deg %>%
  group_by(category) %>%
  summarise(
    n_genes    = n(),
    n_measured = sum(!is.na(log2fc)),
    n_sig      = sum(sig, na.rm = TRUE),
    n_up       = sum(direction == "Up", na.rm = TRUE),
    n_down     = sum(direction == "Down", na.rm = TRUE),
    mean_lfc   = mean(log2fc, na.rm = TRUE),
    median_lfc = median(log2fc, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  arrange(desc(abs(mean_lfc)))

cat_summary_path <- file.path(tables_dir, "mitochondrial_category_summary.csv")
write.csv(category_summary, cat_summary_path, row.names = FALSE)
message(sprintf("Saved: %s", cat_summary_path))

# ==============================================================================
# STEP 9-13: OXPHOS Complex-by-Complex Analysis
# ==============================================================================

message("\n--- OXPHOS Complex-by-Complex Analysis ---")

# Step 9: Filter to OXPHOS genes (complex is not NA)
oxphos_deg <- mito_deg %>%
  filter(!is.na(complex))

message(sprintf("OXPHOS genes: %d total, %d measured.",
                nrow(oxphos_deg), sum(!is.na(oxphos_deg$log2fc))))

# Step 10: Per-complex statistics with one-sample t-test
complex_order <- c("Complex I", "Complex II", "Complex III",
                   "Complex IV", "Complex V")

complex_stats <- oxphos_deg %>%
  filter(!is.na(log2fc)) %>%
  group_by(complex) %>%
  summarise(
    n_measured = n(),
    n_sig      = sum(sig, na.rm = TRUE),
    mean_lfc   = mean(log2fc, na.rm = TRUE),
    median_lfc = median(log2fc, na.rm = TRUE),
    sem        = sd(log2fc, na.rm = TRUE) / sqrt(n()),
    ttest_p    = tryCatch(
      t.test(log2fc, mu = 0)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    sig_star = case_when(
      ttest_p < 0.001 ~ "***",
      ttest_p < 0.01  ~ "**",
      ttest_p < 0.05  ~ "*",
      TRUE            ~ ""
    ),
    complex = factor(complex, levels = complex_order)
  ) %>%
  arrange(complex)

# Step 10b: Permutation testing for complex-level mean LFC shifts
# Validates whether observed mean LFC per complex is more extreme than
# expected by randomly sampling gene sets of the same size from the genome.

message("\n--- Permutation Testing for OXPHOS Complexes ---")

n_perm <- 10000
set.seed(cfg$enrichment$seed %||% 42)

# All measured gene LFCs (genome-wide background for permutation)
all_measured_lfc <- deg_all %>%
  dplyr::filter(!is.na(log2fc)) %>%
  dplyr::pull(log2fc)

if (length(all_measured_lfc) > 100) {

  perm_results <- complex_stats %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      perm_p = {
        obs_mean <- mean_lfc
        n_genes  <- n_measured
        # Sample n_genes from all measured LFCs, n_perm times
        perm_means <- replicate(n_perm, mean(sample(all_measured_lfc, n_genes, replace = FALSE)))
        # Two-sided empirical p-value
        (sum(abs(perm_means) >= abs(obs_mean)) + 1) / (n_perm + 1)
      }
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      perm_sig_star = dplyr::case_when(
        perm_p < 0.001 ~ "***",
        perm_p < 0.01  ~ "**",
        perm_p < 0.05  ~ "*",
        TRUE           ~ ""
      )
    )

  # Add permutation results to complex_stats
  complex_stats <- complex_stats %>%
    dplyr::left_join(
      perm_results %>% dplyr::select(complex, perm_p, perm_sig_star),
      by = "complex"
    )

  message(sprintf("Permutation test (%d iterations) complete.", n_perm))
  message(sprintf("  %-15s %10s %10s %10s %s",
                  "Complex", "t-test P", "Perm P", "MeanLFC", "Perm Sig"))
  message(paste(rep("-", 60), collapse = ""))
  for (i in seq_len(nrow(complex_stats))) {
    row <- complex_stats[i, ]
    message(sprintf("  %-15s %10.4f %10.4f %10.3f %s",
                    as.character(row$complex),
                    row$ttest_p, row$perm_p,
                    row$mean_lfc, row$perm_sig_star))
  }

} else {
  message("Fewer than 100 measured genes in genome — skipping permutation test.")
  complex_stats$perm_p <- NA_real_
  complex_stats$perm_sig_star <- ""
}

# Step 11: Save complex stats (now includes permutation p-values)
complex_stats_path <- file.path(tables_dir, "oxphos_complex_stats.csv")
write.csv(complex_stats, complex_stats_path, row.names = FALSE)
message(sprintf("Saved: %s", complex_stats_path))

# Step 12: Complex summary barplot using make_complex_summary()
# Prepare data in the format expected by make_complex_summary():
# columns: complex, condition, mean_lfc, sem, pval
complex_plot_df <- complex_stats %>%
  mutate(
    condition = contrast_cfg$label,
    pval      = ttest_p
  ) %>%
  dplyr::select(complex, condition, mean_lfc, sem, pval)

colors <- default_plot_colors()

p_complex <- make_complex_summary(complex_plot_df, colors = colors)
p_complex <- p_complex +
  labs(title = sprintf("OXPHOS Complex Mean LFC — %s", contrast_cfg$label))

# Step 13: Save complex summary barplot
complex_plot_base <- file.path(plots_dir, "oxphos_complex_summary")
save_plot(p_complex, complex_plot_base, width = 8, height = 6,
          dpi = cfg$plots$dpi, formats = cfg$plots$format)

# ==============================================================================
# STEP 14-16: Violin Plot by Complex
# ==============================================================================

message("\n--- OXPHOS Violin Plot ---")

oxphos_measured <- oxphos_deg %>%
  filter(!is.na(log2fc)) %>%
  mutate(complex = factor(complex, levels = complex_order))

if (nrow(oxphos_measured) > 0) {

  p_violin <- ggplot(oxphos_measured,
                     aes(x = complex, y = log2fc, fill = complex)) +
    geom_violin(alpha = 0.6, scale = "width", trim = FALSE) +
    geom_jitter(aes(colour = direction),
                width = 0.15, size = 1.5, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40",
               linewidth = 0.5) +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    scale_colour_manual(
      values = c("Up" = colors$up, "Down" = colors$down, "NS" = colors$ns),
      name   = "Direction"
    ) +
    labs(
      title = sprintf("OXPHOS Gene Expression by Complex — %s",
                       contrast_cfg$label),
      x     = NULL,
      y     = expression(log[2] ~ "Fold Change")
    ) +
    get_theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  violin_plot_base <- file.path(plots_dir, "oxphos_violin_by_complex")
  save_plot(p_violin, violin_plot_base, width = 10, height = 7,
            dpi = cfg$plots$dpi, formats = cfg$plots$format)

} else {
  message("No measured OXPHOS genes — skipping violin plot.")
}

# ==============================================================================
# STEP 17-19: Category Heatmaps
# ==============================================================================

message("\n--- Category Heatmaps ---")

# Step 17: Full landscape heatmap — all mito genes with LFC
mito_for_heatmap <- mito_deg %>%
  filter(!is.na(log2fc)) %>%
  arrange(category, desc(abs(log2fc)))

if (nrow(mito_for_heatmap) > 0) {

  # Use make_heatmap with single LFC column
  hm_full <- make_heatmap(
    gene_df      = mito_for_heatmap,
    value_cols   = "log2fc",
    category_col = "category",
    title        = sprintf("Mitochondrial Genes LFC — %s", contrast_cfg$label),
    colors       = colors
  )

  # Determine height based on gene count
  n_genes_hm   <- nrow(mito_for_heatmap)
  hm_height    <- max(10, n_genes_hm * 0.12 + 2)
  hm_full_base <- file.path(plots_dir, "mito_full_heatmap")
  save_plot(hm_full, hm_full_base, width = 8, height = hm_height,
            dpi = cfg$plots$dpi, formats = cfg$plots$format)

} else {
  message("No measured mitochondrial genes — skipping full heatmap.")
}

# Step 18: Focused heatmaps for key categories
focused_categories <- c("fission", "fusion", "mitophagy", "ros_defense",
                         "tca_cycle")

for (focus_cat in focused_categories) {

  focus_df <- mito_deg %>%
    filter(category == focus_cat, !is.na(log2fc)) %>%
    arrange(desc(abs(log2fc)))

  if (nrow(focus_df) < 2) {
    message(sprintf("  %s: fewer than 2 measured genes — skipping.", focus_cat))
    next
  }

  hm_focus <- make_heatmap(
    gene_df      = focus_df,
    value_cols   = "log2fc",
    category_col = "category",
    title        = sprintf("%s Genes LFC — %s",
                           gsub("_", " ", tools::toTitleCase(focus_cat)),
                           contrast_cfg$label),
    colors       = colors
  )

  focus_height <- max(5, nrow(focus_df) * 0.3 + 2)
  focus_base   <- file.path(plots_dir,
                             sprintf("heatmap_%s", focus_cat))
  save_plot(hm_focus, focus_base, width = 7, height = focus_height,
            dpi = cfg$plots$dpi, formats = cfg$plots$format)
}

# ==============================================================================
# STEP 20-21: KEGG Pathview for OXPHOS (optional)
# ==============================================================================

message("\n--- KEGG Pathview (optional) ---")

tryCatch({
  if (!requireNamespace("pathview", quietly = TRUE)) {
    message("pathview package not available — skipping KEGG pathway visualization.")
  } else {

    # Get OXPHOS genes with measured LFC
    oxphos_for_pathview <- oxphos_deg %>%
      filter(!is.na(log2fc))

    if (nrow(oxphos_for_pathview) > 0) {

      # Convert gene symbols to Entrez IDs
      entrez_map <- symbols_to_entrez(oxphos_for_pathview$gene_symbol, orgdb)

      if (length(entrez_map) > 0) {

        # Build named LFC vector keyed by Entrez ID
        lfc_entrez <- oxphos_for_pathview %>%
          filter(gene_symbol %in% names(entrez_map)) %>%
          mutate(entrez_id = entrez_map[gene_symbol]) %>%
          # Deduplicate Entrez IDs, keep highest absolute LFC
          arrange(desc(abs(log2fc))) %>%
          distinct(entrez_id, .keep_all = TRUE)

        lfc_vec <- setNames(lfc_entrez$log2fc, lfc_entrez$entrez_id)

        # Determine pathway ID based on species
        pathway_id <- ifelse(kegg_org == "mmu", "mmu00190", "hsa00190")

        # Run pathview from the plots directory so output lands there
        old_wd <- getwd()
        setwd(plots_dir)

        pathview::pathview(
          gene.data  = lfc_vec,
          pathway.id = pathway_id,
          species    = kegg_org,
          out.suffix = "oxphos_lfc",
          kegg.native = TRUE,
          limit      = list(gene = max(abs(lfc_vec), na.rm = TRUE), cpd = 1)
        )

        setwd(old_wd)
        message("Pathview OXPHOS map saved to plots/mitochondrial/")

      } else {
        message("No OXPHOS genes mapped to Entrez IDs — skipping pathview.")
      }
    } else {
      message("No measured OXPHOS genes — skipping pathview.")
    }
  }
}, error = function(e) {
  message(sprintf("Pathview failed (non-critical): %s", conditionMessage(e)))
  # Restore working directory if pathview failed mid-execution
  tryCatch(setwd(project_root), error = function(e2) NULL)
})

# Ensure we are back in project root
setwd(project_root)

# ==============================================================================
# STEP 22-24: Custom ORA Per Complex
# ==============================================================================

message("\n--- Custom ORA Per Complex ---")

# Step 22: Build TERM2GENE from complex -> gene_symbol
oxphos_t2g <- mito_genes %>%
  filter(!is.na(complex)) %>%
  dplyr::select(complex, gene_symbol) %>%
  dplyr::rename(term = complex)

# Significant DEG symbols (from full DEG set, not just mito)
sig_genes <- deg_all %>%
  filter(abs(log2fc) > lfc_thresh, padj < padj_thresh) %>%
  pull(gene_symbol) %>%
  unique()

background_genes <- unique(deg_all$gene_symbol)

message(sprintf("ORA input: %d significant DEGs, %d background, %d complex gene sets.",
                length(sig_genes), length(background_genes), nrow(oxphos_t2g)))

# Step 23: Run enricher() on sig DEGs vs complex gene sets
ora_complex <- tryCatch({
  if (length(sig_genes) < 3) {
    message("Fewer than 3 significant DEGs — skipping complex ORA.")
    NULL
  } else {
    enricher(
      gene          = sig_genes,
      universe      = background_genes,
      TERM2GENE     = oxphos_t2g,
      pAdjustMethod = "BH",
      pvalueCutoff  = 1,      # Keep all results for reporting
      qvalueCutoff  = 1,
      minGSSize     = 3,
      maxGSSize     = 500
    )
  }
}, error = function(e) {
  message(sprintf("Complex ORA failed: %s", conditionMessage(e)))
  NULL
})

# Step 24: Save ORA results
if (!is.null(ora_complex)) {
  ora_df <- as.data.frame(ora_complex)
  if (nrow(ora_df) > 0) {
    ora_path <- file.path(tables_dir, "oxphos_ora_by_complex.csv")
    write.csv(ora_df, ora_path, row.names = FALSE)
    message(sprintf("Saved: %s (%d terms)", ora_path, nrow(ora_df)))
  } else {
    message("No enriched OXPHOS complexes found.")
  }
} else {
  message("Complex ORA returned NULL — no results to save.")
}

# ==============================================================================
# STEP 25-27: Search Global Results for Mito Terms
# ==============================================================================

message("\n--- Searching Global Enrichment Results for Mito Terms ---")

mito_pattern <- "mitochond|oxidative.phosph|electron.transport|respiratory.chain|OXPHOS|TCA|citrate|succinate"

# Step 25: Load ORA results from tables/ora_*.rds
ora_files <- list.files(tables_dir, pattern = "^ora_.*\\.rds$", full.names = TRUE)

mito_hits_list <- list()

if (length(ora_files) > 0) {

  for (ora_file in ora_files) {

    ora_name <- tools::file_path_sans_ext(basename(ora_file))

    tryCatch({
      ora_res <- readRDS(ora_file)
      ora_df  <- as.data.frame(ora_res)

      if (nrow(ora_df) > 0 && "Description" %in% colnames(ora_df)) {
        # Step 26: Search for mito-related terms
        hits <- ora_df %>%
          filter(grepl(mito_pattern, Description, ignore.case = TRUE)) %>%
          mutate(source = ora_name)

        if (nrow(hits) > 0) {
          mito_hits_list[[ora_name]] <- hits
          message(sprintf("  %s: %d mito-related terms found.", ora_name,
                          nrow(hits)))
        }
      }
    }, error = function(e) {
      message(sprintf("  %s: failed to read — %s", ora_name,
                      conditionMessage(e)))
    })
  }

} else {
  message("No ORA .rds files found in tables directory.")
}

# Step 27: Combine and save
if (length(mito_hits_list) > 0) {
  mito_global_hits <- bind_rows(mito_hits_list)
  global_hits_path <- file.path(tables_dir,
                                 "mitochondrial_global_enrichment_hits.csv")
  write.csv(mito_global_hits, global_hits_path, row.names = FALSE)
  message(sprintf("Saved: %s (%d hits across %d sources)",
                  global_hits_path, nrow(mito_global_hits),
                  length(mito_hits_list)))
} else {
  message("No mitochondrial-related terms found in global ORA results.")
}

# ==============================================================================
# STEP 28: Comprehensive Summary
# ==============================================================================

message("\n========== Mitochondrial Module Summary ==========")
message(sprintf("%-30s %s", "Dataset:", ds))
message(sprintf("%-30s %s (%s)", "Contrast:", ct, contrast_cfg$label))
message(paste(rep("-", 55), collapse = ""))

message(sprintf("%-30s %d", "Total curated mito genes:", nrow(mito_genes)))
message(sprintf("%-30s %d", "Measured in dataset:", n_measured))
message(sprintf("%-30s %d", "Significant DEGs:", n_sig))
message(sprintf("%-30s %d", "  Upregulated:",
                sum(mito_deg$direction == "Up", na.rm = TRUE)))
message(sprintf("%-30s %d", "  Downregulated:",
                sum(mito_deg$direction == "Down", na.rm = TRUE)))
message(paste(rep("-", 55), collapse = ""))

message("\nOXPHOS Complex Statistics:")
has_perm <- "perm_p" %in% colnames(complex_stats) && !all(is.na(complex_stats$perm_p))

if (has_perm) {
  message(sprintf("  %-15s %5s %5s %8s %8s %10s %10s %s",
                  "Complex", "Meas", "Sig", "MeanLFC", "SEM", "t-test P", "Perm P", ""))
  message(paste(rep("-", 80), collapse = ""))
  for (i in seq_len(nrow(complex_stats))) {
    row <- complex_stats[i, ]
    message(sprintf("  %-15s %5d %5d %8.3f %8.3f %10.4f %10.4f %s",
                    as.character(row$complex),
                    row$n_measured, row$n_sig,
                    row$mean_lfc, row$sem,
                    row$ttest_p, row$perm_p,
                    row$perm_sig_star))
  }
} else {
  message(sprintf("  %-15s %5s %5s %8s %8s %10s %s",
                  "Complex", "Meas", "Sig", "MeanLFC", "SEM", "t-test P", ""))
  message(paste(rep("-", 65), collapse = ""))
  for (i in seq_len(nrow(complex_stats))) {
    row <- complex_stats[i, ]
    message(sprintf("  %-15s %5d %5d %8.3f %8.3f %10.4f %s",
                    as.character(row$complex),
                    row$n_measured, row$n_sig,
                    row$mean_lfc, row$sem,
                    row$ttest_p,
                    row$sig_star))
  }
}

message(paste(rep("-", 55), collapse = ""))
message("\nPer-Category Summary:")
message(sprintf("  %-20s %5s %5s %5s %8s",
                "Category", "Meas", "Sig", "Up", "MeanLFC"))
message(paste(rep("-", 50), collapse = ""))

for (i in seq_len(nrow(category_summary))) {
  row <- category_summary[i, ]
  message(sprintf("  %-20s %5d %5d %5d %8.3f",
                  row$category,
                  row$n_measured, row$n_sig, row$n_up,
                  row$mean_lfc))
}

message(paste(rep("-", 55), collapse = ""))

message("\nOutput files:")
message(sprintf("  Tables: %s", tables_dir))
message(sprintf("  Plots:  %s", plots_dir))
message("=================================================\n")

message("07_mitochondrial_module.R completed successfully.")

#!/usr/bin/env Rscript
# ==============================================================================
# 06_purinergic_module.R
# SBMA RNAseq Generalized Pipeline -- Purinergic Signaling Deep-Dive Analysis
#
# Comprehensive analysis of purinergic signaling gene expression changes:
#   - Curated gene-level summaries and per-category statistics
#   - Heatmap and forest plot visualizations
#   - Custom ORA and GSEA on purinergic subcategories
#   - Fisher's exact test for enrichment of purinergic genes in DEGs
#   - P2X7-inflammasome axis focused analysis
#   - Calcium signaling gene identification
#   - Cross-reference with global ORA/GSEA enrichment results
#
# Usage:
#   Rscript scripts/06_purinergic_module.R \
#     --config config/config.yaml \
#     --dataset SBMA_iPSC_Q51 \
#     --contrast disease_vs_control_vehicle
#
# Key Snakemake output:
#   results/<ds>/<ct>/tables/purinergic_gene_summary.csv
# ==============================================================================

# -- Source pipeline libraries -------------------------------------------------
script_dir   <- dirname(sub("--file=", "",
                             grep("--file=", commandArgs(trailingOnly = FALSE),
                                  value = TRUE)))
project_root <- dirname(script_dir)

source(file.path(project_root, "lib", "config_utils.R"))
source(file.path(project_root, "lib", "data_utils.R"))
source(file.path(project_root, "lib", "enrichment_utils.R"))
source(file.path(project_root, "lib", "plot_utils.R"))
source(file.path(project_root, "lib", "gene_lists.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(clusterProfiler)
})

# ==============================================================================
# 1. Parse CLI arguments, load config, resolve paths
# ==============================================================================

args     <- parse_pipeline_args()
cfg      <- load_config(args$config)
contrast <- get_contrast(cfg, args$dataset, args$contrast)
out_dir  <- get_output_dir(cfg, args$dataset, args$contrast)

# Resolve visual settings
colors    <- if (!is.null(cfg$colors)) cfg$colors else default_plot_colors()
base_size <- if (!is.null(cfg$plots$base_size)) cfg$plots$base_size else 14
pub_theme <- get_theme_publication(base_size)
plot_dpi  <- if (!is.null(cfg$plots$dpi)) cfg$plots$dpi else 300
plot_fmts <- if (!is.null(cfg$plots$format)) cfg$plots$format else c("pdf", "png")

# Thresholds
lfc_thresh  <- cfg$thresholds$lfc
padj_thresh <- cfg$thresholds$padj

# Output subdirectories
tables_dir <- file.path(out_dir, "tables")
plots_dir  <- file.path(out_dir, "plots", "purinergic")
data_dir   <- file.path(out_dir, "data")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir,  recursive = TRUE, showWarnings = FALSE)

cat("================================================================\n")
cat(sprintf("  06 Purinergic Signaling Deep-Dive -- %s\n", contrast$label))
cat(sprintf("  Dataset : %s | Contrast: %s\n", args$dataset, args$contrast))
cat("================================================================\n\n")

# ==============================================================================
# 2. Load prepped DEG data
# ==============================================================================

deg_rds <- file.path(data_dir, "deg_all.rds")
if (!file.exists(deg_rds)) {
  stop(sprintf("Prepped DEG file not found: %s\nRun script 01 first.", deg_rds),
       call. = FALSE)
}

deg <- readRDS(deg_rds)
cat(sprintf("Loaded %d genes from %s\n\n", nrow(deg), deg_rds))

# ==============================================================================
# 3. Purinergic Gene Analysis
# ==============================================================================

cat("[1/11] Extracting curated purinergic gene list...\n")
purinergic_genes <- get_purinergic_genes()
cat(sprintf("       Curated genes: %d across %d categories\n",
            nrow(purinergic_genes),
            n_distinct(purinergic_genes$category)))

# -- Merge with DEG data -------------------------------------------------------
cat("[2/11] Merging purinergic genes with DEG results...\n")
purinergic_deg <- merge_with_deg(purinergic_genes, deg,
                                  lfc_thresh  = lfc_thresh,
                                  padj_thresh = padj_thresh)

# Flag genes that were measured (present in DEG data)
purinergic_deg <- purinergic_deg %>%
  mutate(measured = !is.na(log2fc))

n_measured <- sum(purinergic_deg$measured)
n_sig      <- sum(purinergic_deg$sig, na.rm = TRUE)
cat(sprintf("       Measured in data: %d / %d\n", n_measured, nrow(purinergic_deg)))
cat(sprintf("       Significant DEGs: %d\n", n_sig))

# -- Save primary output (Snakemake target) ------------------------------------
gene_summary_csv <- file.path(tables_dir, "purinergic_gene_summary.csv")
write_csv(purinergic_deg, gene_summary_csv)
cat(sprintf("       Saved: %s\n\n", gene_summary_csv))

# ==============================================================================
# 4. Per-Category Summary Statistics
# ==============================================================================

cat("[3/11] Computing per-category summary statistics...\n")

category_summary <- purinergic_deg %>%
  group_by(category) %>%
  summarise(
    n_genes    = n(),
    n_measured = sum(measured),
    n_sig      = sum(sig, na.rm = TRUE),
    n_up       = sum(direction == "Up", na.rm = TRUE),
    n_down     = sum(direction == "Down", na.rm = TRUE),
    mean_lfc   = mean(log2fc, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  arrange(desc(n_sig))

category_csv <- file.path(tables_dir, "purinergic_category_summary.csv")
write_csv(category_summary, category_csv)
cat(sprintf("       Saved: %s\n", category_csv))
cat("\n--- Per-Category Summary ---\n")
print(as.data.frame(category_summary), row.names = FALSE)
cat("\n")

# ==============================================================================
# 5. Heatmap of Purinergic Genes
# ==============================================================================

cat("[4/11] Creating purinergic gene heatmap...\n")

# Prepare data for heatmap: only measured genes, replace NA LFC with 0
heatmap_df <- purinergic_deg %>%
  filter(measured) %>%
  mutate(log2fc = replace_na(log2fc, 0)) %>%
  arrange(category, gene_symbol)

if (nrow(heatmap_df) >= 2) {
  p_heatmap <- make_heatmap(
    gene_df      = heatmap_df,
    value_cols   = "log2fc",
    category_col = "category",
    title        = sprintf("Purinergic Genes -- %s", contrast$label),
    colors       = colors
  )

  # Compute dynamic height based on gene count
  hm_height <- max(6, nrow(heatmap_df) * 0.25 + 2)

  save_plot(
    plot    = p_heatmap,
    path    = file.path(plots_dir, "purinergic_heatmap"),
    width   = 8,
    height  = hm_height,
    dpi     = plot_dpi,
    formats = plot_fmts
  )
} else {
  cat("       Fewer than 2 measured genes -- skipping heatmap.\n")
}

# ==============================================================================
# 6. Forest Plot
# ==============================================================================

cat("[5/11] Creating purinergic forest plot...\n")

if (nrow(heatmap_df) >= 2) {
  p_forest <- make_forest_plot(
    gene_df      = heatmap_df,
    value_cols   = "log2fc",
    category_col = "category",
    sig_col      = "direction",
    colors       = colors
  )

  forest_height <- max(6, nrow(heatmap_df) * 0.3 + 2)

  save_plot(
    plot    = p_forest,
    path    = file.path(plots_dir, "purinergic_forest_plot"),
    width   = 10,
    height  = forest_height,
    dpi     = plot_dpi,
    formats = plot_fmts
  )
} else {
  cat("       Fewer than 2 measured genes -- skipping forest plot.\n")
}

# ==============================================================================
# 7. Custom ORA on Purinergic Subcategories
# ==============================================================================

cat("[6/11] Running custom ORA on purinergic subcategories...\n")

# Build TERM2GENE: category -> gene_symbol
purinergic_t2g <- purinergic_genes %>%
  dplyr::select(term = category, gene = gene_symbol)

# Significant DEGs from the full dataset
if (!"diffexpressed" %in% colnames(deg)) {
  deg <- classify_deg(deg, lfc_thresh = lfc_thresh, padj_thresh = padj_thresh)
}
deg_sig_symbols <- deg %>%
  filter(diffexpressed != "NS") %>%
  pull(gene_symbol)

ora_result <- NULL
if (length(deg_sig_symbols) >= 3) {
  ora_result <- tryCatch({
    enricher(
      gene         = deg_sig_symbols,
      universe     = deg$gene_symbol,
      TERM2GENE    = purinergic_t2g,
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      minGSSize    = 2,
      maxGSSize    = 500
    )
  }, error = function(e) {
    message(sprintf("  Custom ORA failed: %s", conditionMessage(e)))
    NULL
  })
}

if (!is.null(ora_result) && nrow(as.data.frame(ora_result)) > 0) {
  ora_df <- as.data.frame(ora_result)
  ora_csv <- file.path(tables_dir, "purinergic_ora.csv")
  write_csv(ora_df, ora_csv)
  cat(sprintf("       Saved: %s (%d terms)\n", ora_csv, nrow(ora_df)))

  # Dotplot
  tryCatch({
    n_show <- min(20, nrow(ora_df))
    p_ora <- dotplot(ora_result, showCategory = n_show) +
      labs(title = sprintf("Purinergic Subcategory ORA -- %s", contrast$label)) +
      pub_theme

    save_plot(
      plot    = p_ora,
      path    = file.path(plots_dir, "purinergic_ora_dotplot"),
      width   = 10,
      height  = max(5, n_show * 0.4 + 2),
      dpi     = plot_dpi,
      formats = plot_fmts
    )
  }, error = function(e) {
    message(sprintf("  ORA dotplot failed: %s", conditionMessage(e)))
  })
} else {
  cat("       No enriched purinergic subcategories found.\n")
}

# ==============================================================================
# 8. Custom GSEA on Purinergic Subcategories (if full_ranked)
# ==============================================================================

cat("[7/11] Running custom GSEA on purinergic subcategories...\n")

ranked_rds <- file.path(data_dir, "ranked_symbol.rds")

if (file.exists(ranked_rds)) {
  ranked_symbol <- readRDS(ranked_rds)
  cat(sprintf("       Loaded ranked list: %d genes\n", length(ranked_symbol)))

  gsea_result <- tryCatch({
    set.seed(if (!is.null(cfg$enrichment$seed)) cfg$enrichment$seed else 42)
    GSEA(
      geneList      = ranked_symbol,
      TERM2GENE     = purinergic_t2g,
      pAdjustMethod = "BH",
      pvalueCutoff  = 1,
      minGSSize     = 3,
      maxGSSize     = 500,
      nPermSimple   = if (!is.null(cfg$enrichment$n_perm)) cfg$enrichment$n_perm else 10000
    )
  }, error = function(e) {
    message(sprintf("  Custom GSEA failed: %s", conditionMessage(e)))
    NULL
  })

  if (!is.null(gsea_result) && nrow(as.data.frame(gsea_result)) > 0) {
    gsea_df <- as.data.frame(gsea_result)
    gsea_csv <- file.path(tables_dir, "purinergic_gsea.csv")
    write_csv(gsea_df, gsea_csv)
    cat(sprintf("       Saved: %s (%d terms)\n", gsea_csv, nrow(gsea_df)))

    # Running score plots for each significant term
    tryCatch({
      for (i in seq_len(min(nrow(gsea_df), 10))) {
        term_id <- gsea_df$ID[i]
        safe_name <- gsub("[^A-Za-z0-9_]", "_", term_id)

        p_rs <- enrichplot::gseaplot2(gsea_result, geneSetID = i,
                                       title = term_id)

        save_plot(
          plot    = p_rs,
          path    = file.path(plots_dir, paste0("gsea_running_", safe_name)),
          width   = 8,
          height  = 6,
          dpi     = plot_dpi,
          formats = plot_fmts
        )
      }
    }, error = function(e) {
      message(sprintf("  GSEA running score plots failed: %s", conditionMessage(e)))
    })
  } else {
    cat("       No GSEA results for purinergic subcategories.\n")
  }
} else {
  cat("       ranked_symbol.rds not found -- skipping GSEA.\n")
}

# ==============================================================================
# 9. Fisher's Exact Test
# ==============================================================================

cat("[8/11] Running Fisher's exact test for purinergic gene enrichment...\n")

# Build 2x2 contingency table
all_genes       <- deg$gene_symbol
purinergic_set  <- purinergic_genes$gene_symbol
de_genes        <- deg_sig_symbols

is_purinergic <- all_genes %in% purinergic_set
is_de         <- all_genes %in% de_genes

# Contingency counts
a <- sum(is_purinergic & is_de)       # Purinergic & DE
b <- sum(is_purinergic & !is_de)      # Purinergic & not DE
c <- sum(!is_purinergic & is_de)      # Not purinergic & DE
d <- sum(!is_purinergic & !is_de)     # Not purinergic & not DE

contingency <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
                       dimnames = list(c("Purinergic", "Non_purinergic"),
                                       c("DE", "Not_DE")))

fisher_res <- fisher.test(contingency)

fisher_df <- tibble(
  group             = c("Purinergic_DE", "Purinergic_notDE",
                         "NonPurinergic_DE", "NonPurinergic_notDE"),
  count             = c(a, b, c, d),
  odds_ratio        = fisher_res$estimate,
  p_value           = fisher_res$p.value,
  conf_int_low      = fisher_res$conf.int[1],
  conf_int_high     = fisher_res$conf.int[2],
  alternative       = fisher_res$alternative
)

fisher_csv <- file.path(tables_dir, "purinergic_fisher_test.csv")
write_csv(fisher_df, fisher_csv)
cat(sprintf("       Contingency: Purinergic&DE=%d, Purinergic&notDE=%d, Other&DE=%d, Other&notDE=%d\n",
            a, b, c, d))
cat(sprintf("       Odds ratio: %.3f | P-value: %.4g\n",
            fisher_res$estimate, fisher_res$p.value))
cat(sprintf("       Saved: %s\n\n", fisher_csv))

# ==============================================================================
# 10. P2X7-Inflammasome Axis
# ==============================================================================

cat("[9/11] Analyzing P2X7-inflammasome axis...\n")

inflammasome_symbols <- get_inflammasome_genes()
inflammasome_deg <- deg %>%
  filter(gene_symbol %in% inflammasome_symbols) %>%
  mutate(
    sig = (abs(log2fc) > lfc_thresh & padj < padj_thresh),
    direction = case_when(
      sig & log2fc > 0 ~ "Up",
      sig & log2fc < 0 ~ "Down",
      TRUE             ~ "NS"
    ),
    sig_label = case_when(
      padj < 0.001 ~ "***",
      padj < 0.01  ~ "**",
      padj < 0.05  ~ "*",
      TRUE         ~ ""
    )
  )

if (nrow(inflammasome_deg) > 0) {
  # Save table
  inflammasome_csv <- file.path(tables_dir, "purinergic_inflammasome_axis.csv")
  write_csv(inflammasome_deg, inflammasome_csv)
  cat(sprintf("       Found %d / %d inflammasome genes in data\n",
              nrow(inflammasome_deg), length(inflammasome_symbols)))
  cat(sprintf("       Saved: %s\n", inflammasome_csv))

  # Barplot
  bar_df <- inflammasome_deg %>%
    mutate(
      bar_color = case_when(
        direction == "Up"   ~ "pos",
        direction == "Down" ~ "neg",
        TRUE                ~ "ns"
      ),
      gene_symbol = factor(gene_symbol, levels = gene_symbol[order(log2fc)])
    )

  p_inflamm <- ggplot(bar_df, aes(x = gene_symbol, y = log2fc, fill = bar_color)) +
    geom_col(width = 0.7, alpha = 0.9) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.5) +
    geom_hline(yintercept = c(-lfc_thresh, lfc_thresh),
               linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    geom_text(aes(label = sig_label,
                  y = log2fc + sign(log2fc) * max(abs(log2fc), na.rm = TRUE) * 0.05),
              vjust = ifelse(bar_df$log2fc >= 0, 0, 1),
              size = 4, show.legend = FALSE) +
    scale_fill_manual(
      values = c("pos" = colors$up, "neg" = colors$down, "ns" = colors$ns),
      guide  = "none"
    ) +
    labs(
      title = sprintf("P2X7-Inflammasome Axis -- %s", contrast$label),
      x     = NULL,
      y     = expression(log[2]~"Fold Change")
    ) +
    pub_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_plot(
    plot    = p_inflamm,
    path    = file.path(plots_dir, "inflammasome_axis_barplot"),
    width   = max(6, nrow(bar_df) * 0.6 + 2),
    height  = 6,
    dpi     = plot_dpi,
    formats = plot_fmts
  )
} else {
  cat("       No inflammasome genes detected in DEG data.\n")
}

# ==============================================================================
# 11. Calcium Signaling Genes
# ==============================================================================

cat("[10/11] Identifying calcium signaling genes...\n")

calcium_patterns <- get_calcium_patterns()

# Search DEG data for genes matching any calcium pattern
calcium_hits <- deg %>%
  filter(
    Reduce(`|`, lapply(calcium_patterns, function(pat) grepl(pat, gene_symbol)))
  ) %>%
  mutate(
    calcium_category = {
      cats <- sapply(gene_symbol, function(gs) {
        matched <- names(calcium_patterns)[
          sapply(calcium_patterns, function(pat) grepl(pat, gs))
        ]
        if (length(matched) == 0) NA_character_ else paste(matched, collapse = ";")
      })
      cats
    },
    sig = (abs(log2fc) > lfc_thresh & padj < padj_thresh),
    direction = case_when(
      sig & log2fc > 0 ~ "Up",
      sig & log2fc < 0 ~ "Down",
      TRUE             ~ "NS"
    )
  )

calcium_csv <- file.path(tables_dir, "purinergic_calcium_genes.csv")
write_csv(calcium_hits, calcium_csv)
cat(sprintf("       Found %d calcium-related genes in DEG data\n", nrow(calcium_hits)))
cat(sprintf("       Significant: %d\n", sum(calcium_hits$sig, na.rm = TRUE)))
cat(sprintf("       Saved: %s\n\n", calcium_csv))

# ==============================================================================
# 12. Search Global ORA/GSEA for Purinergic Terms
# ==============================================================================

cat("[11/11] Searching global enrichment results for purinergic terms...\n")

purinergic_pattern <- "purin|ATP|adenosine|P2X|P2Y|nucleotide"

# Scan for ORA/GSEA result RDS files
ora_rds_files <- list.files(tables_dir, pattern = "^ora_.*\\.rds$",
                             full.names = TRUE)
gsea_rds_files <- list.files(tables_dir, pattern = "^gsea_.*\\.rds$",
                              full.names = TRUE)
all_enrich_files <- c(ora_rds_files, gsea_rds_files)

global_hits <- tibble()

if (length(all_enrich_files) > 0) {
  for (rds_file in all_enrich_files) {
    tryCatch({
      res_obj <- readRDS(rds_file)
      if (is.null(res_obj)) next

      res_df <- as.data.frame(res_obj)
      if (nrow(res_df) == 0) next

      # Search Description column for purinergic-related terms
      desc_col <- if ("Description" %in% colnames(res_df)) "Description"
                  else if ("ID" %in% colnames(res_df)) "ID"
                  else NULL

      if (is.null(desc_col)) next

      matched <- res_df %>%
        filter(grepl(purinergic_pattern, .data[[desc_col]], ignore.case = TRUE)) %>%
        mutate(source_file = basename(rds_file))

      if (nrow(matched) > 0) {
        global_hits <- bind_rows(global_hits, matched)
      }
    }, error = function(e) {
      message(sprintf("  Could not process %s: %s",
                      basename(rds_file), conditionMessage(e)))
    })
  }
}

if (nrow(global_hits) > 0) {
  global_csv <- file.path(tables_dir, "purinergic_global_enrichment_hits.csv")
  write_csv(global_hits, global_csv)
  cat(sprintf("       Found %d purinergic-related terms across %d enrichment files\n",
              nrow(global_hits), n_distinct(global_hits$source_file)))
  cat(sprintf("       Saved: %s\n", global_csv))
} else {
  # Write empty file so Snakemake does not fail on missing output
  global_csv <- file.path(tables_dir, "purinergic_global_enrichment_hits.csv")
  write_csv(tibble(note = "No purinergic terms found in global enrichment results"),
            global_csv)
  cat("       No purinergic-related terms found in global enrichment results.\n")
  cat(sprintf("       Saved (empty): %s\n", global_csv))
}

# ==============================================================================
# Comprehensive Summary
# ==============================================================================

cat("\n================================================================\n")
cat("  Purinergic Module Summary\n")
cat("================================================================\n")
cat(sprintf("  Contrast           : %s (%s)\n", args$contrast, contrast$label))
cat(sprintf("  Dataset            : %s\n", args$dataset))
cat("----------------------------------------------------------------\n")
cat(sprintf("  Curated genes      : %d\n", nrow(purinergic_genes)))
cat(sprintf("  Measured in data   : %d\n", n_measured))
cat(sprintf("  Significant DEGs   : %d (|LFC| > %g, padj < %g)\n",
            n_sig, lfc_thresh, padj_thresh))
cat(sprintf("    - Upregulated    : %d\n",
            sum(purinergic_deg$direction == "Up", na.rm = TRUE)))
cat(sprintf("    - Downregulated  : %d\n",
            sum(purinergic_deg$direction == "Down", na.rm = TRUE)))
cat("----------------------------------------------------------------\n")
cat(sprintf("  Fisher's exact OR  : %.3f (P = %.4g)\n",
            fisher_res$estimate, fisher_res$p.value))
cat(sprintf("  Inflammasome genes : %d / %d measured\n",
            nrow(inflammasome_deg), length(inflammasome_symbols)))
cat(sprintf("  Calcium genes      : %d found (%d significant)\n",
            nrow(calcium_hits), sum(calcium_hits$sig, na.rm = TRUE)))
if (!is.null(ora_result) && nrow(as.data.frame(ora_result)) > 0) {
  cat(sprintf("  ORA subcategories  : %d enriched\n",
              nrow(as.data.frame(ora_result))))
}
if (exists("gsea_result") && !is.null(gsea_result) &&
    nrow(as.data.frame(gsea_result)) > 0) {
  cat(sprintf("  GSEA subcategories : %d terms\n",
              nrow(as.data.frame(gsea_result))))
}
cat(sprintf("  Global enrich hits : %d\n", nrow(global_hits)))
cat("----------------------------------------------------------------\n")
cat("  Output files:\n")
cat(sprintf("    %s\n", gene_summary_csv))
cat(sprintf("    %s\n", category_csv))
cat(sprintf("    %s\n", fisher_csv))
if (nrow(inflammasome_deg) > 0)
  cat(sprintf("    %s\n", inflammasome_csv))
cat(sprintf("    %s\n", calcium_csv))
cat(sprintf("    %s\n", global_csv))
cat("================================================================\n\n")

cat("Script 06 completed successfully.\n")

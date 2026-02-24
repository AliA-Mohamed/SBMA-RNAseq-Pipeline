#!/usr/bin/env Rscript
# ==============================================================================
# 09_integration_figures.R
# SBMA RNAseq Generalized Pipeline — Integration, Master Table, Summary Figures
#
# Aggregates all upstream analysis outputs for a single contrast into:
#   1. A master supplementary table (CSV + XLSX)
#   2. Multi-panel summary Figure 1 (volcano, hallmark, purinergic, OXPHOS)
#   3. Focused pathway Figure 2 (forest, heatmap, crosstalk, inflammasome)
#   4. Auto-generated plain-text analysis summary report
#
# Usage:
#   Rscript scripts/09_integration_figures.R \
#     --config config/config.yaml \
#     --dataset SBMA_iPSC_Q51 \
#     --contrast disease_vs_control_vehicle
# ==============================================================================

# ── Source pipeline libraries ────────────────────────────────────────────────
script_dir   <- dirname(sub("--file=", "",
                             grep("--file=", commandArgs(trailingOnly = FALSE),
                                  value = TRUE)))
project_root <- dirname(script_dir)
setwd(project_root)

source(file.path(project_root, "lib", "config_utils.R"))
source(file.path(project_root, "lib", "data_utils.R"))
source(file.path(project_root, "lib", "plot_utils.R"))
source(file.path(project_root, "lib", "gene_lists.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(cowplot)
  library(readr)
  library(purrr)
})

# ── Parse CLI arguments ─────────────────────────────────────────────────────
args <- parse_pipeline_args()

# ── Load config and contrast info ────────────────────────────────────────────
cfg      <- load_config(args$config)
contrast <- get_contrast(cfg, args$dataset, args$contrast)
out_dir  <- get_output_dir(cfg, args$dataset, args$contrast)

# ── Resolve plotting parameters ──────────────────────────────────────────────
colors    <- if (!is.null(cfg$colors)) cfg$colors else default_plot_colors()
base_size <- if (!is.null(cfg$plots$base_size)) cfg$plots$base_size else 14
pub_theme <- get_theme_publication(base_size)
plot_dpi  <- if (!is.null(cfg$plots$dpi)) cfg$plots$dpi else 300
plot_fmts <- if (!is.null(cfg$plots$format)) cfg$plots$format else c("pdf", "png")

lfc_thresh  <- cfg$thresholds$lfc
padj_thresh <- cfg$thresholds$padj

# ── Derive paths ─────────────────────────────────────────────────────────────
data_dir    <- file.path(out_dir, "data")
tables_dir  <- file.path(out_dir, "tables")
plots_dir   <- file.path(out_dir, "plots")
reports_dir <- file.path(out_dir, "reports")

dir.create(tables_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(reports_dir, recursive = TRUE, showWarnings = FALSE)

cat("================================================================\n")
cat(sprintf("  09 Integration & Summary Figures — %s\n", contrast$label))
cat(sprintf("  Dataset : %s | Contrast: %s\n", args$dataset, args$contrast))
cat("================================================================\n\n")

# ==============================================================================
# HELPER: safe loader for optional upstream outputs
# ==============================================================================

safe_load_rds <- function(path, label = basename(path)) {
  tryCatch({
    if (file.exists(path)) {
      obj <- readRDS(path)
      message(sprintf("  Loaded: %s", label))
      obj
    } else {
      message(sprintf("  Not found (skipping): %s", label))
      NULL
    }
  }, error = function(e) {
    message(sprintf("  Failed to load %s: %s", label, conditionMessage(e)))
    NULL
  })
}

safe_load_csv <- function(path, label = basename(path)) {
  tryCatch({
    if (file.exists(path)) {
      obj <- readr::read_csv(path, show_col_types = FALSE)
      message(sprintf("  Loaded: %s", label))
      obj
    } else {
      message(sprintf("  Not found (skipping): %s", label))
      NULL
    }
  }, error = function(e) {
    message(sprintf("  Failed to load %s: %s", label, conditionMessage(e)))
    NULL
  })
}

safe_load_enrichment <- function(path, label = basename(path)) {
  res <- safe_load_rds(path, label)
  if (is.null(res)) return(NULL)
  tryCatch({
    df <- as.data.frame(res)
    if (nrow(df) == 0) return(NULL)
    df
  }, error = function(e) {
    message(sprintf("  Could not coerce %s to data.frame: %s",
                    label, conditionMessage(e)))
    NULL
  })
}

# ==============================================================================
# 3. Load all previous outputs
# ==============================================================================

cat("Loading upstream outputs ...\n")

# DEG data (required)
deg_rds <- file.path(data_dir, "deg_all.rds")
if (!file.exists(deg_rds)) {
  stop(sprintf("Required DEG file not found: %s\nRun upstream scripts first.",
               deg_rds), call. = FALSE)
}
deg <- readRDS(deg_rds)
cat(sprintf("  Loaded %d genes from deg_all.rds\n", nrow(deg)))

# Ensure classification
if (!"diffexpressed" %in% colnames(deg)) {
  deg <- classify_deg(deg, lfc_thresh = lfc_thresh, padj_thresh = padj_thresh)
}

# Purinergic and mitochondrial gene summaries
puri_summary <- safe_load_csv(
  file.path(tables_dir, "purinergic_gene_summary.csv"),
  "purinergic_gene_summary.csv"
)

mito_summary <- safe_load_csv(
  file.path(tables_dir, "mitochondrial_gene_summary.csv"),
  "mitochondrial_gene_summary.csv"
)

# ORA results (Hallmarks)
ora_hallmark_up <- safe_load_enrichment(
  file.path(tables_dir, "ora_hallmarks_up.rds"),
  "ora_hallmarks_up.rds"
)

ora_hallmark_down <- safe_load_enrichment(
  file.path(tables_dir, "ora_hallmarks_down.rds"),
  "ora_hallmarks_down.rds"
)

# ORA results (GO BP)
ora_gobp_up <- safe_load_enrichment(
  file.path(tables_dir, "ora_go_bp_up.rds"),
  "ora_go_bp_up.rds"
)

ora_gobp_down <- safe_load_enrichment(
  file.path(tables_dir, "ora_go_bp_down.rds"),
  "ora_go_bp_down.rds"
)

# ORA results (KEGG)
ora_kegg_up <- safe_load_enrichment(
  file.path(tables_dir, "ora_kegg_up.rds"),
  "ora_kegg_up.rds"
)

ora_kegg_down <- safe_load_enrichment(
  file.path(tables_dir, "ora_kegg_down.rds"),
  "ora_kegg_down.rds"
)

# Crosstalk data
crosstalk_csv <- safe_load_csv(
  file.path(tables_dir, "crosstalk_gene_summary.csv"),
  "crosstalk_gene_summary.csv"
)

cat("\n")

# ==============================================================================
# 4-7. Master Supplementary Table
# ==============================================================================

cat("Building master supplementary table ...\n")

# Start with the full DEG table
master <- deg %>%
  select(gene_symbol, log2fc, pvalue, padj, diffexpressed)

# --- Purinergic annotation ---
puri_genes <- get_purinergic_genes()
if (!is.null(puri_summary) && "category" %in% colnames(puri_summary)) {
  puri_annot <- puri_summary %>%
    select(gene_symbol, purinergic_category = category) %>%
    distinct(gene_symbol, .keep_all = TRUE) %>%
    mutate(is_purinergic = TRUE)
} else {
  puri_annot <- puri_genes %>%
    select(gene_symbol, purinergic_category = category) %>%
    distinct(gene_symbol, .keep_all = TRUE) %>%
    mutate(is_purinergic = TRUE)
}

master <- master %>%
  left_join(puri_annot, by = "gene_symbol") %>%
  mutate(is_purinergic = ifelse(is.na(is_purinergic), FALSE, is_purinergic))

# --- Mitochondrial annotation ---
mito_genes <- get_mitochondrial_genes()
if (!is.null(mito_summary) && all(c("category", "complex") %in% colnames(mito_summary))) {
  mito_annot <- mito_summary %>%
    select(gene_symbol, mito_category = category, mito_complex = complex) %>%
    distinct(gene_symbol, .keep_all = TRUE) %>%
    mutate(is_mitochondrial = TRUE)
} else {
  mito_annot <- mito_genes %>%
    select(gene_symbol, mito_category = category, mito_complex = complex) %>%
    distinct(gene_symbol, .keep_all = TRUE) %>%
    mutate(is_mitochondrial = TRUE)
}

master <- master %>%
  left_join(mito_annot, by = "gene_symbol") %>%
  mutate(is_mitochondrial = ifelse(is.na(is_mitochondrial), FALSE, is_mitochondrial))

# --- Crosstalk annotation ---
crosstalk_genes <- get_crosstalk_genes()
master <- master %>%
  mutate(is_crosstalk = gene_symbol %in% crosstalk_genes$gene_symbol)

# --- Top hallmark term assignment ---
master$top_hallmark <- NA_character_

if (!is.null(ora_hallmark_up) || !is.null(ora_hallmark_down)) {
  # Combine hallmark results from both directions
  hallmark_combined <- bind_rows(
    if (!is.null(ora_hallmark_up))   ora_hallmark_up   %>% mutate(direction = "UP"),
    if (!is.null(ora_hallmark_down)) ora_hallmark_down %>% mutate(direction = "DOWN")
  )

  if (nrow(hallmark_combined) > 0 && "geneID" %in% colnames(hallmark_combined)) {
    # Parse geneID column (slash-separated gene lists) and map each gene
    # to the most significant hallmark term it belongs to
    hallmark_gene_map <- hallmark_combined %>%
      arrange(p.adjust) %>%
      mutate(term_rank = row_number()) %>%
      select(Description, geneID, term_rank) %>%
      mutate(genes = strsplit(geneID, "/")) %>%
      unnest(genes) %>%
      rename(gene_symbol = genes) %>%
      arrange(term_rank) %>%
      distinct(gene_symbol, .keep_all = TRUE) %>%
      select(gene_symbol, top_hallmark = Description)

    # Clean hallmark names: remove "HALLMARK_" prefix and replace underscores
    hallmark_gene_map <- hallmark_gene_map %>%
      mutate(top_hallmark = gsub("^HALLMARK_", "", top_hallmark),
             top_hallmark = gsub("_", " ", top_hallmark))

    # Merge into master
    master <- master %>%
      select(-top_hallmark) %>%
      left_join(hallmark_gene_map, by = "gene_symbol")
  }
}

# Reorder columns for clarity
master <- master %>%
  select(gene_symbol, log2fc, pvalue, padj, diffexpressed,
         is_purinergic, purinergic_category,
         is_mitochondrial, mito_category, mito_complex,
         is_crosstalk, top_hallmark)

# Save CSV
master_csv <- file.path(tables_dir, "master_supplementary_table.csv")
readr::write_csv(master, master_csv)
cat(sprintf("  Saved: %s (%d rows)\n", master_csv, nrow(master)))

# Save XLSX if writexl is available
tryCatch({
  if (requireNamespace("writexl", quietly = TRUE)) {
    master_xlsx <- file.path(tables_dir, "master_supplementary_table.xlsx")
    writexl::write_xlsx(master, master_xlsx)
    cat(sprintf("  Saved: %s\n", master_xlsx))
  } else {
    message("  writexl not available — XLSX export skipped.")
  }
}, error = function(e) {
  message(sprintf("  XLSX export failed: %s", conditionMessage(e)))
})

cat("\n")

# ==============================================================================
# 8-9. Multi-Panel Summary Figure (Figure 1)
# ==============================================================================

cat("Generating Figure 1: Multi-panel summary ...\n")

# --- Panel A: Volcano plot ---
p_volcano <- tryCatch({
  make_volcano(
    df          = deg,
    title       = "A. Volcano Plot",
    lfc_thresh  = lfc_thresh,
    padj_thresh = padj_thresh,
    colors      = colors,
    top_n_label = 10
  ) +
    theme(legend.position = "none",
          plot.title = element_text(size = 12, face = "bold"))
}, error = function(e) {
  message(sprintf("  Volcano plot failed: %s", conditionMessage(e)))
  ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Volcano plot\nnot available") +
    theme_void() + ggtitle("A. Volcano Plot")
})

# --- Panel B: Top 10 Hallmark pathways barplot ---
p_hallmark <- tryCatch({
  hallmark_frames <- list()

  if (!is.null(ora_hallmark_up) && nrow(ora_hallmark_up) > 0) {
    hallmark_frames[["up"]] <- ora_hallmark_up %>%
      arrange(p.adjust) %>%
      head(5) %>%
      mutate(Direction = "UP",
             neglog10p = -log10(pmax(p.adjust, 1e-300)))
  }

  if (!is.null(ora_hallmark_down) && nrow(ora_hallmark_down) > 0) {
    hallmark_frames[["down"]] <- ora_hallmark_down %>%
      arrange(p.adjust) %>%
      head(5) %>%
      mutate(Direction = "DOWN",
             neglog10p = -log10(pmax(p.adjust, 1e-300)))
  }

  if (length(hallmark_frames) == 0) stop("No hallmark results available.")

  hallmark_df <- bind_rows(hallmark_frames) %>%
    mutate(
      term_label = gsub("^HALLMARK_", "", Description),
      term_label = gsub("_", " ", term_label),
      term_label = ifelse(nchar(term_label) > 45,
                          paste0(substr(term_label, 1, 42), "..."),
                          term_label),
      Direction = factor(Direction, levels = c("UP", "DOWN"))
    ) %>%
    arrange(Direction, desc(neglog10p)) %>%
    mutate(term_label = factor(term_label, levels = rev(unique(term_label))))

  ggplot(hallmark_df, aes(x = neglog10p, y = term_label, fill = Direction)) +
    geom_col(alpha = 0.85, width = 0.7) +
    scale_fill_manual(values = c("UP" = colors$up, "DOWN" = colors$down)) +
    labs(title = "B. Top Hallmark Pathways",
         x = expression(-log[10]~"adj. P"),
         y = NULL, fill = NULL) +
    get_theme_publication(base_size = 10) +
    theme(axis.text.y = element_text(size = 8),
          plot.title = element_text(size = 12, face = "bold"),
          legend.position = "bottom")
}, error = function(e) {
  message(sprintf("  Hallmark barplot failed: %s", conditionMessage(e)))
  ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Hallmark results\nnot available") +
    theme_void() + ggtitle("B. Top Hallmark Pathways")
})

# --- Panel C: Purinergic category summary ---
p_purinergic_cat <- tryCatch({
  puri_merged <- merge_with_deg(puri_genes, deg,
                                lfc_thresh = lfc_thresh,
                                padj_thresh = padj_thresh)

  puri_cat_summary <- puri_merged %>%
    group_by(category) %>%
    summarise(
      n_total = n(),
      n_sig   = sum(sig, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(n_sig > 0 | n_total > 0) %>%
    arrange(desc(n_sig)) %>%
    mutate(category = factor(category, levels = rev(category)))

  ggplot(puri_cat_summary, aes(x = n_sig, y = category)) +
    geom_col(fill = colors$down, alpha = 0.85, width = 0.7) +
    geom_text(aes(label = sprintf("%d/%d", n_sig, n_total)),
              hjust = -0.1, size = 3) +
    labs(title = "C. Purinergic Signaling",
         x = "Significant Genes",
         y = NULL) +
    get_theme_publication(base_size = 10) +
    theme(plot.title = element_text(size = 12, face = "bold")) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.3)))
}, error = function(e) {
  message(sprintf("  Purinergic category plot failed: %s", conditionMessage(e)))
  ggplot() + annotate("text", x = 0.5, y = 0.5,
                       label = "Purinergic summary\nnot available") +
    theme_void() + ggtitle("C. Purinergic Signaling")
})

# --- Panel D: OXPHOS complex summary (mean LFC per complex with error bars) ---
p_oxphos <- tryCatch({
  mito_merged <- merge_with_deg(mito_genes, deg,
                                lfc_thresh = lfc_thresh,
                                padj_thresh = padj_thresh)

  # Filter to OXPHOS genes only (those with a complex assignment)
  oxphos_df <- mito_merged %>%
    filter(!is.na(complex), !is.na(log2fc))

  if (nrow(oxphos_df) == 0) stop("No OXPHOS genes with LFC data.")

  complex_stats <- oxphos_df %>%
    group_by(complex) %>%
    summarise(
      mean_lfc = mean(log2fc, na.rm = TRUE),
      sem      = sd(log2fc, na.rm = TRUE) / sqrt(n()),
      n_genes  = n(),
      n_sig    = sum(sig, na.rm = TRUE),
      pval     = tryCatch(t.test(log2fc, mu = 0)$p.value,
                           error = function(e) NA_real_),
      .groups  = "drop"
    ) %>%
    mutate(
      bar_color = ifelse(mean_lfc >= 0, "pos", "neg"),
      sig_star = case_when(
        pval < 0.001 ~ "***",
        pval < 0.01  ~ "**",
        pval < 0.05  ~ "*",
        TRUE         ~ ""
      )
    )

  ggplot(complex_stats, aes(x = complex, y = mean_lfc, fill = bar_color)) +
    geom_col(alpha = 0.9, width = 0.7) +
    geom_errorbar(aes(ymin = mean_lfc - sem, ymax = mean_lfc + sem),
                  width = 0.25, linewidth = 0.4) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.5) +
    geom_text(aes(label = sig_star,
                  y = mean_lfc + sign(mean_lfc) * (sem + max(abs(mean_lfc)) * 0.08)),
              size = 5, show.legend = FALSE) +
    scale_fill_manual(values = c("pos" = colors$up, "neg" = colors$down),
                      guide = "none") +
    labs(title = "D. OXPHOS Complex Summary",
         x = NULL,
         y = expression("Mean" ~ log[2] ~ "FC")) +
    get_theme_publication(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 12, face = "bold"))
}, error = function(e) {
  message(sprintf("  OXPHOS complex plot failed: %s", conditionMessage(e)))
  ggplot() + annotate("text", x = 0.5, y = 0.5,
                       label = "OXPHOS summary\nnot available") +
    theme_void() + ggtitle("D. OXPHOS Complex Summary")
})

# --- Assemble Figure 1 ---
fig1 <- cowplot::plot_grid(
  p_volcano, p_hallmark,
  p_purinergic_cat, p_oxphos,
  ncol   = 2,
  nrow   = 2,
  labels = NULL,
  align  = "hv",
  axis   = "tblr",
  rel_widths  = c(1, 1),
  rel_heights = c(1, 1)
)

fig1_title <- cowplot::ggdraw() +
  cowplot::draw_label(
    sprintf("Figure 1: %s — Multi-Panel Summary", contrast$label),
    fontface = "bold", size = 14, x = 0.5, hjust = 0.5
  )

fig1_final <- cowplot::plot_grid(fig1_title, fig1,
                                  ncol = 1, rel_heights = c(0.04, 1))

save_plot(
  plot    = fig1_final,
  path    = file.path(plots_dir, "summary_figure_1"),
  width   = 16,
  height  = 12,
  dpi     = plot_dpi,
  formats = plot_fmts
)

cat("  Figure 1 saved.\n\n")

# ==============================================================================
# 10-11. Focused Pathway Figure (Figure 2)
# ==============================================================================

cat("Generating Figure 2: Focused pathway panels ...\n")

# --- Panel A: Purinergic forest plot (top 20 significant genes) ---
p_puri_forest <- tryCatch({
  puri_merged <- merge_with_deg(puri_genes, deg,
                                lfc_thresh = lfc_thresh,
                                padj_thresh = padj_thresh)

  puri_sig <- puri_merged %>%
    filter(!is.na(log2fc)) %>%
    arrange(padj) %>%
    head(20) %>%
    mutate(
      gene_symbol = factor(gene_symbol, levels = rev(gene_symbol)),
      color_group = case_when(
        sig & log2fc > 0 ~ "Up",
        sig & log2fc < 0 ~ "Down",
        TRUE             ~ "NS"
      )
    )

  if (nrow(puri_sig) == 0) stop("No purinergic genes with expression data.")

  ggplot(puri_sig, aes(x = log2fc, y = gene_symbol, colour = color_group)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40",
               linewidth = 0.5) +
    geom_point(size = 2.5, alpha = 0.85) +
    scale_colour_manual(
      values = c("Up" = colors$up, "Down" = colors$down, "NS" = colors$ns),
      name = NULL
    ) +
    labs(title = "A. Purinergic Genes (Top 20)",
         x = expression(log[2]~"Fold Change"),
         y = NULL) +
    get_theme_publication(base_size = 10) +
    theme(axis.text.y = element_text(size = 8, face = "italic"),
          plot.title = element_text(size = 12, face = "bold"),
          legend.position = "bottom")
}, error = function(e) {
  message(sprintf("  Purinergic forest plot failed: %s", conditionMessage(e)))
  ggplot() + annotate("text", x = 0.5, y = 0.5,
                       label = "Purinergic forest plot\nnot available") +
    theme_void() + ggtitle("A. Purinergic Genes (Top 20)")
})

# --- Panel B: OXPHOS heatmap (top genes per complex) ---
p_oxphos_heat <- tryCatch({
  mito_merged <- merge_with_deg(mito_genes, deg,
                                lfc_thresh = lfc_thresh,
                                padj_thresh = padj_thresh)

  oxphos_heat <- mito_merged %>%
    filter(!is.na(complex), !is.na(log2fc)) %>%
    group_by(complex) %>%
    arrange(padj) %>%
    slice_head(n = 8) %>%
    ungroup() %>%
    arrange(complex, padj) %>%
    mutate(gene_symbol = factor(gene_symbol, levels = rev(unique(gene_symbol))))

  if (nrow(oxphos_heat) == 0) stop("No OXPHOS genes with expression data.")

  lfc_max <- max(abs(oxphos_heat$log2fc), na.rm = TRUE)

  ggplot(oxphos_heat, aes(x = complex, y = gene_symbol, fill = log2fc)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_gradient2(
      low = colors$heatmap_low, mid = colors$heatmap_mid,
      high = colors$heatmap_high, midpoint = 0,
      limits = c(-lfc_max, lfc_max),
      name = expression(log[2]~FC)
    ) +
    labs(title = "B. OXPHOS Complexes I-V",
         x = NULL, y = NULL) +
    get_theme_publication(base_size = 10) +
    theme(axis.text.y = element_text(size = 7, face = "italic"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          plot.title = element_text(size = 12, face = "bold"),
          legend.position = "right",
          legend.key.height = unit(0.8, "cm"))
}, error = function(e) {
  message(sprintf("  OXPHOS heatmap failed: %s", conditionMessage(e)))
  ggplot() + annotate("text", x = 0.5, y = 0.5,
                       label = "OXPHOS heatmap\nnot available") +
    theme_void() + ggtitle("B. OXPHOS Complexes I-V")
})

# --- Panel C: Crosstalk gene barplot (28 genes, colored by function_group) ---
p_crosstalk <- tryCatch({
  ct_genes <- get_crosstalk_genes()
  ct_merged <- ct_genes %>%
    left_join(
      deg %>% select(gene_symbol, log2fc, padj),
      by = "gene_symbol"
    ) %>%
    filter(!is.na(log2fc)) %>%
    mutate(
      gene_symbol = factor(gene_symbol,
                           levels = gene_symbol[order(function_group, log2fc)])
    )

  if (nrow(ct_merged) == 0) stop("No crosstalk genes with expression data.")

  # Generate a discrete palette for function groups
  n_groups <- length(unique(ct_merged$function_group))
  group_colors <- setNames(
    RColorBrewer::brewer.pal(max(3, min(n_groups, 8)), "Set2")[seq_len(n_groups)],
    sort(unique(ct_merged$function_group))
  )

  ggplot(ct_merged, aes(x = gene_symbol, y = log2fc, fill = function_group)) +
    geom_col(alpha = 0.85, width = 0.7) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.5) +
    scale_fill_manual(values = group_colors, name = "Function Group") +
    labs(title = "C. Crosstalk Genes",
         x = NULL,
         y = expression(log[2]~"FC")) +
    get_theme_publication(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                      size = 7, face = "italic"),
          plot.title = element_text(size = 12, face = "bold"),
          legend.position = "bottom",
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 8)) +
    guides(fill = guide_legend(nrow = 2))
}, error = function(e) {
  message(sprintf("  Crosstalk barplot failed: %s", conditionMessage(e)))
  ggplot() + annotate("text", x = 0.5, y = 0.5,
                       label = "Crosstalk barplot\nnot available") +
    theme_void() + ggtitle("C. Crosstalk Genes")
})

# --- Panel D: Inflammasome axis barplot ---
p_inflammasome <- tryCatch({
  inflamm_genes <- c("P2RX7", "NLRP3", "CASP1", "IL1B", "IL18")

  inflamm_df <- deg %>%
    filter(gene_symbol %in% inflamm_genes) %>%
    select(gene_symbol, log2fc, padj) %>%
    mutate(
      gene_symbol = factor(gene_symbol, levels = inflamm_genes),
      sig_label = case_when(
        padj < 0.001 ~ "***",
        padj < 0.01  ~ "**",
        padj < 0.05  ~ "*",
        TRUE         ~ ""
      ),
      bar_color = ifelse(log2fc >= 0, "pos", "neg")
    )

  if (nrow(inflamm_df) == 0) stop("No inflammasome genes found in DEG data.")

  ggplot(inflamm_df, aes(x = gene_symbol, y = log2fc, fill = bar_color)) +
    geom_col(alpha = 0.9, width = 0.7) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.5) +
    geom_text(aes(label = sig_label,
                  y = log2fc + sign(log2fc) * max(abs(log2fc), na.rm = TRUE) * 0.08),
              size = 5, show.legend = FALSE) +
    scale_fill_manual(values = c("pos" = colors$up, "neg" = colors$down),
                      guide = "none") +
    labs(title = "D. Inflammasome Axis",
         x = NULL,
         y = expression(log[2]~"Fold Change")) +
    get_theme_publication(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
          plot.title = element_text(size = 12, face = "bold"))
}, error = function(e) {
  message(sprintf("  Inflammasome barplot failed: %s", conditionMessage(e)))
  ggplot() + annotate("text", x = 0.5, y = 0.5,
                       label = "Inflammasome axis\nnot available") +
    theme_void() + ggtitle("D. Inflammasome Axis")
})

# --- Assemble Figure 2 ---
fig2 <- cowplot::plot_grid(
  p_puri_forest, p_oxphos_heat,
  p_crosstalk, p_inflammasome,
  ncol   = 2,
  nrow   = 2,
  labels = NULL,
  align  = "hv",
  axis   = "tblr",
  rel_widths  = c(1, 1),
  rel_heights = c(1, 1)
)

fig2_title <- cowplot::ggdraw() +
  cowplot::draw_label(
    sprintf("Figure 2: %s — Focused Pathway Analysis", contrast$label),
    fontface = "bold", size = 14, x = 0.5, hjust = 0.5
  )

fig2_final <- cowplot::plot_grid(fig2_title, fig2,
                                  ncol = 1, rel_heights = c(0.04, 1))

save_plot(
  plot    = fig2_final,
  path    = file.path(plots_dir, "summary_figure_2"),
  width   = 16,
  height  = 12,
  dpi     = plot_dpi,
  formats = plot_fmts
)

cat("  Figure 2 saved.\n\n")

# ==============================================================================
# 12-13. Analysis Summary Report
# ==============================================================================

cat("Generating analysis summary report ...\n")

# -- Compute summary statistics ---
n_total <- nrow(deg)
n_up    <- sum(deg$diffexpressed == "Upregulated",   na.rm = TRUE)
n_down  <- sum(deg$diffexpressed == "Downregulated", na.rm = TRUE)
n_ns    <- sum(deg$diffexpressed == "NS",            na.rm = TRUE)

# -- Helper: extract top N terms from enrichment result ---
.top_terms <- function(df, n = 5) {
  if (is.null(df) || nrow(df) == 0) return("  (no significant terms)")
  top <- df %>%
    arrange(p.adjust) %>%
    head(n)
  desc_col <- if ("Description" %in% colnames(top)) "Description" else "ID"
  paste0("  - ", top[[desc_col]], " (padj=",
         formatC(top$p.adjust, format = "e", digits = 2), ")",
         collapse = "\n")
}

# -- Purinergic summary ---
puri_for_report <- merge_with_deg(puri_genes, deg,
                                  lfc_thresh = lfc_thresh,
                                  padj_thresh = padj_thresh)
n_puri_total   <- sum(!is.na(puri_for_report$log2fc))
n_puri_sig     <- sum(puri_for_report$sig, na.rm = TRUE)
n_puri_defined <- nrow(puri_genes)

# -- Mitochondrial summary ---
mito_for_report <- merge_with_deg(mito_genes, deg,
                                  lfc_thresh = lfc_thresh,
                                  padj_thresh = padj_thresh)
n_mito_total   <- sum(!is.na(mito_for_report$log2fc))
n_mito_sig     <- sum(mito_for_report$sig, na.rm = TRUE)
n_mito_defined <- nrow(mito_genes)

# -- OXPHOS complex summary ---
oxphos_report <- mito_for_report %>%
  filter(!is.na(complex), !is.na(log2fc)) %>%
  group_by(complex) %>%
  summarise(
    n_genes  = n(),
    n_sig    = sum(sig, na.rm = TRUE),
    mean_lfc = mean(log2fc, na.rm = TRUE),
    pval     = tryCatch(t.test(log2fc, mu = 0)$p.value,
                         error = function(e) NA_real_),
    .groups  = "drop"
  )

oxphos_lines <- if (nrow(oxphos_report) > 0) {
  apply(oxphos_report, 1, function(row) {
    shift <- if (!is.na(as.numeric(row["pval"])) && as.numeric(row["pval"]) < 0.05) {
      "SIGNIFICANT shift"
    } else {
      "no significant shift"
    }
    sprintf("  - %s: %s sig of %s genes, mean LFC=%.3f (%s)",
            row["complex"], row["n_sig"], row["n_genes"],
            as.numeric(row["mean_lfc"]), shift)
  })
} else {
  "  (no OXPHOS data available)"
}

# -- Crosstalk summary ---
ct_for_report <- get_crosstalk_genes() %>%
  left_join(deg %>% select(gene_symbol, log2fc, padj), by = "gene_symbol") %>%
  mutate(sig = !is.na(padj) & padj < padj_thresh & abs(log2fc) > lfc_thresh)
n_ct_total   <- nrow(get_crosstalk_genes())
n_ct_measured <- sum(!is.na(ct_for_report$log2fc))
n_ct_sig     <- sum(ct_for_report$sig, na.rm = TRUE)

# -- Key findings (auto-derived) ---
key_findings <- character(0)

if (n_up > n_down * 2) {
  key_findings <- c(key_findings,
    sprintf("Strong upregulation bias: %d up vs %d down genes.", n_up, n_down))
} else if (n_down > n_up * 2) {
  key_findings <- c(key_findings,
    sprintf("Strong downregulation bias: %d down vs %d up genes.", n_down, n_up))
}

if (n_puri_sig > 0) {
  key_findings <- c(key_findings,
    sprintf("Purinergic signaling is significantly dysregulated: %d of %d measured genes are significant.",
            n_puri_sig, n_puri_total))
}

if (n_mito_sig > 0) {
  key_findings <- c(key_findings,
    sprintf("Mitochondrial gene dysregulation detected: %d of %d measured genes are significant.",
            n_mito_sig, n_mito_total))
}

sig_complexes <- oxphos_report %>%
  filter(!is.na(pval), pval < 0.05)
if (nrow(sig_complexes) > 0) {
  key_findings <- c(key_findings,
    sprintf("OXPHOS complexes with significant mean LFC shift: %s.",
            paste(sig_complexes$complex, collapse = ", ")))
}

if (n_ct_sig > 0) {
  key_findings <- c(key_findings,
    sprintf("Mitochondrial-purinergic crosstalk evidence: %d of %d crosstalk genes are significant.",
            n_ct_sig, n_ct_measured))
}

# Check inflammasome axis
inflamm_check <- deg %>% filter(gene_symbol %in% c("P2RX7", "NLRP3", "CASP1", "IL1B", "IL18"))
inflamm_sig   <- inflamm_check %>%
  filter(padj < padj_thresh, abs(log2fc) > lfc_thresh)
if (nrow(inflamm_sig) > 0) {
  key_findings <- c(key_findings,
    sprintf("Inflammasome axis engagement: %s are significantly altered.",
            paste(inflamm_sig$gene_symbol, collapse = ", ")))
}

if (length(key_findings) == 0) {
  key_findings <- "No outstanding patterns detected at current thresholds."
}

# -- Assemble report text ---
report_lines <- c(
  "==============================================================================",
  "ANALYSIS SUMMARY REPORT",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "==============================================================================",
  "",
  "--- Dataset and Contrast ---",
  sprintf("  Dataset:     %s", args$dataset),
  sprintf("  Contrast:    %s", contrast$label),
  sprintf("  Tissue:      %s", contrast$tissue_type),
  sprintf("  Species:     %s", cfg$species),
  "",
  "--- Thresholds ---",
  sprintf("  LFC threshold:  %g", lfc_thresh),
  sprintf("  padj threshold: %g", padj_thresh),
  "",
  "--- Differential Expression Summary ---",
  sprintf("  Total genes tested: %d", n_total),
  sprintf("  Upregulated:        %d", n_up),
  sprintf("  Downregulated:      %d", n_down),
  sprintf("  Not significant:    %d", n_ns),
  "",
  "--- Top Enriched Pathways ---",
  "",
  "  GO Biological Process (UP):",
  .top_terms(ora_gobp_up),
  "",
  "  GO Biological Process (DOWN):",
  .top_terms(ora_gobp_down),
  "",
  "  KEGG (UP):",
  .top_terms(ora_kegg_up),
  "",
  "  KEGG (DOWN):",
  .top_terms(ora_kegg_down),
  "",
  "  Hallmark (UP):",
  .top_terms(ora_hallmark_up),
  "",
  "  Hallmark (DOWN):",
  .top_terms(ora_hallmark_down),
  "",
  "--- Purinergic Signaling Summary ---",
  sprintf("  Defined purinergic genes:   %d", n_puri_defined),
  sprintf("  Measured in this dataset:   %d", n_puri_total),
  sprintf("  Significantly dysregulated: %d", n_puri_sig),
  "",
  "--- Mitochondrial Function Summary ---",
  sprintf("  Defined mitochondrial genes: %d", n_mito_defined),
  sprintf("  Measured in this dataset:    %d", n_mito_total),
  sprintf("  Significantly dysregulated:  %d", n_mito_sig),
  "",
  "--- OXPHOS Complex Summary ---",
  paste(oxphos_lines, collapse = "\n"),
  "",
  "--- Crosstalk Summary ---",
  sprintf("  Defined crosstalk genes:    %d", n_ct_total),
  sprintf("  Measured in this dataset:   %d", n_ct_measured),
  sprintf("  Significantly dysregulated: %d", n_ct_sig),
  "",
  "--- Key Findings ---",
  paste(paste0("  * ", key_findings), collapse = "\n"),
  "",
  "=============================================================================="
)

# Write report
report_path <- file.path(reports_dir, "analysis_summary.txt")
writeLines(report_lines, report_path)
cat(sprintf("  Saved: %s\n", report_path))

# ==============================================================================
# 14. Completion
# ==============================================================================

cat("\n")
cat("================================================================\n")
cat(sprintf("  09 Integration & Summary Figures — COMPLETE\n"))
cat(sprintf("  Dataset:  %s\n", args$dataset))
cat(sprintf("  Contrast: %s (%s)\n", args$contrast, contrast$label))
cat("================================================================\n")
cat("Outputs:\n")
cat(sprintf("  Master table: %s\n", master_csv))
cat(sprintf("  Figure 1:     %s/summary_figure_1.{pdf,png}\n", plots_dir))
cat(sprintf("  Figure 2:     %s/summary_figure_2.{pdf,png}\n", plots_dir))
cat(sprintf("  Report:       %s\n", report_path))
cat("================================================================\n")
cat("Script 09 completed successfully.\n")

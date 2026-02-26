#!/usr/bin/env Rscript
# ══════════════════════════════════════════════════════════════════════════════
# 10_cross_comparison.R — Cross-comparison analyses across contrasts/datasets
# ══════════════════════════════════════════════════════════════════════════════
#
# Usage:
#   Rscript scripts/10_cross_comparison.R --config config/config.yaml
#
# Description:
#   Performs within-dataset and across-dataset comparative analyses as
#   specified in the cross_comparisons section of config.yaml.
#
#   Within-dataset: UpSet plots, LFC scatter plots, compareCluster enrichment,
#   and condition-specific gene lists (if contrasts differ by androgen treatment).
#
#   Across-datasets: DEG overlap, pathway Jaccard similarity, purinergic and
#   mitochondrial signature comparisons, and direction concordance heatmaps.
#
# Inputs:
#   results/<dataset>/<contrast>/data/deg_all.rds  — classified DEG tables
#   results/<dataset>/<contrast>/tables/ora_*.rds  — ORA result objects
#   results/<dataset>/<contrast>/tables/purinergic_gene_summary.csv
#   results/<dataset>/<contrast>/tables/mitochondrial_gene_summary.csv
#
# Outputs (under results/cross_comparison/):
#   <dataset>/upset_plot.{pdf,png}
#   <dataset>/lfc_scatter.{pdf,png}
#   <dataset>/compareCluster_*.{pdf,png,csv}
#   <dataset>/condition_specific_*.{csv,pdf,png}
#   <comparison_name>/deg_overlap_upset.{pdf,png}
#   <comparison_name>/pathway_jaccard_heatmap.{pdf,png}
#   <comparison_name>/purinergic_comparison.{pdf,png}
#   <comparison_name>/mitochondrial_comparison.{pdf,png}
#   <comparison_name>/direction_concordance_heatmap.{pdf,png}
#   done.flag
# ══════════════════════════════════════════════════════════════════════════════

# ── Load libraries ───────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(UpSetR)
  library(clusterProfiler)
  library(ReactomePA)
  library(enrichplot)
})

# ── Resolve project root and source utility modules ──────────────────────────
script_dir   <- dirname(sub("--file=", "",
                             grep("--file=", commandArgs(trailingOnly = FALSE),
                                  value = TRUE)))
project_root <- dirname(script_dir)
setwd(project_root)

source(file.path(project_root, "lib", "config_utils.R"))
source(file.path(project_root, "lib", "data_utils.R"))
source(file.path(project_root, "lib", "plot_utils.R"))
source(file.path(project_root, "lib", "gene_lists.R"))

# ── Parse CLI arguments (only --config) ──────────────────────────────────────
option_list <- list(
  optparse::make_option(
    c("-c", "--config"),
    type    = "character",
    default = "config/config.yaml",
    help    = "Path to the pipeline config.yaml [default: %default]"
  )
)

parser <- optparse::OptionParser(
  usage       = "usage: %prog [options]",
  option_list = option_list,
  description = "SBMA RNAseq Pipeline — Step 10: Cross-Comparison Analyses"
)

args <- optparse::parse_args(parser)

cat("════════════════════════════════════════════════════════════════════\n")
cat("  SBMA RNAseq Pipeline — Step 10: Cross-Comparison Analyses\n")
cat("════════════════════════════════════════════════════════════════════\n")
cat(sprintf("  Config : %s\n", args$config))
cat("────────────────────────────────────────────────────────────────────\n\n")

# ── 1. Load and validate configuration ───────────────────────────────────────
cat("[1] Loading configuration...\n")
cfg <- load_config(args$config)

sp_params     <- get_species_params(cfg)
orgdb         <- sp_params$orgdb
kegg_org      <- sp_params$kegg_organism
reactome_org  <- sp_params$reactome_organism

# Load annotation package
if (!requireNamespace(orgdb, quietly = TRUE)) {
  stop(sprintf("Annotation package '%s' is not installed.", orgdb), call. = FALSE)
}
library(orgdb, character.only = TRUE)

species_full <- switch(tolower(cfg$species),
  human = "Homo sapiens",
  mouse = "Mus musculus",
  stop(sprintf("Unsupported species: '%s'", cfg$species), call. = FALSE)
)

lfc_thresh  <- cfg$thresholds$lfc
padj_thresh <- cfg$thresholds$padj
params      <- cfg$enrichment
plot_fmts   <- cfg$plots$format %||% c("pdf", "png")
colors      <- if (!is.null(cfg$colors)) cfg$colors else default_plot_colors()

cat(sprintf("       Species        : %s\n", cfg$species))
cat(sprintf("       LFC threshold  : %s\n", lfc_thresh))
cat(sprintf("       padj threshold : %s\n", padj_thresh))

# ── 2. Create base output directory ──────────────────────────────────────────
cat("\n[2] Creating output directories...\n")
cross_dir <- file.path("results", "cross_comparison")
dir.create(cross_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("       Base dir: %s\n", cross_dir))

# ══════════════════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

#' Load deg_all.rds for a given dataset/contrast pair
#' @return A tibble or NULL on failure
.load_deg_rds <- function(dataset_id, contrast_id) {
  path <- file.path("results", dataset_id, contrast_id, "data", "deg_all.rds")
  if (!file.exists(path)) {
    warning(sprintf("DEG file not found: %s", path), call. = FALSE)
    return(NULL)
  }
  tryCatch(readRDS(path), error = function(e) {
    warning(sprintf("Failed to read %s: %s", path, conditionMessage(e)),
            call. = FALSE)
    NULL
  })
}

#' Extract significant genes (up/down) from a classified DEG table
#' @return A list with up and down character vectors
.extract_sig <- function(deg_df, lfc_t = lfc_thresh, padj_t = padj_thresh) {
  if (!"diffexpressed" %in% colnames(deg_df)) {
    deg_df <- classify_deg(deg_df, lfc_thresh = lfc_t, padj_thresh = padj_t)
  }
  list(
    up   = deg_df$gene_symbol[deg_df$diffexpressed == "Upregulated"],
    down = deg_df$gene_symbol[deg_df$diffexpressed == "Downregulated"]
  )
}

#' Save a ggplot to multiple formats
.save_gg <- function(p, path_base, width, height, dpi = 300) {
  for (fmt in plot_fmts) {
    filepath <- paste0(path_base, ".", fmt)
    tryCatch({
      ggsave(filepath, plot = p, width = width, height = height,
             dpi = dpi, bg = "white")
      message(sprintf("  Saved: %s", filepath))
    }, error = function(e) {
      warning(sprintf("Failed to save %s: %s", filepath, conditionMessage(e)),
              call. = FALSE)
    })
  }
}

#' Save an UpSetR plot to PDF (UpSetR uses base graphics)
.save_upset <- function(upset_data, sets, path_base, width = 12, height = 8) {
  for (fmt in plot_fmts) {
    filepath <- paste0(path_base, ".", fmt)
    tryCatch({
      if (fmt == "pdf") {
        pdf(filepath, width = width, height = height)
      } else if (fmt == "png") {
        png(filepath, width = width, height = height, units = "in", res = 300)
      } else {
        next
      }
      print(UpSetR::upset(
        upset_data,
        sets            = sets,
        order.by        = "freq",
        keep.order      = TRUE,
        text.scale      = 1.3,
        point.size      = 3,
        line.size       = 1,
        mainbar.y.label = "Intersection Size",
        sets.x.label    = "Set Size"
      ))
      dev.off()
      message(sprintf("  Saved: %s", filepath))
    }, error = function(e) {
      tryCatch(dev.off(), error = function(x) NULL)
      warning(sprintf("Failed to save upset %s: %s", filepath,
                      conditionMessage(e)), call. = FALSE)
    })
  }
}


# ══════════════════════════════════════════════════════════════════════════════
# WITHIN-DATASET COMPARISONS
# ══════════════════════════════════════════════════════════════════════════════

within_comparisons <- cfg$cross_comparisons$within_dataset

if (!is.null(within_comparisons) && length(within_comparisons) > 0) {

  cat("\n════════════════════════════════════════════════════════════════════\n")
  cat("  WITHIN-DATASET COMPARISONS\n")
  cat("════════════════════════════════════════════════════════════════════\n")

  for (wc in within_comparisons) {

    dataset_id <- wc$dataset
    contrast_ids <- unlist(wc$contrasts)

    cat(sprintf("\n── Dataset: %s ──\n", dataset_id))
    cat(sprintf("   Contrasts: %s\n", paste(contrast_ids, collapse = ", ")))

    out_within <- file.path(cross_dir, dataset_id)
    dir.create(out_within, recursive = TRUE, showWarnings = FALSE)

    # ── Step 4: Load DEG tables for each contrast ────────────────────────────
    cat("   [4] Loading DEG tables...\n")

    deg_list <- list()
    sig_list <- list()
    all_loaded <- TRUE

    for (ct_id in contrast_ids) {
      deg_df <- .load_deg_rds(dataset_id, ct_id)
      if (is.null(deg_df)) {
        cat(sprintf("       WARNING: Could not load DEG data for %s/%s — skipping dataset.\n",
                    dataset_id, ct_id))
        all_loaded <- FALSE
        break
      }
      # Ensure classification
      if (!"diffexpressed" %in% colnames(deg_df)) {
        deg_df <- classify_deg(deg_df, lfc_thresh = lfc_thresh, padj_thresh = padj_thresh)
      }
      deg_list[[ct_id]] <- deg_df
      sig_list[[ct_id]] <- .extract_sig(deg_df)
      cat(sprintf("       %s: %d genes (%d UP, %d DOWN)\n",
                  ct_id, nrow(deg_df),
                  length(sig_list[[ct_id]]$up),
                  length(sig_list[[ct_id]]$down)))
    }

    if (!all_loaded) next

    # ── Step 5-6: UpSet Plot ─────────────────────────────────────────────────
    cat("   [5-6] Creating UpSet plot...\n")
    tryCatch({
      # Build binary membership matrix
      all_genes <- unique(unlist(lapply(deg_list, function(df) df$gene_symbol)))

      upset_sets <- list()
      set_names  <- character(0)

      for (ct_id in contrast_ids) {
        up_name   <- paste0("UP_", ct_id)
        down_name <- paste0("DOWN_", ct_id)
        upset_sets[[up_name]]   <- sig_list[[ct_id]]$up
        upset_sets[[down_name]] <- sig_list[[ct_id]]$down
        set_names <- c(set_names, up_name, down_name)
      }

      # Build fromList-compatible data
      upset_data <- UpSetR::fromList(upset_sets)

      if (nrow(upset_data) > 0 && ncol(upset_data) > 0) {
        .save_upset(upset_data, sets = set_names,
                    path_base = file.path(out_within, "upset_plot"))
      } else {
        cat("       No genes in any set — UpSet plot skipped.\n")
      }
    }, error = function(e) {
      warning(sprintf("UpSet plot failed for %s: %s", dataset_id,
                      conditionMessage(e)), call. = FALSE)
    })

    # ── Step 7: LFC Scatter Plot ─────────────────────────────────────────────
    if (length(contrast_ids) == 2) {
      cat("   [7] Creating LFC scatter plot...\n")
      tryCatch({
        ct1 <- contrast_ids[1]
        ct2 <- contrast_ids[2]

        # Get contrast labels from config
        ct1_label <- tryCatch(
          get_contrast(cfg, dataset_id, ct1)$label,
          error = function(e) ct1
        )
        ct2_label <- tryCatch(
          get_contrast(cfg, dataset_id, ct2)$label,
          error = function(e) ct2
        )

        df1 <- deg_list[[ct1]] %>%
          dplyr::select(gene_symbol, log2fc, padj, diffexpressed) %>%
          dplyr::rename(lfc_1 = log2fc, padj_1 = padj, de_1 = diffexpressed)

        df2 <- deg_list[[ct2]] %>%
          dplyr::select(gene_symbol, log2fc, padj, diffexpressed) %>%
          dplyr::rename(lfc_2 = log2fc, padj_2 = padj, de_2 = diffexpressed)

        scatter_df <- inner_join(df1, df2, by = "gene_symbol")

        # Classify quadrants
        sig1_up   <- scatter_df$de_1 == "Upregulated"
        sig1_down <- scatter_df$de_1 == "Downregulated"
        sig2_up   <- scatter_df$de_2 == "Upregulated"
        sig2_down <- scatter_df$de_2 == "Downregulated"
        sig1      <- sig1_up | sig1_down
        sig2      <- sig2_up | sig2_down

        scatter_df <- scatter_df %>%
          mutate(
            quadrant = case_when(
              (sig1_up & sig2_up) | (sig1_down & sig2_down) ~ "Concordant UP/DOWN",
              (sig1_up & sig2_down) | (sig1_down & sig2_up) ~ "Discordant",
              sig1 & !sig2 ~ paste0(ct1, " only"),
              !sig1 & sig2 ~ paste0(ct2, " only"),
              TRUE ~ "Neither"
            ),
            # Refine concordant
            quadrant = case_when(
              sig1_up & sig2_up     ~ "Concordant UP",
              sig1_down & sig2_down ~ "Concordant DOWN",
              TRUE                  ~ quadrant
            )
          )

        # Pearson correlation
        cor_val <- cor(scatter_df$lfc_1, scatter_df$lfc_2,
                       use = "complete.obs", method = "pearson")

        # Identify top 30 outlier genes (largest LFC difference)
        scatter_df <- scatter_df %>%
          mutate(lfc_diff = abs(lfc_1 - lfc_2))

        top_outliers <- scatter_df %>%
          filter(de_1 != "NS" | de_2 != "NS") %>%
          arrange(desc(lfc_diff)) %>%
          head(30)

        # Quadrant color mapping
        quad_colors <- c(
          "Concordant UP"   = colors$up,
          "Concordant DOWN" = colors$down,
          "Discordant"      = "#E69F00",
          "Neither"         = "grey75"
        )
        # Add dynamic labels for contrast-specific
        ct1_only_label <- paste0(ct1, " only")
        ct2_only_label <- paste0(ct2, " only")
        quad_colors[[ct1_only_label]] <- "#7570B3"
        quad_colors[[ct2_only_label]] <- "#66A61E"

        p_scatter <- ggplot(scatter_df, aes(x = lfc_1, y = lfc_2, colour = quadrant)) +
          geom_point(alpha = 0.5, size = 1) +
          geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                      colour = "grey40", linewidth = 0.5) +
          geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                      colour = "black", linewidth = 0.6, linetype = "solid",
                      inherit.aes = FALSE,
                      aes(x = lfc_1, y = lfc_2),
                      data = scatter_df) +
          geom_text_repel(
            data         = top_outliers,
            aes(label = gene_symbol),
            size         = 2.5,
            max.overlaps = 25,
            show.legend  = FALSE,
            fontface     = "italic",
            segment.color = "grey50",
            segment.size = 0.2
          ) +
          scale_colour_manual(values = quad_colors, name = "Category") +
          annotate("text",
                   x = min(scatter_df$lfc_1, na.rm = TRUE) * 0.9,
                   y = max(scatter_df$lfc_2, na.rm = TRUE) * 0.95,
                   label = sprintf("r = %.3f", cor_val),
                   fontface = "bold", size = 4, hjust = 0) +
          labs(
            title = sprintf("LFC Comparison: %s", dataset_id),
            x     = bquote(log[2]~FC ~ ": " ~ .(ct1_label)),
            y     = bquote(log[2]~FC ~ ": " ~ .(ct2_label))
          ) +
          get_theme_publication() +
          theme(legend.position = "right")

        .save_gg(p_scatter, file.path(out_within, "lfc_scatter"),
                 width = 10, height = 9)

        # Save quadrant counts
        quad_summary <- scatter_df %>%
          count(quadrant, name = "n_genes") %>%
          arrange(desc(n_genes))
        write.csv(quad_summary,
                  file.path(out_within, "lfc_scatter_quadrant_counts.csv"),
                  row.names = FALSE)
        cat(sprintf("       Pearson r = %.3f, %d shared genes\n",
                    cor_val, nrow(scatter_df)))

      }, error = function(e) {
        warning(sprintf("LFC scatter failed for %s: %s", dataset_id,
                        conditionMessage(e)), call. = FALSE)
      })
    } else {
      cat("   [7] LFC scatter skipped (requires exactly 2 contrasts).\n")
    }

    # ── Step 8: compareCluster enrichment ────────────────────────────────────
    cat("   [8] Running compareCluster enrichment...\n")

    tryCatch({
      # Build named gene clusters: UP_ct1, DOWN_ct1, UP_ct2, DOWN_ct2
      gene_clusters <- list()
      for (ct_id in contrast_ids) {
        up_genes   <- sig_list[[ct_id]]$up
        down_genes <- sig_list[[ct_id]]$down
        if (length(up_genes) >= 5)   gene_clusters[[paste0("UP_", ct_id)]]   <- up_genes
        if (length(down_genes) >= 5) gene_clusters[[paste0("DOWN_", ct_id)]] <- down_genes
      }

      if (length(gene_clusters) >= 2) {

        # --- GO_BP ---
        tryCatch({
          cat("       compareCluster: GO_BP...\n")
          cc_gobp <- compareCluster(
            geneClusters  = gene_clusters,
            fun           = "enrichGO",
            OrgDb         = orgdb,
            ont           = "BP",
            keyType       = "SYMBOL",
            pAdjustMethod = "BH",
            pvalueCutoff  = params$pvalue_cutoff,
            qvalueCutoff  = params$qvalue_cutoff,
            minGSSize     = params$min_gs_size,
            maxGSSize     = params$max_gs_size,
            readable      = FALSE
          )
          if (!is.null(cc_gobp) && nrow(as.data.frame(cc_gobp)) > 0) {
            p_cc <- dotplot(cc_gobp, showCategory = 10) +
              labs(title = sprintf("compareCluster GO:BP — %s", dataset_id)) +
              theme_minimal(base_size = 11)
            .save_gg(p_cc, file.path(out_within, "compareCluster_GO_BP"),
                     width = 14, height = 10)
            write.csv(as.data.frame(cc_gobp),
                      file.path(out_within, "compareCluster_GO_BP.csv"),
                      row.names = FALSE)
          } else {
            cat("       GO_BP: no significant terms.\n")
          }
        }, error = function(e) {
          warning(sprintf("compareCluster GO_BP: %s", conditionMessage(e)),
                  call. = FALSE)
        })

        # --- KEGG ---
        tryCatch({
          cat("       compareCluster: KEGG...\n")
          # Convert gene clusters to Entrez IDs
          entrez_clusters <- lapply(gene_clusters, function(genes) {
            mapped <- suppressWarnings(
              clusterProfiler::bitr(genes, fromType = "SYMBOL",
                                    toType = "ENTREZID", OrgDb = orgdb)
            )
            if (nrow(mapped) > 0) unique(mapped$ENTREZID) else character(0)
          })
          entrez_clusters <- entrez_clusters[sapply(entrez_clusters, length) >= 5]

          if (length(entrez_clusters) >= 2) {
            cc_kegg <- compareCluster(
              geneClusters  = entrez_clusters,
              fun           = "enrichKEGG",
              organism      = kegg_org,
              pAdjustMethod = "BH",
              pvalueCutoff  = params$pvalue_cutoff,
              qvalueCutoff  = params$qvalue_cutoff,
              minGSSize     = params$min_gs_size,
              maxGSSize     = params$max_gs_size
            )
            if (!is.null(cc_kegg) && nrow(as.data.frame(cc_kegg)) > 0) {
              p_cc <- dotplot(cc_kegg, showCategory = 10) +
                labs(title = sprintf("compareCluster KEGG — %s", dataset_id)) +
                theme_minimal(base_size = 11)
              .save_gg(p_cc, file.path(out_within, "compareCluster_KEGG"),
                       width = 14, height = 10)
              write.csv(as.data.frame(cc_kegg),
                        file.path(out_within, "compareCluster_KEGG.csv"),
                        row.names = FALSE)
            } else {
              cat("       KEGG: no significant terms.\n")
            }
          }
        }, error = function(e) {
          warning(sprintf("compareCluster KEGG: %s", conditionMessage(e)),
                  call. = FALSE)
        })

        # --- Reactome ---
        tryCatch({
          cat("       compareCluster: Reactome...\n")
          # Reuse entrez_clusters from KEGG block if available
          if (!exists("entrez_clusters") || length(entrez_clusters) < 2) {
            entrez_clusters <- lapply(gene_clusters, function(genes) {
              mapped <- suppressWarnings(
                clusterProfiler::bitr(genes, fromType = "SYMBOL",
                                      toType = "ENTREZID", OrgDb = orgdb)
              )
              if (nrow(mapped) > 0) unique(mapped$ENTREZID) else character(0)
            })
            entrez_clusters <- entrez_clusters[sapply(entrez_clusters, length) >= 5]
          }

          if (length(entrez_clusters) >= 2) {
            cc_react <- compareCluster(
              geneClusters  = entrez_clusters,
              fun           = "enrichPathway",
              organism      = reactome_org,
              pAdjustMethod = "BH",
              pvalueCutoff  = params$pvalue_cutoff,
              qvalueCutoff  = params$qvalue_cutoff,
              minGSSize     = params$min_gs_size,
              maxGSSize     = params$max_gs_size,
              readable      = TRUE
            )
            if (!is.null(cc_react) && nrow(as.data.frame(cc_react)) > 0) {
              p_cc <- dotplot(cc_react, showCategory = 10) +
                labs(title = sprintf("compareCluster Reactome — %s", dataset_id)) +
                theme_minimal(base_size = 11)
              .save_gg(p_cc, file.path(out_within, "compareCluster_Reactome"),
                       width = 14, height = 10)
              write.csv(as.data.frame(cc_react),
                        file.path(out_within, "compareCluster_Reactome.csv"),
                        row.names = FALSE)
            } else {
              cat("       Reactome: no significant terms.\n")
            }
          }
        }, error = function(e) {
          warning(sprintf("compareCluster Reactome: %s", conditionMessage(e)),
                  call. = FALSE)
        })

      } else {
        cat("       Fewer than 2 gene clusters with >= 5 genes — skipping compareCluster.\n")
      }
    }, error = function(e) {
      warning(sprintf("compareCluster failed for %s: %s", dataset_id,
                      conditionMessage(e)), call. = FALSE)
    })

    # ── Step 9: Condition-specific gene lists ────────────────────────────────
    if (length(contrast_ids) == 2) {
      cat("   [9] Condition-specific gene lists...\n")
      tryCatch({
        ct1 <- contrast_ids[1]
        ct2 <- contrast_ids[2]

        ct1_cfg <- get_contrast(cfg, dataset_id, ct1)
        ct2_cfg <- get_contrast(cfg, dataset_id, ct2)

        # Check if contrasts have different androgen_treatment flags
        andro_1 <- isTRUE(ct1_cfg$androgen_treatment)
        andro_2 <- isTRUE(ct2_cfg$androgen_treatment)

        if (andro_1 != andro_2) {
          # Identify which is the androgen contrast
          if (andro_2) {
            andro_ct   <- ct2
            noandro_ct <- ct1
            andro_label   <- ct2_cfg$label
            noandro_label <- ct1_cfg$label
          } else {
            andro_ct   <- ct1
            noandro_ct <- ct2
            andro_label   <- ct1_cfg$label
            noandro_label <- ct2_cfg$label
          }

          sig_andro   <- unique(c(sig_list[[andro_ct]]$up, sig_list[[andro_ct]]$down))
          sig_noandro <- unique(c(sig_list[[noandro_ct]]$up, sig_list[[noandro_ct]]$down))

          # Categories
          core_genes     <- intersect(sig_andro, sig_noandro)
          baseline_genes <- setdiff(sig_noandro, sig_andro)
          androgen_genes <- setdiff(sig_andro, sig_noandro)

          cat(sprintf("       Core genes (both):              %d\n", length(core_genes)))
          cat(sprintf("       Baseline genes (non-androgen):  %d\n", length(baseline_genes)))
          cat(sprintf("       Androgen-dependent genes:       %d\n", length(androgen_genes)))

          # Save gene lists as CSVs
          gene_category_df <- bind_rows(
            tibble(gene_symbol = core_genes,     category = "Core"),
            tibble(gene_symbol = baseline_genes, category = "Baseline"),
            tibble(gene_symbol = androgen_genes, category = "Androgen-dependent")
          )
          write.csv(gene_category_df,
                    file.path(out_within, "condition_specific_gene_lists.csv"),
                    row.names = FALSE)

          # Run GO enrichment on each category
          for (cat_name in c("Core", "Baseline", "Androgen-dependent")) {
            tryCatch({
              cat_genes <- gene_category_df$gene_symbol[gene_category_df$category == cat_name]
              if (length(cat_genes) < 5) {
                cat(sprintf("       %s: fewer than 5 genes — skipping enrichment.\n", cat_name))
                next
              }

              # Background: all genes present in both DEG tables
              bg_genes <- intersect(
                deg_list[[andro_ct]]$gene_symbol,
                deg_list[[noandro_ct]]$gene_symbol
              )

              ego <- enrichGO(
                gene         = cat_genes,
                universe     = bg_genes,
                OrgDb        = orgdb,
                ont          = "BP",
                keyType      = "SYMBOL",
                pAdjustMethod = "BH",
                pvalueCutoff = params$pvalue_cutoff,
                qvalueCutoff = params$qvalue_cutoff,
                minGSSize    = params$min_gs_size,
                maxGSSize    = params$max_gs_size,
                readable     = TRUE
              )

              safe_cat <- gsub("[^A-Za-z0-9_]", "_", cat_name)

              if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
                write.csv(as.data.frame(ego),
                          file.path(out_within, sprintf("condition_specific_%s_GO_BP.csv",
                                                        safe_cat)),
                          row.names = FALSE)

                n_show <- min(20, nrow(as.data.frame(ego)))
                p_ego <- dotplot(ego, showCategory = n_show) +
                  labs(title = sprintf("GO:BP — %s genes (%s)", cat_name, dataset_id)) +
                  theme_minimal(base_size = 11)
                .save_gg(p_ego,
                         file.path(out_within, sprintf("condition_specific_%s_GO_BP",
                                                       safe_cat)),
                         width = 11, height = 9)
                cat(sprintf("       %s enrichment: %d terms\n",
                            cat_name, nrow(as.data.frame(ego))))
              } else {
                cat(sprintf("       %s enrichment: no significant terms.\n", cat_name))
              }
            }, error = function(e) {
              warning(sprintf("Enrichment for %s genes: %s", cat_name,
                              conditionMessage(e)), call. = FALSE)
            })
          }

        } else {
          cat("       Both contrasts have same androgen_treatment flag.\n")
          cat("       Generating generic condition-specific gene lists instead.\n")

          # Generic condition-specific analysis (not androgen-dependent)
          sig_ct1 <- unique(c(sig_list[[ct1]]$up, sig_list[[ct1]]$down))
          sig_ct2 <- unique(c(sig_list[[ct2]]$up, sig_list[[ct2]]$down))

          shared_genes   <- intersect(sig_ct1, sig_ct2)
          only_ct1_genes <- setdiff(sig_ct1, sig_ct2)
          only_ct2_genes <- setdiff(sig_ct2, sig_ct1)

          ct1_label <- tryCatch(get_contrast(cfg, dataset_id, ct1)$label,
                                error = function(e) ct1)
          ct2_label <- tryCatch(get_contrast(cfg, dataset_id, ct2)$label,
                                error = function(e) ct2)

          cat(sprintf("       Shared DEGs:         %d\n", length(shared_genes)))
          cat(sprintf("       Unique to %s: %d\n", ct1, length(only_ct1_genes)))
          cat(sprintf("       Unique to %s: %d\n", ct2, length(only_ct2_genes)))

          if (length(shared_genes) > 0 || length(only_ct1_genes) > 0 || length(only_ct2_genes) > 0) {
            gene_category_df <- bind_rows(
              tibble(gene_symbol = shared_genes,   category = "Shared"),
              tibble(gene_symbol = only_ct1_genes, category = paste0("Unique_", ct1)),
              tibble(gene_symbol = only_ct2_genes, category = paste0("Unique_", ct2))
            )
            write.csv(gene_category_df,
                      file.path(out_within, "condition_specific_gene_lists.csv"),
                      row.names = FALSE)
          } else {
            cat("       No significant DEGs in either condition — no gene list written.\n")
          }
        }
      }, error = function(e) {
        warning(sprintf("Condition-specific analysis failed for %s: %s",
                        dataset_id, conditionMessage(e)), call. = FALSE)
      })
    } else {
      cat("   [9] Condition-specific analysis skipped (requires exactly 2 contrasts).\n")
    }

    cat(sprintf("\n   Completed within-dataset comparison for %s.\n", dataset_id))
  }

} else {
  cat("\n   No within-dataset comparisons configured — skipping.\n")
}


# ══════════════════════════════════════════════════════════════════════════════
# ACROSS-DATASET COMPARISONS
# ══════════════════════════════════════════════════════════════════════════════

across_comparisons <- cfg$cross_comparisons$across_datasets

if (!is.null(across_comparisons) && length(across_comparisons) > 0) {

  cat("\n════════════════════════════════════════════════════════════════════\n")
  cat("  ACROSS-DATASET COMPARISONS\n")
  cat("════════════════════════════════════════════════════════════════════\n")

  for (ac in across_comparisons) {

    comp_name <- ac$name
    pairs     <- ac$pairs

    cat(sprintf("\n── Comparison: %s ──\n", comp_name))

    out_across <- file.path(cross_dir, comp_name)
    dir.create(out_across, recursive = TRUE, showWarnings = FALSE)

    # ── Step 10: Load DEG tables for each pair ────────────────────────────
    cat("   [10] Loading DEG tables for each pair...\n")

    pair_degs  <- list()
    pair_sigs  <- list()
    pair_labels <- character(0)
    all_loaded  <- TRUE

    for (pr in pairs) {
      ds_id <- pr$dataset
      ct_id <- pr$contrast
      pair_key <- paste(ds_id, ct_id, sep = "/")

      deg_df <- .load_deg_rds(ds_id, ct_id)
      if (is.null(deg_df)) {
        cat(sprintf("       WARNING: Could not load %s — skipping comparison.\n", pair_key))
        all_loaded <- FALSE
        break
      }
      if (!"diffexpressed" %in% colnames(deg_df)) {
        deg_df <- classify_deg(deg_df, lfc_thresh = lfc_thresh, padj_thresh = padj_thresh)
      }
      pair_degs[[pair_key]]  <- deg_df
      pair_sigs[[pair_key]]  <- .extract_sig(deg_df)

      ct_label <- tryCatch(
        get_contrast(cfg, ds_id, ct_id)$label,
        error = function(e) pair_key
      )
      pair_labels <- c(pair_labels, setNames(ct_label, pair_key))

      cat(sprintf("       %s: %d genes (%d UP, %d DOWN)\n",
                  pair_key, nrow(deg_df),
                  length(pair_sigs[[pair_key]]$up),
                  length(pair_sigs[[pair_key]]$down)))
    }

    if (!all_loaded) next

    pair_keys <- names(pair_degs)

    # ── Step 11: DEG Overlap (UpSet/Venn) ────────────────────────────────
    cat("   [11] DEG overlap analysis...\n")
    tryCatch({
      upset_sets <- list()
      set_names  <- character(0)

      for (pk in pair_keys) {
        # Use short label for set names
        short_label <- gsub("/", "_", pk)
        up_name   <- paste0("UP_", short_label)
        down_name <- paste0("DOWN_", short_label)
        upset_sets[[up_name]]   <- pair_sigs[[pk]]$up
        upset_sets[[down_name]] <- pair_sigs[[pk]]$down
        set_names <- c(set_names, up_name, down_name)
      }

      upset_data <- UpSetR::fromList(upset_sets)

      if (nrow(upset_data) > 0) {
        .save_upset(upset_data, sets = set_names,
                    path_base = file.path(out_across, "deg_overlap_upset"))
      } else {
        cat("       No overlapping genes — UpSet skipped.\n")
      }

      # Also save overlap counts
      all_sig_sets <- lapply(pair_sigs, function(s) unique(c(s$up, s$down)))
      overlap_genes <- Reduce(intersect, all_sig_sets)
      cat(sprintf("       Overlap: %d genes significant in all pairs\n",
                  length(overlap_genes)))

      if (length(overlap_genes) > 0) {
        write.csv(data.frame(gene_symbol = sort(overlap_genes)),
                  file.path(out_across, "deg_overlap_shared_genes.csv"),
                  row.names = FALSE)
      }
    }, error = function(e) {
      warning(sprintf("DEG overlap failed for %s: %s", comp_name,
                      conditionMessage(e)), call. = FALSE)
    })

    # ── Step 12: Pathway Overlap (Jaccard similarity) ────────────────────
    cat("   [12] Pathway overlap (Jaccard similarity)...\n")
    tryCatch({
      # Load ORA results for each pair
      ora_term_sets <- list()

      for (pr in pairs) {
        ds_id <- pr$dataset
        ct_id <- pr$contrast
        pair_key <- paste(ds_id, ct_id, sep = "/")
        tables_path <- file.path("results", ds_id, ct_id, "tables")

        pair_terms <- character(0)

        # Try loading all ORA RDS files
        ora_files <- list.files(tables_path, pattern = "^ora_.*\\.rds$",
                                full.names = TRUE)
        for (ora_f in ora_files) {
          tryCatch({
            ora_res <- readRDS(ora_f)
            if (!is.null(ora_res)) {
              ora_df <- as.data.frame(ora_res)
              if (nrow(ora_df) > 0 && "Description" %in% colnames(ora_df)) {
                pair_terms <- c(pair_terms, ora_df$Description)
              } else if (nrow(ora_df) > 0 && "ID" %in% colnames(ora_df)) {
                pair_terms <- c(pair_terms, ora_df$ID)
              }
            }
          }, error = function(e) NULL)
        }

        ora_term_sets[[pair_key]] <- unique(pair_terms)
        cat(sprintf("       %s: %d enriched terms loaded\n",
                    pair_key, length(ora_term_sets[[pair_key]])))
      }

      # Compute Jaccard similarity matrix
      if (all(sapply(ora_term_sets, length) > 0)) {
        n_pairs <- length(ora_term_sets)
        jaccard_mat <- matrix(0, nrow = n_pairs, ncol = n_pairs,
                              dimnames = list(names(ora_term_sets),
                                              names(ora_term_sets)))

        for (i in seq_len(n_pairs)) {
          for (j in seq_len(n_pairs)) {
            set_i <- ora_term_sets[[i]]
            set_j <- ora_term_sets[[j]]
            inter <- length(intersect(set_i, set_j))
            union <- length(union(set_i, set_j))
            jaccard_mat[i, j] <- if (union > 0) inter / union else 0
          }
        }

        # Plot heatmap
        jaccard_long <- as.data.frame(as.table(jaccard_mat)) %>%
          dplyr::rename(Pair1 = Var1, Pair2 = Var2, Jaccard = Freq)

        p_jaccard <- ggplot(jaccard_long, aes(x = Pair1, y = Pair2, fill = Jaccard)) +
          geom_tile(colour = "white", linewidth = 0.5) +
          geom_text(aes(label = sprintf("%.2f", Jaccard)), size = 4) +
          scale_fill_gradient2(
            low      = "white",
            mid      = "#FEE08B",
            high     = "#D73027",
            midpoint = 0.25,
            limits   = c(0, 1),
            name     = "Jaccard\nSimilarity"
          ) +
          labs(
            title = sprintf("Pathway Overlap — %s", comp_name),
            x = NULL, y = NULL
          ) +
          get_theme_publication() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
                axis.text.y = element_text(size = 9))

        .save_gg(p_jaccard, file.path(out_across, "pathway_jaccard_heatmap"),
                 width = 8, height = 7)

        # Save matrix
        write.csv(as.data.frame(jaccard_mat),
                  file.path(out_across, "pathway_jaccard_matrix.csv"))
      } else {
        cat("       Some pairs have no ORA results — Jaccard heatmap skipped.\n")
      }
    }, error = function(e) {
      warning(sprintf("Pathway overlap failed for %s: %s", comp_name,
                      conditionMessage(e)), call. = FALSE)
    })

    # ── Step 13: Purinergic Signature Comparison ─────────────────────────
    cat("   [13] Purinergic signature comparison...\n")
    tryCatch({
      purinergic_frames <- list()

      for (pr in pairs) {
        ds_id <- pr$dataset
        ct_id <- pr$contrast
        pair_key <- paste(ds_id, ct_id, sep = "/")
        pur_path <- file.path("results", ds_id, ct_id, "tables",
                              "purinergic_gene_summary.csv")

        if (file.exists(pur_path)) {
          pur_df <- read.csv(pur_path, stringsAsFactors = FALSE)
          if ("category" %in% colnames(pur_df)) {
            # Count significant genes per category
            sig_col_name <- intersect(c("sig", "significant"), colnames(pur_df))
            if (length(sig_col_name) > 0) {
              pur_summary <- pur_df %>%
                filter(!!sym(sig_col_name[1]) == TRUE) %>%
                count(category, name = "n_sig") %>%
                mutate(pair = pair_key)
            } else {
              # Fallback: count rows per category
              pur_summary <- pur_df %>%
                count(category, name = "n_sig") %>%
                mutate(pair = pair_key)
            }
            purinergic_frames[[pair_key]] <- pur_summary
          }
          cat(sprintf("       %s: purinergic data loaded\n", pair_key))
        } else {
          cat(sprintf("       %s: purinergic_gene_summary.csv not found — skipping.\n",
                      pair_key))
        }
      }

      if (length(purinergic_frames) >= 2) {
        pur_combined <- bind_rows(purinergic_frames)

        p_pur <- ggplot(pur_combined, aes(x = category, y = n_sig, fill = pair)) +
          geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.85) +
          labs(
            title = sprintf("Purinergic Sig. Genes — %s", comp_name),
            x     = NULL,
            y     = "Number of Significant Genes",
            fill  = "Dataset/Contrast"
          ) +
          get_theme_publication() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
                legend.position = "right")

        .save_gg(p_pur, file.path(out_across, "purinergic_comparison"),
                 width = 12, height = 7)
      } else {
        cat("       Insufficient purinergic data for comparison.\n")
      }
    }, error = function(e) {
      warning(sprintf("Purinergic comparison failed for %s: %s", comp_name,
                      conditionMessage(e)), call. = FALSE)
    })

    # ── Step 14: Mitochondrial Signature Comparison ──────────────────────
    cat("   [14] Mitochondrial signature comparison...\n")
    tryCatch({
      mito_frames <- list()

      for (pr in pairs) {
        ds_id <- pr$dataset
        ct_id <- pr$contrast
        pair_key <- paste(ds_id, ct_id, sep = "/")
        mito_path <- file.path("results", ds_id, ct_id, "tables",
                               "mitochondrial_gene_summary.csv")

        if (file.exists(mito_path)) {
          mito_df <- read.csv(mito_path, stringsAsFactors = FALSE)
          if ("category" %in% colnames(mito_df)) {
            sig_col_name <- intersect(c("sig", "significant"), colnames(mito_df))
            if (length(sig_col_name) > 0) {
              mito_summary <- mito_df %>%
                filter(!!sym(sig_col_name[1]) == TRUE) %>%
                count(category, name = "n_sig") %>%
                mutate(pair = pair_key)
            } else {
              mito_summary <- mito_df %>%
                count(category, name = "n_sig") %>%
                mutate(pair = pair_key)
            }
            mito_frames[[pair_key]] <- mito_summary
          }
          cat(sprintf("       %s: mitochondrial data loaded\n", pair_key))
        } else {
          cat(sprintf("       %s: mitochondrial_gene_summary.csv not found — skipping.\n",
                      pair_key))
        }
      }

      if (length(mito_frames) >= 2) {
        mito_combined <- bind_rows(mito_frames)

        # Main category barplot
        p_mito <- ggplot(mito_combined, aes(x = category, y = n_sig, fill = pair)) +
          geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.85) +
          labs(
            title = sprintf("Mitochondrial Sig. Genes — %s", comp_name),
            x     = NULL,
            y     = "Number of Significant Genes",
            fill  = "Dataset/Contrast"
          ) +
          get_theme_publication() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
                legend.position = "right")

        .save_gg(p_mito, file.path(out_across, "mitochondrial_comparison"),
                 width = 14, height = 7)

        # OXPHOS complex comparison (filter to complex_ categories)
        oxphos_combined <- mito_combined %>%
          filter(grepl("^complex_", category, ignore.case = TRUE))

        if (nrow(oxphos_combined) > 0) {
          p_oxphos <- ggplot(oxphos_combined, aes(x = category, y = n_sig, fill = pair)) +
            geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.85) +
            labs(
              title = sprintf("OXPHOS Complex Comparison — %s", comp_name),
              x     = NULL,
              y     = "Number of Significant Genes",
              fill  = "Dataset/Contrast"
            ) +
            get_theme_publication() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
                  legend.position = "right")

          .save_gg(p_oxphos, file.path(out_across, "oxphos_complex_comparison"),
                   width = 10, height = 7)
        }
      } else {
        cat("       Insufficient mitochondrial data for comparison.\n")
      }
    }, error = function(e) {
      warning(sprintf("Mitochondrial comparison failed for %s: %s", comp_name,
                      conditionMessage(e)), call. = FALSE)
    })

    # ── Step 15: Direction Concordance Heatmap ───────────────────────────
    cat("   [15] Direction concordance heatmap...\n")
    tryCatch({
      # Get curated gene lists
      purinergic_genes <- get_purinergic_genes()
      mito_genes       <- get_mitochondrial_genes()

      curated_genes <- bind_rows(
        purinergic_genes %>% dplyr::select(gene_symbol, category) %>%
          mutate(gene_class = "Purinergic"),
        mito_genes %>% dplyr::select(gene_symbol, category) %>%
          mutate(gene_class = "Mitochondrial")
      ) %>%
        distinct(gene_symbol, .keep_all = TRUE)

      # For each pair, get direction of curated genes
      direction_frames <- list()

      for (pk in pair_keys) {
        deg_df <- pair_degs[[pk]]
        merged <- curated_genes %>%
          left_join(
            deg_df %>% dplyr::select(gene_symbol, log2fc, padj),
            by = "gene_symbol"
          ) %>%
          mutate(
            sig = !is.na(padj) & padj < padj_thresh & abs(log2fc) > lfc_thresh,
            direction = case_when(
              sig & log2fc > 0 ~ 1,    # Up
              sig & log2fc < 0 ~ -1,   # Down
              TRUE             ~ 0      # NS / not detected
            ),
            pair = pk
          )
        direction_frames[[pk]] <- merged %>%
          dplyr::select(gene_symbol, category, gene_class, direction, pair)
      }

      dir_combined <- bind_rows(direction_frames)

      # Pivot to wide: genes x pairs
      dir_wide <- dir_combined %>%
        dplyr::select(gene_symbol, category, gene_class, pair, direction) %>%
        pivot_wider(names_from = pair, values_from = direction, values_fill = 0)

      # Filter to genes that are significant in at least one pair
      dir_cols <- pair_keys
      dir_wide <- dir_wide %>%
        filter(if_any(all_of(dir_cols), ~ . != 0))

      if (nrow(dir_wide) > 0) {
        # Pivot back to long for ggplot
        dir_long <- dir_wide %>%
          pivot_longer(
            cols = all_of(dir_cols),
            names_to = "pair",
            values_to = "direction"
          ) %>%
          mutate(
            direction_label = case_when(
              direction ==  1 ~ "Up",
              direction == -1 ~ "Down",
              TRUE            ~ "NS"
            ),
            direction_label = factor(direction_label, levels = c("Down", "NS", "Up"))
          )

        # Order genes by gene_class then category
        gene_order <- dir_wide %>%
          arrange(gene_class, category) %>%
          pull(gene_symbol)
        dir_long$gene_symbol <- factor(dir_long$gene_symbol,
                                       levels = rev(gene_order))

        p_dir <- ggplot(dir_long, aes(x = pair, y = gene_symbol,
                                       fill = direction_label)) +
          geom_tile(colour = "white", linewidth = 0.2) +
          scale_fill_manual(
            values = c("Up" = colors$up, "Down" = colors$down, "NS" = "grey90"),
            name   = "Direction"
          ) +
          facet_grid(gene_class ~ ., scales = "free_y", space = "free_y") +
          labs(
            title = sprintf("Direction Concordance — %s", comp_name),
            x     = NULL,
            y     = NULL
          ) +
          get_theme_publication(base_size = 10) +
          theme(
            axis.text.x  = element_text(angle = 45, hjust = 1, size = 8),
            axis.text.y  = element_text(size = 6),
            strip.text.y = element_text(angle = 0, size = 10),
            legend.position = "bottom"
          )

        # Dynamic height
        n_genes <- length(unique(dir_long$gene_symbol))
        fig_height <- max(8, n_genes * 0.15 + 3)

        .save_gg(p_dir, file.path(out_across, "direction_concordance_heatmap"),
                 width = 10, height = fig_height)

        # Save table
        write.csv(dir_wide,
                  file.path(out_across, "direction_concordance_table.csv"),
                  row.names = FALSE)

        cat(sprintf("       %d curated genes with differential expression in >= 1 pair\n",
                    n_genes))
      } else {
        cat("       No curated genes are significant in any pair — heatmap skipped.\n")
      }
    }, error = function(e) {
      warning(sprintf("Direction concordance failed for %s: %s", comp_name,
                      conditionMessage(e)), call. = FALSE)
    })

    cat(sprintf("\n   Completed across-dataset comparison: %s\n", comp_name))
  }

} else {
  cat("\n   No across-dataset comparisons configured — skipping.\n")
}


# ══════════════════════════════════════════════════════════════════════════════
# STEP 16: Write completion flag
# ══════════════════════════════════════════════════════════════════════════════

cat("\n[16] Writing completion flag...\n")
flag_path <- file.path(cross_dir, "done.flag")
writeLines(
  sprintf("cross_comparison completed: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  flag_path
)
cat(sprintf("       Written: %s\n", flag_path))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 17: Print summary
# ══════════════════════════════════════════════════════════════════════════════

cat("\n════════════════════════════════════════════════════════════════════\n")
cat("  Step 10 complete — Cross-Comparison Analyses\n")
cat("────────────────────────────────────────────────────────────────────\n")

n_within  <- if (!is.null(within_comparisons)) length(within_comparisons) else 0
n_across  <- if (!is.null(across_comparisons)) length(across_comparisons) else 0

cat(sprintf("  Within-dataset comparisons : %d\n", n_within))
if (n_within > 0) {
  for (wc in within_comparisons) {
    cat(sprintf("    - %s (%s)\n", wc$dataset,
                paste(unlist(wc$contrasts), collapse = " vs ")))
  }
}

cat(sprintf("  Across-dataset comparisons : %d\n", n_across))
if (n_across > 0) {
  for (ac in across_comparisons) {
    pair_desc <- sapply(ac$pairs, function(p) paste(p$dataset, p$contrast, sep = "/"))
    cat(sprintf("    - %s (%s)\n", ac$name, paste(pair_desc, collapse = " vs ")))
  }
}

cat(sprintf("  Output directory : %s\n", cross_dir))
cat("════════════════════════════════════════════════════════════════════\n")

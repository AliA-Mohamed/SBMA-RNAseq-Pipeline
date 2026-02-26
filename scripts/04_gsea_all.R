#!/usr/bin/env Rscript
# ══════════════════════════════════════════════════════════════════════════════
# 04_gsea_all.R — Run GSEA across all configured databases for one contrast
# ══════════════════════════════════════════════════════════════════════════════
#
# Usage:
#   Rscript scripts/04_gsea_all.R \
#     --config config/config.yaml \
#     --dataset SBMA_iPSC_Q51 \
#     --contrast disease_vs_control_vehicle
#
# Description:
#   Runs Gene Set Enrichment Analysis (GSEA) for a single dataset/contrast
#   pair across all databases configured in config.yaml.  Requires that
#   step 01 has already produced the ranked list RDS files.  If the contrast
#   is not full_ranked, the script exits gracefully with a skip notice.
#
# Inputs (from results/<dataset>/<contrast>/data/):
#   ranked_symbol.rds  — Named numeric vector (symbol -> metric), descending
#   ranked_entrez.rds  — Named numeric vector (Entrez -> metric), descending
#
# Outputs (under results/<dataset>/<contrast>/):
#   tables/gsea_<DB>.rds         — gseaResult objects
#   tables/gsea_<DB>.csv         — Tabular GSEA results
#   plots/gsea/gsea_<DB>_dotplot.{pdf,png}    — Dotplots (top 20 by NES)
#   plots/gsea/gsea_<DB>_ridgeplot.{pdf,png}  — Ridge plots
#   plots/gsea/gsea_<DB>_running_<i>.{pdf,png} — Running score plots (top 5)
#   plots/gsea/gsea_nes_heatmap.{pdf,png}     — NES heatmap across databases
#   tables/gsea_done.flag                      — Completion flag
#   tables/gsea_SKIPPED.txt                    — Written if full_ranked=FALSE
# ══════════════════════════════════════════════════════════════════════════════

# ── Load libraries ───────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(ggplot2)
  library(clusterProfiler)
  library(enrichplot)
  library(ReactomePA)
  library(msigdbr)
})

# ── Resolve project root and source utility modules ──────────────────────────
script_dir   <- dirname(sub("--file=", "",
                             grep("--file=", commandArgs(trailingOnly = FALSE),
                                  value = TRUE)))
project_root <- dirname(script_dir)

source(file.path(project_root, "lib", "config_utils.R"))
source(file.path(project_root, "lib", "data_utils.R"))
source(file.path(project_root, "lib", "enrichment_utils.R"))
source(file.path(project_root, "lib", "plot_utils.R"))

# ── Parse CLI arguments ──────────────────────────────────────────────────────
args <- parse_pipeline_args()

cat("════════════════════════════════════════════════════════════════════\n")
cat("  SBMA RNAseq Pipeline — Step 04: GSEA (All Databases)\n")
cat("════════════════════════════════════════════════════════════════════\n")
cat(sprintf("  Config   : %s\n", args$config))
cat(sprintf("  Dataset  : %s\n", args$dataset))
cat(sprintf("  Contrast : %s\n", args$contrast))
cat("────────────────────────────────────────────────────────────────────\n\n")

# ── 1. Load and validate configuration ───────────────────────────────────────
cat("[1/8] Loading configuration...\n")
cfg <- load_config(args$config)
cat(sprintf("       Species        : %s\n", cfg$species))
cat(sprintf("       OrgDb          : %s\n", cfg$orgdb))
cat(sprintf("       Databases      : %s\n", paste(cfg$databases, collapse = ", ")))

# ── 2. Extract contrast info and build output paths ──────────────────────────
cat("\n[2/8] Extracting contrast configuration...\n")
contrast_cfg <- get_contrast(cfg, args$dataset, args$contrast)
out_dir      <- get_output_dir(cfg, args$dataset, args$contrast)
data_dir     <- file.path(out_dir, "data")
tables_dir   <- file.path(out_dir, "tables")
plots_dir    <- file.path(out_dir, "plots", "gsea")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir,  recursive = TRUE, showWarnings = FALSE)

cat(sprintf("       Label        : %s\n", contrast_cfg$label))
cat(sprintf("       Full ranked  : %s\n", contrast_cfg$full_ranked))
cat(sprintf("       Output base  : %s\n", out_dir))

# ── 3. GUARD: Check full_ranked ──────────────────────────────────────────────
cat("\n[3/8] Checking full_ranked flag...\n")

if (!isTRUE(contrast_cfg$full_ranked)) {

  skip_msg <- paste0(
    "GSEA SKIPPED for ", args$dataset, "/", args$contrast, "\n\n",
    "Reason: contrast_cfg$full_ranked is FALSE.\n",
    "GSEA requires a full ranked gene list (all genes with log2FC and p-values),\n",
    "not a pre-filtered subset. This contrast provides only filtered DEGs.\n\n",
    "To enable GSEA for this contrast, set 'full_ranked: true' in config.yaml\n",
    "and ensure the DEG file contains the complete unfiltered results table.\n"
  )

  # Write skip notice
  writeLines(skip_msg, file.path(tables_dir, "gsea_SKIPPED.txt"))
  cat(sprintf("       WARNING: %s\n", "full_ranked = FALSE — GSEA requires full ranked list"))
  cat(sprintf("       Written: %s\n", file.path(tables_dir, "gsea_SKIPPED.txt")))

  # Write completion flag
  writeLines(paste("gsea_done — SKIPPED at", Sys.time()),
             file.path(tables_dir, "gsea_done.flag"))
  cat(sprintf("       Written: %s\n", file.path(tables_dir, "gsea_done.flag")))

  cat("\n════════════════════════════════════════════════════════════════════\n")
  cat("  Step 04 complete (SKIPPED — full_ranked = FALSE).\n")
  cat("════════════════════════════════════════════════════════════════════\n")
  quit(save = "no", status = 0)
}

cat("       full_ranked = TRUE — proceeding with GSEA\n")

# ── 4. Load ranked lists ─────────────────────────────────────────────────────
cat("\n[4/8] Loading ranked lists...\n")

ranked_symbol_path <- file.path(data_dir, "ranked_symbol.rds")
ranked_entrez_path <- file.path(data_dir, "ranked_entrez.rds")

if (!file.exists(ranked_symbol_path)) {
  stop(sprintf("Ranked symbol list not found: %s\nRun step 01 first.", ranked_symbol_path),
       call. = FALSE)
}
if (!file.exists(ranked_entrez_path)) {
  stop(sprintf("Ranked Entrez list not found: %s\nRun step 01 first.", ranked_entrez_path),
       call. = FALSE)
}

ranked_symbol <- readRDS(ranked_symbol_path)
ranked_entrez <- readRDS(ranked_entrez_path)

cat(sprintf("       Ranked symbols : %d genes\n", length(ranked_symbol)))
cat(sprintf("       Ranked Entrez  : %d genes\n", length(ranked_entrez)))
cat(sprintf("       Metric range   : [%.2f, %.2f]\n",
            min(ranked_symbol), max(ranked_symbol)))

# ── 5. Run GSEA across all databases ────────────────────────────────────────
cat("\n[5/8] Running GSEA analyses...\n")

params    <- cfg$enrichment
orgdb     <- cfg$orgdb
kegg_org  <- cfg$kegg_organism
react_org <- cfg$reactome_organism
plot_fmts <- cfg$plots$format %||% c("pdf", "png")

# Resolve msigdbr species for Hallmarks / WikiPathways
species_full <- switch(tolower(cfg$species),
  human = "Homo sapiens",
  mouse = "Mus musculus",
  stop(sprintf("Unsupported species for msigdbr: '%s'", cfg$species), call. = FALSE)
)

# Container for all GSEA results
gsea_results <- list()

for (db in cfg$databases) {

  cat(sprintf("\n  --- %s ---\n", db))

  set.seed(cfg$enrichment$seed)

  res <- tryCatch(
    switch(db,
      GO_BP        = run_gsea_go(ranked_symbol, "BP", orgdb, params),
      GO_CC        = run_gsea_go(ranked_symbol, "CC", orgdb, params),
      GO_MF        = run_gsea_go(ranked_symbol, "MF", orgdb, params),
      KEGG         = run_gsea_kegg(ranked_entrez, kegg_org, params),
      Reactome     = run_gsea_reactome(ranked_entrez, react_org, params),
      Hallmarks    = run_gsea_hallmarks(ranked_symbol, species_full, params),
      WikiPathways = {
        cat("       Attempting WikiPathways GSEA via msigdbr...\n")
        tryCatch({
          wp_sets <- msigdbr(species = species_full, category = "C2",
                             subcategory = "CP:WIKIPATHWAYS")
          wp_t2g  <- wp_sets[, c("gs_name", "gene_symbol")]
          set.seed(params$seed)
          GSEA(
            geneList      = ranked_symbol,
            TERM2GENE     = wp_t2g,
            pAdjustMethod = "BH",
            pvalueCutoff  = params$pvalue_cutoff,
            minGSSize     = params$min_gs_size,
            maxGSSize     = params$max_gs_size,
            nPermSimple   = params$n_perm
          )
        }, error = function(e) {
          message(sprintf("       WikiPathways GSEA: %s", conditionMessage(e)))
          NULL
        })
      },
      {
        message(sprintf("       Unknown database '%s' — skipping.", db))
        NULL
      }
    ),
    error = function(e) {
      message(sprintf("       %s GSEA failed: %s", db, conditionMessage(e)))
      NULL
    }
  )

  gsea_results[[db]] <- res

  if (is.null(res)) {
    cat(sprintf("       Result: NULL (skipped or no significant terms)\n"))
  } else {
    n_terms <- nrow(as.data.frame(res))
    cat(sprintf("       Result: %d enriched terms\n", n_terms))
  }
}

# ── 5b. Annotate artifact gene families in GSEA results ────────────────────
cat("\n[5b/8] Annotating artifact gene families...\n")

artifact_patterns <- c("^OR[0-9]", "^TAS[12]R", "^KRT[0-9]", "^KRTAP")
artifact_regex <- paste(artifact_patterns, collapse = "|")

for (db in names(gsea_results)) {
  res <- gsea_results[[db]]
  if (is.null(res)) next

  res_df <- as.data.frame(res)
  if (nrow(res_df) == 0) next

  # GSEA uses core_enrichment column (slash-separated gene list)
  gene_col <- if ("core_enrichment" %in% colnames(res_df)) "core_enrichment"
              else if ("geneID" %in% colnames(res_df)) "geneID"
              else next

  res_df$artifact_gene_count <- sapply(res_df[[gene_col]], function(gid) {
    genes <- unlist(strsplit(gid, "/"))
    sum(grepl(artifact_regex, genes, ignore.case = FALSE))
  })
  res_df$artifact_gene_fraction <- round(
    res_df$artifact_gene_count / sapply(res_df[[gene_col]], function(gid) {
      length(unlist(strsplit(gid, "/")))
    }), 3
  )
  res_df$artifact_flag <- ifelse(
    res_df$artifact_gene_fraction >= 0.5,
    "ARTIFACT_DOMINATED",
    ifelse(res_df$artifact_gene_count > 0, "HAS_ARTIFACTS", "")
  )

  n_dominated <- sum(res_df$artifact_flag == "ARTIFACT_DOMINATED")
  if (n_dominated > 0) {
    cat(sprintf("       %s: %d artifact-dominated terms flagged\n", db, n_dominated))
  }

  # Store annotated df for saving later
  attr(gsea_results[[db]], "annotated_df") <- res_df
}

# ── 6. Save results and generate plots ──────────────────────────────────────
cat("\n[6/8] Saving results and generating plots...\n")

for (db in names(gsea_results)) {

  res <- gsea_results[[db]]
  if (is.null(res)) next

  res_df <- as.data.frame(res)
  if (nrow(res_df) == 0) {
    cat(sprintf("       %s: 0 terms — skipping save/plots\n", db))
    next
  }

  db_label <- gsub("_", "", db)
  rds_name <- paste0("gsea_", db)

  # -- Save RDS + CSV via save_enrichment --
  save_enrichment(res, tables_dir, rds_name, formats = plot_fmts)

  # Overwrite CSV with artifact-annotated version if available
  annotated_df <- attr(res, "annotated_df")
  if (!is.null(annotated_df) && nrow(annotated_df) > 0) {
    csv_path <- file.path(tables_dir, paste0(rds_name, ".csv"))
    write.csv(annotated_df, csv_path, row.names = FALSE)
  }

  cat(sprintf("       Saved: %s (.rds, .csv, dotplot)\n", rds_name))

  # -- Dotplot (top 20 by NES) --
  tryCatch({
    n_show <- min(20, nrow(res_df))
    p_dot <- dotplot(res, showCategory = n_show, orderBy = "NES") +
      labs(title = paste("GSEA", db, "-", contrast_cfg$label)) +
      theme_minimal(base_size = 12)

    for (fmt in plot_fmts) {
      dot_path <- file.path(plots_dir, paste0("gsea_", db, "_dotplot.", fmt))
      ggsave(dot_path, plot = p_dot, width = 12, height = 10, dpi = 300)
    }
    cat(sprintf("       Saved: gsea_%s_dotplot\n", db))
  }, error = function(e) {
    message(sprintf("       Dotplot [%s]: %s", db, conditionMessage(e)))
  })

  # -- Ridge plot --
  tryCatch({
    n_show <- min(20, nrow(res_df))
    p_ridge <- ridgeplot(res, showCategory = n_show) +
      labs(title = paste("GSEA Ridge:", db, "-", contrast_cfg$label)) +
      theme_minimal(base_size = 11)

    for (fmt in plot_fmts) {
      ridge_path <- file.path(plots_dir, paste0("gsea_", db, "_ridgeplot.", fmt))
      ggsave(ridge_path, plot = p_ridge, width = 12, height = 10, dpi = 300)
    }
    cat(sprintf("       Saved: gsea_%s_ridgeplot\n", db))
  }, error = function(e) {
    message(sprintf("       Ridgeplot [%s]: %s", db, conditionMessage(e)))
  })

  # -- Running score plots (top 5 by absolute NES) --
  tryCatch({
    res_df_sorted <- res_df %>%
      arrange(desc(abs(NES)))

    n_running <- min(5, nrow(res_df_sorted))

    for (i in seq_len(n_running)) {
      term_idx <- which(res_df$ID == res_df_sorted$ID[i])
      term_label <- res_df_sorted$Description[i]
      if (is.na(term_label) || term_label == "") {
        term_label <- res_df_sorted$ID[i]
      }

      p_run <- gseaplot2(res, geneSetID = term_idx,
                          title = term_label,
                          pvalue_table = TRUE)

      # Sanitize term label for filename
      safe_label <- gsub("[^A-Za-z0-9_]", "_", substr(term_label, 1, 60))
      safe_label <- gsub("_+", "_", safe_label)
      safe_label <- gsub("^_|_$", "", safe_label)

      for (fmt in plot_fmts) {
        run_path <- file.path(plots_dir,
                              paste0("gsea_", db, "_running_", i, "_",
                                     safe_label, ".", fmt))
        ggsave(run_path, plot = p_run, width = 10, height = 7, dpi = 300)
      }
    }
    cat(sprintf("       Saved: %d running score plots for %s\n", n_running, db))
  }, error = function(e) {
    message(sprintf("       Running score [%s]: %s", db, conditionMessage(e)))
  })
}

# ── 7. NES heatmap across all databases ──────────────────────────────────────
cat("\n[7/8] Creating NES heatmap across databases...\n")

tryCatch({
  # Collect top 10 terms from each database
  nes_frames <- list()

  for (db in names(gsea_results)) {
    res <- gsea_results[[db]]
    if (is.null(res)) next

    res_df <- as.data.frame(res)
    if (nrow(res_df) == 0) next

    top_n <- min(10, nrow(res_df))
    top_terms <- res_df %>%
      arrange(desc(abs(NES))) %>%
      head(top_n) %>%
      mutate(
        Database = db,
        Term     = ifelse(is.na(Description) | Description == "", ID, Description)
      ) %>%
      dplyr::select(Database, Term, NES, p.adjust)

    nes_frames[[db]] <- top_terms
  }

  if (length(nes_frames) > 0) {
    nes_combined <- bind_rows(nes_frames)

    # Truncate long term names for display
    nes_combined <- nes_combined %>%
      mutate(
        Term_short = ifelse(nchar(Term) > 55,
                            paste0(substr(Term, 1, 52), "..."),
                            Term),
        Term_label = paste0("[", Database, "] ", Term_short),
        sig_star   = case_when(
          p.adjust < 0.001 ~ "***",
          p.adjust < 0.01  ~ "**",
          p.adjust < 0.05  ~ "*",
          TRUE             ~ ""
        )
      )

    # Order terms by NES within each database group
    nes_combined <- nes_combined %>%
      arrange(Database, desc(NES)) %>%
      mutate(Term_label = factor(Term_label, levels = rev(unique(Term_label))))

    # Determine symmetric color limits
    nes_max <- max(abs(nes_combined$NES), na.rm = TRUE)

    colors <- default_plot_colors()

    p_heatmap <- ggplot(nes_combined, aes(x = Database, y = Term_label, fill = NES)) +
      geom_tile(colour = "white", linewidth = 0.3) +
      geom_text(aes(label = sig_star), size = 3, vjust = 0.8) +
      scale_fill_gradient2(
        low      = colors$heatmap_low,
        mid      = colors$heatmap_mid,
        high     = colors$heatmap_high,
        midpoint = 0,
        limits   = c(-nes_max, nes_max),
        name     = "NES"
      ) +
      labs(
        title = paste("GSEA NES Heatmap -", contrast_cfg$label),
        x     = NULL,
        y     = NULL
      ) +
      get_theme_publication(base_size = 11) +
      theme(
        axis.text.x  = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y  = element_text(size = 7),
        panel.grid    = element_blank(),
        plot.title    = element_text(size = 13)
      )

    # Dynamic height based on number of terms
    fig_height <- max(8, nrow(nes_combined) * 0.28 + 2)
    fig_width  <- max(8, length(unique(nes_combined$Database)) * 1.5 + 4)

    for (fmt in plot_fmts) {
      hm_path <- file.path(plots_dir, paste0("gsea_nes_heatmap.", fmt))
      ggsave(hm_path, plot = p_heatmap, width = fig_width, height = fig_height,
             dpi = 300, bg = "white")
    }
    cat(sprintf("       Saved: gsea_nes_heatmap (%d terms across %d databases)\n",
                nrow(nes_combined), length(unique(nes_combined$Database))))
  } else {
    cat("       No databases returned significant terms — heatmap skipped.\n")
  }
}, error = function(e) {
  message(sprintf("       NES heatmap: %s", conditionMessage(e)))
})

# ── 8. Write flag file and print summary ─────────────────────────────────────
cat("\n[8/8] Writing completion flag...\n")
writeLines(paste("gsea_done at", Sys.time()),
           file.path(tables_dir, "gsea_done.flag"))
cat(sprintf("       Written: %s\n", file.path(tables_dir, "gsea_done.flag")))

# ── Summary ──────────────────────────────────────────────────────────────────
cat("\n════════════════════════════════════════════════════════════════════\n")
cat("  Step 04 complete — GSEA (All Databases)\n")
cat("────────────────────────────────────────────────────────────────────\n")
cat(sprintf("  Dataset  : %s\n", args$dataset))
cat(sprintf("  Contrast : %s (%s)\n", args$contrast, contrast_cfg$label))
cat("  Results by database:\n")

for (db in names(gsea_results)) {
  res <- gsea_results[[db]]
  if (is.null(res)) {
    cat(sprintf("    %-15s : NULL (skipped/failed)\n", db))
  } else {
    n_terms <- nrow(as.data.frame(res))
    n_sig   <- sum(as.data.frame(res)$p.adjust < 0.05, na.rm = TRUE)
    cat(sprintf("    %-15s : %d terms (%d with padj < 0.05)\n", db, n_terms, n_sig))
  }
}

cat(sprintf("  Output   : %s\n", out_dir))
cat("════════════════════════════════════════════════════════════════════\n")

#!/usr/bin/env Rscript
# ==============================================================================
# 03_ora_all.R
# Over-Representation Analysis (ORA) across all configured databases
#
# Runs ORA on UP and DOWN regulated gene sets separately for each database
# specified in the pipeline config, saves results and generates combined
# barplots showing top enriched terms from both directions.
#
# Usage:
#   Rscript scripts/03_ora_all.R --config config/config.yaml \
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

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# ---- Parse CLI arguments and load config -------------------------------------

args <- parse_pipeline_args()

cfg <- load_config(args$config)

contrast_cfg <- get_contrast(cfg, args$dataset, args$contrast)
out_dir      <- get_output_dir(cfg, args$dataset, args$contrast)

ds <- args$dataset
ct <- args$contrast

message("=== 03_ora_all.R ===")
message(sprintf("Dataset:  %s", ds))
message(sprintf("Contrast: %s (%s)", ct, contrast_cfg$label))

# ---- Derive paths ------------------------------------------------------------

data_dir   <- file.path(out_dir, "data")
tables_dir <- file.path(out_dir, "tables")
plots_dir  <- file.path(out_dir, "plots", "ora")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir,  recursive = TRUE, showWarnings = FALSE)

# ---- Load prepped data -------------------------------------------------------

deg_splits_path <- file.path(data_dir, "deg_splits.rds")
deg_all_path    <- file.path(data_dir, "deg_all.rds")

if (!file.exists(deg_splits_path)) {
  stop(sprintf("deg_splits.rds not found: %s\nRun the prep step first.",
               deg_splits_path), call. = FALSE)
}
if (!file.exists(deg_all_path)) {
  stop(sprintf("deg_all.rds not found: %s\nRun the prep step first.",
               deg_all_path), call. = FALSE)
}

deg_splits <- readRDS(deg_splits_path)
deg_all    <- readRDS(deg_all_path)

message(sprintf("Loaded %d UP genes, %d DOWN genes, %d total genes.",
                length(deg_splits$up), length(deg_splits$down),
                nrow(deg_all)))

# ---- Species parameters ------------------------------------------------------

sp_params   <- get_species_params(cfg)
orgdb       <- sp_params$orgdb
kegg_org    <- sp_params$kegg_organism
reactome_org <- sp_params$reactome_organism

# Load the OrgDb annotation package
if (!requireNamespace(orgdb, quietly = TRUE)) {
  stop(sprintf("Annotation package '%s' is not installed. "
               , "Install with BiocManager::install('%s').",
               orgdb, orgdb), call. = FALSE)
}
library(orgdb, character.only = TRUE)

# Map config species to msigdbr species name
species_full <- switch(tolower(cfg$species),
  human = "Homo sapiens",
  mouse = "Mus musculus",
  stop(sprintf("Unsupported species '%s' for msigdbr. "
               , "Expected 'human' or 'mouse'.", cfg$species), call. = FALSE)
)

# ---- Background genes --------------------------------------------------------

background_genes <- unique(deg_all$gene_symbol)
message(sprintf("Background universe: %d gene symbols.", length(background_genes)))

# ---- Convert symbols to Entrez IDs ------------------------------------------

needs_entrez <- any(c("KEGG", "Reactome") %in% cfg$databases)

entrez_up    <- character(0)
entrez_down  <- character(0)
entrez_univ  <- character(0)

if (needs_entrez) {
  message("Converting gene symbols to Entrez IDs ...")
  entrez_up   <- .symbols_to_entrez(deg_splits$up,   orgdb)
  entrez_down <- .symbols_to_entrez(deg_splits$down,  orgdb)
  entrez_univ <- .symbols_to_entrez(background_genes, orgdb)
  message(sprintf("Entrez mapping: UP=%d, DOWN=%d, universe=%d",
                  length(entrez_up), length(entrez_down), length(entrez_univ)))
}

# ---- Enrichment params -------------------------------------------------------

params <- cfg$enrichment

# ---- Run ORA per database and direction --------------------------------------

# Collector for all results: key = "ora_{db_lower}_{direction}"
all_results <- list()

# Summary tracker: data.frame with database, direction, n_terms
summary_rows <- list()

for (db in cfg$databases) {

  db_lower <- tolower(gsub("_", "_", db))

  for (direction in c("up", "down")) {

    genes <- deg_splits[[direction]]
    name  <- sprintf("ora_%s_%s", db_lower, direction)

    message(sprintf("Running ORA: %s [%s] (%d genes) ...", db, toupper(direction),
                    length(genes)))

    res <- switch(db,

      GO_BP = run_ora_go(genes, background_genes, "BP", orgdb, params),
      GO_CC = run_ora_go(genes, background_genes, "CC", orgdb, params),
      GO_MF = run_ora_go(genes, background_genes, "MF", orgdb, params),

      KEGG = {
        ids <- if (direction == "up") entrez_up else entrez_down
        run_ora_kegg(ids, entrez_univ, kegg_org, params)
      },

      Reactome = {
        ids <- if (direction == "up") entrez_up else entrez_down
        run_ora_reactome(ids, entrez_univ, reactome_org, params)
      },

      Hallmarks    = run_ora_hallmarks(genes, background_genes, species_full, params),
      WikiPathways = run_ora_wp(genes, background_genes, species_full, params),

      {
        message(sprintf("  Unknown database '%s' -- skipping.", db))
        NULL
      }
    )

    # Count significant terms
    n_terms <- 0L
    if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
      n_terms <- nrow(as.data.frame(res))
    }

    message(sprintf("  -> %d significant terms.", n_terms))

    # Save enrichment result (RDS, CSV, dotplot)
    save_enrichment(res, tables_dir, name, cfg$plots$format)

    all_results[[name]] <- res

    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      database  = db,
      direction = toupper(direction),
      n_terms   = n_terms,
      stringsAsFactors = FALSE
    )
  }
}

# ---- Annotate artifact gene families in results ------------------------------

message("\nAnnotating enrichment results for artifact gene families ...")

# Known artifact gene families that frequently produce spurious enrichment
# in non-relevant tissues due to large family size and low/zero expression
artifact_patterns <- c(
  "^OR[0-9]"    ,  # Olfactory receptors (~400 genes)
  "^TAS[12]R"   ,  # Taste receptors
  "^KRT[0-9]"   ,  # Keratins (in non-epithelial tissues)
  "^KRTAP"         # Keratin-associated proteins
)
artifact_regex <- paste(artifact_patterns, collapse = "|")

for (name in names(all_results)) {
  res <- all_results[[name]]
  if (is.null(res)) next

  res_df <- as.data.frame(res)
  if (nrow(res_df) == 0 || !"geneID" %in% colnames(res_df)) next

  # For each term, count how many genes match artifact families
  res_df$artifact_gene_count <- sapply(res_df$geneID, function(gid) {
    genes <- unlist(strsplit(gid, "/"))
    sum(grepl(artifact_regex, genes, ignore.case = FALSE))
  })
  res_df$artifact_gene_fraction <- round(
    res_df$artifact_gene_count / sapply(res_df$geneID, function(gid) {
      length(unlist(strsplit(gid, "/")))
    }), 3
  )
  res_df$artifact_flag <- ifelse(
    res_df$artifact_gene_fraction >= 0.5,
    "ARTIFACT_DOMINATED",
    ifelse(res_df$artifact_gene_count > 0, "HAS_ARTIFACTS", "")
  )

  # Count flagged terms
  n_dominated <- sum(res_df$artifact_flag == "ARTIFACT_DOMINATED")
  n_has       <- sum(res_df$artifact_flag == "HAS_ARTIFACTS")
  if (n_dominated > 0 || n_has > 0) {
    message(sprintf("  %s: %d artifact-dominated, %d with some artifacts (of %d terms)",
                    name, n_dominated, n_has, nrow(res_df)))
  }

  # Re-save CSV with artifact annotations
  csv_path <- file.path(tables_dir, paste0(name, ".csv"))
  if (file.exists(csv_path)) {
    write.csv(res_df, csv_path, row.names = FALSE)
  }
}

# ---- Combined barplots per database ------------------------------------------

message("Generating combined UP/DOWN barplots per database ...")

colors <- default_plot_colors()

for (db in cfg$databases) {

  db_lower <- tolower(gsub("_", "_", db))

  key_up   <- sprintf("ora_%s_up",   db_lower)
  key_down <- sprintf("ora_%s_down", db_lower)

  res_up   <- all_results[[key_up]]
  res_down <- all_results[[key_down]]

  # Extract top 10 terms from each direction
  frames <- list()

  if (!is.null(res_up) && nrow(as.data.frame(res_up)) > 0) {
    df_up <- as.data.frame(res_up) %>%
      arrange(.data$p.adjust) %>%
      head(10) %>%
      mutate(Direction = "UP",
             neglog10p = -log10(pmax(.data$p.adjust, 1e-300)))
    frames[["up"]] <- df_up
  }

  if (!is.null(res_down) && nrow(as.data.frame(res_down)) > 0) {
    df_down <- as.data.frame(res_down) %>%
      arrange(.data$p.adjust) %>%
      head(10) %>%
      mutate(Direction = "DOWN",
             neglog10p = -log10(pmax(.data$p.adjust, 1e-300)))
    frames[["down"]] <- df_down
  }

  if (length(frames) == 0) {
    message(sprintf("  %s: no significant terms in either direction -- skipping plot.", db))
    next
  }

  combined_df <- bind_rows(frames)

  # Truncate long descriptions for readability
  combined_df <- combined_df %>%
    mutate(
      term_label = ifelse(nchar(.data$Description) > 60,
                          paste0(substr(.data$Description, 1, 57), "..."),
                          .data$Description),
      Direction  = factor(.data$Direction, levels = c("UP", "DOWN"))
    )

  # Order terms by significance within each direction
  combined_df <- combined_df %>%
    arrange(.data$Direction, desc(.data$neglog10p)) %>%
    mutate(term_label = factor(.data$term_label,
                               levels = rev(unique(.data$term_label))))

  p <- ggplot(combined_df, aes(x = .data$neglog10p,
                               y = .data$term_label,
                               fill = .data$Direction)) +
    geom_col(alpha = 0.85, width = 0.7) +
    facet_wrap(~ Direction, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = c("UP" = colors$up, "DOWN" = colors$down),
                      guide = "none") +
    labs(
      title = sprintf("ORA: %s — %s", db, contrast_cfg$label),
      x     = expression(-log[10]~"adjusted P"),
      y     = NULL
    ) +
    get_theme_publication(base_size = 11) +
    theme(
      axis.text.y  = element_text(size = 9),
      strip.text   = element_text(face = "bold", size = 12),
      plot.title   = element_text(size = 13)
    )

  # Determine height based on number of terms
  n_terms_total <- nrow(combined_df)
  plot_height   <- max(5, 1.5 + n_terms_total * 0.35)

  plot_base <- file.path(plots_dir, sprintf("ora_%s_combined", db_lower))
  save_plot(p, plot_base, width = 11, height = plot_height,
            dpi = cfg$plots$dpi, formats = cfg$plots$format)
}

# ---- Write flag file ---------------------------------------------------------

flag_path <- file.path(tables_dir, "ora_done.flag")
writeLines(
  sprintf("ORA completed: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  flag_path
)
message(sprintf("Flag file written: %s", flag_path))

# ---- Print summary ----------------------------------------------------------

summary_df <- bind_rows(summary_rows)

message("\n========== ORA Summary ==========")
message(sprintf("%-15s %-10s %s", "Database", "Direction", "Terms"))
message(paste(rep("-", 40), collapse = ""))

for (i in seq_len(nrow(summary_df))) {
  message(sprintf("%-15s %-10s %d",
                  summary_df$database[i],
                  summary_df$direction[i],
                  summary_df$n_terms[i]))
}

total_terms <- sum(summary_df$n_terms)
message(paste(rep("-", 40), collapse = ""))
message(sprintf("Total significant terms: %d", total_terms))
message("=================================\n")

message("03_ora_all.R completed successfully.")

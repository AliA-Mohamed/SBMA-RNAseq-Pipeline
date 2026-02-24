#!/usr/bin/env Rscript
# ══════════════════════════════════════════════════════════════════════════════
# 05_tf_ppi_disease_hpo.R — TF target, PPI, Disease, and HPO enrichment
# ══════════════════════════════════════════════════════════════════════════════
#
# Usage:
#   Rscript scripts/05_tf_ppi_disease_hpo.R \
#     --config config/config.yaml \
#     --dataset SBMA_iPSC_Q51 \
#     --contrast disease_vs_control_vehicle
#
# Description:
#   Performs transcription factor target enrichment (MSigDB C3 TFT:GTRD),
#   exports gene lists for STRINGdb PPI network analysis, runs disease
#   enrichment (DisGeNET/NCG via DOSE), and HPO term enrichment (MSigDB C5
#   HPO).  Neuro-related terms are flagged across all enrichment results.
#
# Outputs (under results/<dataset>/<contrast>/):
#   tables/tf_ora_up.csv, tf_ora_down.csv    — TF target ORA results
#   tables/tf_gsea.csv                       — TF target GSEA results (if ranked)
#   tables/ppi_up_genes.txt, ppi_down_genes.txt — Gene lists for STRINGdb
#   tables/ppi_network.csv                   — PPI interactions (if STRINGdb available)
#   tables/ppi_done.flag                     — PPI completion flag
#   tables/disease_dgn_up.csv, disease_dgn_down.csv — DisGeNET ORA results
#   tables/disease_ncg_up.csv, disease_ncg_down.csv — NCG ORA results
#   tables/hpo_ora_up.csv, hpo_ora_down.csv  — HPO ORA results
#   tables/neuro_flagged_terms.csv           — Neuro-related terms across analyses
#   tables/tf_done.flag                      — Pipeline step completion flag
#   plots/ora/tf_*, plots/ora/disease_*, plots/ora/hpo_* — Dotplots
# ══════════════════════════════════════════════════════════════════════════════

# ── Load libraries ───────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
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
cat("  SBMA RNAseq Pipeline — Step 05: TF, PPI, Disease & HPO\n")
cat("════════════════════════════════════════════════════════════════════\n")
cat(sprintf("  Config   : %s\n", args$config))
cat(sprintf("  Dataset  : %s\n", args$dataset))
cat(sprintf("  Contrast : %s\n", args$contrast))
cat("────────────────────────────────────────────────────────────────────\n\n")

# ── 1. Load configuration and contrast info ──────────────────────────────────
cat("[1/19] Loading configuration...\n")
cfg          <- load_config(args$config)
contrast_cfg <- get_contrast(cfg, args$dataset, args$contrast)
out_dir      <- get_output_dir(cfg, args$dataset, args$contrast)
params       <- cfg$enrichment
plot_formats <- cfg$plots$format %||% c("pdf", "png")

cat(sprintf("       Label       : %s\n", contrast_cfg$label))
cat(sprintf("       Full ranked : %s\n", contrast_cfg$full_ranked))
cat(sprintf("       Output      : %s\n", out_dir))

# ── 2. Set up output directories ─────────────────────────────────────────────
cat("\n[2/19] Creating output directories...\n")
tables_dir <- file.path(out_dir, "tables")
data_dir   <- file.path(out_dir, "data")
ora_dir    <- file.path(out_dir, "plots", "ora")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ora_dir,    recursive = TRUE, showWarnings = FALSE)

# ── 3. Load prepared DEG data ────────────────────────────────────────────────
cat("\n[3/19] Loading DEG data...\n")

deg_splits_path <- file.path(data_dir, "deg_splits.rds")
deg_all_path    <- file.path(data_dir, "deg_all.rds")

if (!file.exists(deg_splits_path) || !file.exists(deg_all_path)) {
  stop(sprintf("Required input files not found in %s. Run step 01 first.", data_dir),
       call. = FALSE)
}

deg_splits <- readRDS(deg_splits_path)
deg_all    <- readRDS(deg_all_path)

cat(sprintf("       Up genes         : %d\n", length(deg_splits$up)))
cat(sprintf("       Down genes       : %d\n", length(deg_splits$down)))
cat(sprintf("       Background genes : %d\n", length(deg_splits$background)))

# Load ranked lists if available
ranked_symbol <- NULL
ranked_symbol_path <- file.path(data_dir, "ranked_symbol.rds")
if (isTRUE(contrast_cfg$full_ranked) && file.exists(ranked_symbol_path)) {
  ranked_symbol <- readRDS(ranked_symbol_path)
  cat(sprintf("       Ranked genes     : %d\n", length(ranked_symbol)))
}

# ── Resolve species for msigdbr ──────────────────────────────────────────────
species_full <- switch(tolower(cfg$species),
  human = "Homo sapiens",
  mouse = "Mus musculus",
  stop(sprintf("Unsupported species for msigdbr: '%s'", cfg$species), call. = FALSE)
)
cat(sprintf("       Species (msigdbr): %s\n", species_full))

# Resolve NCBI taxonomy ID for STRINGdb
string_species <- switch(tolower(cfg$species),
  human = 9606L,
  mouse = 10090L,
  9606L
)

# ── Collector for all enrichment results (for neuro flagging) ────────────────
all_enrichment_results <- list()

# =============================================================================
# TF TARGET ENRICHMENT (MSigDB C3 TFT:GTRD)
# =============================================================================

cat("\n[4/19] Loading TF target gene sets from MSigDB C3 (TFT:GTRD)...\n")

tf_t2g <- tryCatch({
  tf_sets <- msigdbr(species = species_full, category = "C3", subcategory = "TFT:GTRD")
  tf_t2g  <- tf_sets[, c("gs_name", "gene_symbol")]
  cat(sprintf("       TF gene sets : %d\n", length(unique(tf_t2g$gs_name))))
  cat(sprintf("       Total mappings: %d\n", nrow(tf_t2g)))
  tf_t2g
}, error = function(e) {
  message(sprintf("       WARNING: Failed to load TF targets — %s", conditionMessage(e)))
  NULL
})

# ── 5. TF ORA — Upregulated genes ───────────────────────────────────────────
cat("\n[5/19] Running TF target ORA on upregulated genes...\n")

tf_ora_up <- NULL
if (!is.null(tf_t2g) && length(deg_splits$up) >= 5L) {
  tf_ora_up <- tryCatch({
    res <- enricher(
      gene          = deg_splits$up,
      universe      = deg_splits$background,
      TERM2GENE     = tf_t2g,
      pAdjustMethod = "BH",
      pvalueCutoff  = params$pvalue_cutoff,
      qvalueCutoff  = params$qvalue_cutoff,
      minGSSize     = params$min_gs_size,
      maxGSSize     = params$max_gs_size
    )
    if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
      cat(sprintf("       Significant TF sets (up): %d\n", nrow(as.data.frame(res))))
      res
    } else {
      cat("       No significant TF enrichments for upregulated genes.\n")
      NULL
    }
  }, error = function(e) {
    message(sprintf("       WARNING: TF ORA (up) failed — %s", conditionMessage(e)))
    NULL
  })
} else {
  cat("       Skipping: insufficient upregulated genes or TF data unavailable.\n")
}

# ── 6. TF ORA — Downregulated genes ─────────────────────────────────────────
cat("\n[6/19] Running TF target ORA on downregulated genes...\n")

tf_ora_down <- NULL
if (!is.null(tf_t2g) && length(deg_splits$down) >= 5L) {
  tf_ora_down <- tryCatch({
    res <- enricher(
      gene          = deg_splits$down,
      universe      = deg_splits$background,
      TERM2GENE     = tf_t2g,
      pAdjustMethod = "BH",
      pvalueCutoff  = params$pvalue_cutoff,
      qvalueCutoff  = params$qvalue_cutoff,
      minGSSize     = params$min_gs_size,
      maxGSSize     = params$max_gs_size
    )
    if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
      cat(sprintf("       Significant TF sets (down): %d\n", nrow(as.data.frame(res))))
      res
    } else {
      cat("       No significant TF enrichments for downregulated genes.\n")
      NULL
    }
  }, error = function(e) {
    message(sprintf("       WARNING: TF ORA (down) failed — %s", conditionMessage(e)))
    NULL
  })
} else {
  cat("       Skipping: insufficient downregulated genes or TF data unavailable.\n")
}

# ── 7. TF GSEA (if full_ranked) ─────────────────────────────────────────────
cat("\n[7/19] Running TF target GSEA...\n")

tf_gsea <- NULL
if (!is.null(tf_t2g) && !is.null(ranked_symbol) && length(ranked_symbol) >= 10L) {
  tf_gsea <- tryCatch({
    set.seed(params$seed)
    res <- GSEA(
      geneList      = ranked_symbol,
      TERM2GENE     = tf_t2g,
      pAdjustMethod = "BH",
      pvalueCutoff  = params$pvalue_cutoff,
      minGSSize     = params$min_gs_size,
      maxGSSize     = params$max_gs_size,
      nPermSimple   = params$n_perm
    )
    if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
      cat(sprintf("       Significant TF GSEA terms: %d\n", nrow(as.data.frame(res))))
      res
    } else {
      cat("       No significant TF GSEA enrichments.\n")
      NULL
    }
  }, error = function(e) {
    message(sprintf("       WARNING: TF GSEA failed — %s", conditionMessage(e)))
    NULL
  })
} else {
  cat("       Skipping: ranked list not available or TF data unavailable.\n")
}

# ── 8. Save TF results ──────────────────────────────────────────────────────
cat("\n[8/19] Saving TF target enrichment results...\n")

save_enrichment(tf_ora_up,   tables_dir, "tf_ora_up",   formats = plot_formats)
save_enrichment(tf_ora_down, tables_dir, "tf_ora_down", formats = plot_formats)
save_enrichment(tf_gsea,     tables_dir, "tf_gsea",     formats = plot_formats)

# Save dotplots to plots/ora/ directory
if (!is.null(tf_ora_up) && nrow(as.data.frame(tf_ora_up)) > 0) {
  tryCatch({
    n_show <- min(20, nrow(as.data.frame(tf_ora_up)))
    p <- dotplot(tf_ora_up, showCategory = n_show) +
      ggtitle("TF Target Enrichment — Upregulated") +
      theme_minimal(base_size = 12)
    for (fmt in plot_formats) {
      ggsave(file.path(ora_dir, paste0("tf_ora_up_dotplot.", fmt)),
             plot = p, width = 12, height = 8, dpi = 300)
    }
    cat("       Saved TF ORA (up) dotplot.\n")
  }, error = function(e) {
    message(sprintf("       WARNING: TF ORA (up) dotplot failed — %s", conditionMessage(e)))
  })
  all_enrichment_results[["tf_ora_up"]] <- as.data.frame(tf_ora_up)
}

if (!is.null(tf_ora_down) && nrow(as.data.frame(tf_ora_down)) > 0) {
  tryCatch({
    n_show <- min(20, nrow(as.data.frame(tf_ora_down)))
    p <- dotplot(tf_ora_down, showCategory = n_show) +
      ggtitle("TF Target Enrichment — Downregulated") +
      theme_minimal(base_size = 12)
    for (fmt in plot_formats) {
      ggsave(file.path(ora_dir, paste0("tf_ora_down_dotplot.", fmt)),
             plot = p, width = 12, height = 8, dpi = 300)
    }
    cat("       Saved TF ORA (down) dotplot.\n")
  }, error = function(e) {
    message(sprintf("       WARNING: TF ORA (down) dotplot failed — %s", conditionMessage(e)))
  })
  all_enrichment_results[["tf_ora_down"]] <- as.data.frame(tf_ora_down)
}

if (!is.null(tf_gsea) && nrow(as.data.frame(tf_gsea)) > 0) {
  all_enrichment_results[["tf_gsea"]] <- as.data.frame(tf_gsea)
}

# =============================================================================
# PPI NETWORK EXPORT
# =============================================================================

cat("\n[9/19] Exporting gene lists for PPI analysis...\n")

# Write gene lists (one gene per line)
up_genes_path   <- file.path(tables_dir, "ppi_up_genes.txt")
down_genes_path <- file.path(tables_dir, "ppi_down_genes.txt")

writeLines(deg_splits$up,   up_genes_path)
writeLines(deg_splits$down, down_genes_path)
cat(sprintf("       Saved: %s (%d genes)\n", up_genes_path, length(deg_splits$up)))
cat(sprintf("       Saved: %s (%d genes)\n", down_genes_path, length(deg_splits$down)))

# ── 10. Attempt automated STRINGdb analysis ──────────────────────────────────
cat("\n[10/19] Attempting STRINGdb PPI network construction...\n")

ppi_network <- NULL
if (requireNamespace("STRINGdb", quietly = TRUE)) {
  tryCatch({
    string_db <- STRINGdb::STRINGdb$new(
      version         = "11.5",
      species         = string_species,
      score_threshold = 400
    )

    # Prepare gene data frame for mapping
    all_sig_genes <- unique(c(deg_splits$up, deg_splits$down))
    gene_df <- data.frame(gene = all_sig_genes, stringsAsFactors = FALSE)

    # Map genes to STRING IDs
    mapped <- string_db$map(gene_df, "gene", removeUnmappedRows = TRUE)
    cat(sprintf("       Mapped to STRING: %d / %d genes\n",
                nrow(mapped), length(all_sig_genes)))

    if (nrow(mapped) > 1) {
      # Get interactions
      interactions <- string_db$get_interactions(mapped$STRING_id)

      if (nrow(interactions) > 0) {
        ppi_network <- interactions
        write.csv(ppi_network,
                  file.path(tables_dir, "ppi_network.csv"),
                  row.names = FALSE)
        cat(sprintf("       PPI interactions: %d\n", nrow(ppi_network)))
        cat(sprintf("       Saved: %s\n", file.path(tables_dir, "ppi_network.csv")))

        # Attempt network plot
        tryCatch({
          png_path <- file.path(ora_dir, "ppi_network.png")
          string_db$plot_network(mapped$STRING_id, payload_id = NULL)
          cat("       STRINGdb network plot generated.\n")
        }, error = function(e) {
          message(sprintf("       WARNING: PPI network plot failed — %s", conditionMessage(e)))
        })
      } else {
        cat("       No PPI interactions found above threshold.\n")
      }
    } else {
      cat("       Too few genes mapped to STRING for network analysis.\n")
    }
  }, error = function(e) {
    message(sprintf("       WARNING: STRINGdb analysis failed — %s", conditionMessage(e)))
  })
} else {
  cat("       STRINGdb package not available. Gene lists exported for manual use.\n")
  cat("       Install with: BiocManager::install('STRINGdb')\n")
}

# Write PPI completion flag
writeLines(
  sprintf("PPI export completed: %s", Sys.time()),
  file.path(tables_dir, "ppi_done.flag")
)
cat(sprintf("       Flag: %s\n", file.path(tables_dir, "ppi_done.flag")))

# =============================================================================
# DISEASE ENRICHMENT (DOSE)
# =============================================================================

cat("\n[11/19] Running disease enrichment analysis...\n")

# Convert significant genes to Entrez IDs for DOSE functions
entrez_up   <- tryCatch(
  symbols_to_entrez(deg_splits$up, cfg$orgdb),
  error = function(e) { message("       Entrez mapping (up) failed: ", conditionMessage(e)); character(0) }
)
entrez_down <- tryCatch(
  symbols_to_entrez(deg_splits$down, cfg$orgdb),
  error = function(e) { message("       Entrez mapping (down) failed: ", conditionMessage(e)); character(0) }
)

cat(sprintf("       Entrez IDs (up)   : %d\n", length(entrez_up)))
cat(sprintf("       Entrez IDs (down) : %d\n", length(entrez_down)))

disease_dgn_up   <- NULL
disease_dgn_down <- NULL
disease_ncg_up   <- NULL
disease_ncg_down <- NULL

if (requireNamespace("DOSE", quietly = TRUE)) {

  # ── DisGeNET enrichment ────────────────────────────────────────────────────
  cat("\n[12/19] Running DisGeNET enrichment...\n")

  if (length(entrez_up) >= 5L) {
    disease_dgn_up <- tryCatch({
      res <- DOSE::enrichDGN(
        gene          = unname(entrez_up),
        pAdjustMethod = "BH",
        pvalueCutoff  = params$pvalue_cutoff,
        qvalueCutoff  = params$qvalue_cutoff,
        minGSSize     = params$min_gs_size,
        maxGSSize     = params$max_gs_size,
        readable      = TRUE
      )
      if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
        cat(sprintf("       DisGeNET (up): %d terms\n", nrow(as.data.frame(res))))
        res
      } else {
        cat("       DisGeNET (up): no significant terms.\n")
        NULL
      }
    }, error = function(e) {
      message(sprintf("       WARNING: DisGeNET (up) failed — %s", conditionMessage(e)))
      NULL
    })
  } else {
    cat("       Skipping DisGeNET (up): insufficient Entrez IDs.\n")
  }

  if (length(entrez_down) >= 5L) {
    disease_dgn_down <- tryCatch({
      res <- DOSE::enrichDGN(
        gene          = unname(entrez_down),
        pAdjustMethod = "BH",
        pvalueCutoff  = params$pvalue_cutoff,
        qvalueCutoff  = params$qvalue_cutoff,
        minGSSize     = params$min_gs_size,
        maxGSSize     = params$max_gs_size,
        readable      = TRUE
      )
      if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
        cat(sprintf("       DisGeNET (down): %d terms\n", nrow(as.data.frame(res))))
        res
      } else {
        cat("       DisGeNET (down): no significant terms.\n")
        NULL
      }
    }, error = function(e) {
      message(sprintf("       WARNING: DisGeNET (down) failed — %s", conditionMessage(e)))
      NULL
    })
  } else {
    cat("       Skipping DisGeNET (down): insufficient Entrez IDs.\n")
  }

  # ── NCG enrichment ────────────────────────────────────────────────────────
  cat("\n[13/19] Running NCG (Network of Cancer Genes) enrichment...\n")

  if (length(entrez_up) >= 5L) {
    disease_ncg_up <- tryCatch({
      res <- DOSE::enrichNCG(
        gene          = unname(entrez_up),
        pAdjustMethod = "BH",
        pvalueCutoff  = params$pvalue_cutoff,
        qvalueCutoff  = params$qvalue_cutoff,
        minGSSize     = params$min_gs_size,
        maxGSSize     = params$max_gs_size,
        readable      = TRUE
      )
      if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
        cat(sprintf("       NCG (up): %d terms\n", nrow(as.data.frame(res))))
        res
      } else {
        cat("       NCG (up): no significant terms.\n")
        NULL
      }
    }, error = function(e) {
      message(sprintf("       WARNING: NCG (up) failed — %s", conditionMessage(e)))
      NULL
    })
  } else {
    cat("       Skipping NCG (up): insufficient Entrez IDs.\n")
  }

  if (length(entrez_down) >= 5L) {
    disease_ncg_down <- tryCatch({
      res <- DOSE::enrichNCG(
        gene          = unname(entrez_down),
        pAdjustMethod = "BH",
        pvalueCutoff  = params$pvalue_cutoff,
        qvalueCutoff  = params$qvalue_cutoff,
        minGSSize     = params$min_gs_size,
        maxGSSize     = params$max_gs_size,
        readable      = TRUE
      )
      if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
        cat(sprintf("       NCG (down): %d terms\n", nrow(as.data.frame(res))))
        res
      } else {
        cat("       NCG (down): no significant terms.\n")
        NULL
      }
    }, error = function(e) {
      message(sprintf("       WARNING: NCG (down) failed — %s", conditionMessage(e)))
      NULL
    })
  } else {
    cat("       Skipping NCG (down): insufficient Entrez IDs.\n")
  }

} else {
  cat("       DOSE package not available. Skipping disease enrichment.\n")
  cat("       Install with: BiocManager::install('DOSE')\n")
  cat("\n[12/19] Skipping DisGeNET...\n")
  cat("\n[13/19] Skipping NCG...\n")
}

# ── 14. Save disease enrichment results ──────────────────────────────────────
cat("\n[14/19] Saving disease enrichment results...\n")

save_enrichment(disease_dgn_up,   tables_dir, "disease_dgn_up",   formats = plot_formats)
save_enrichment(disease_dgn_down, tables_dir, "disease_dgn_down", formats = plot_formats)
save_enrichment(disease_ncg_up,   tables_dir, "disease_ncg_up",   formats = plot_formats)
save_enrichment(disease_ncg_down, tables_dir, "disease_ncg_down", formats = plot_formats)

# Save dotplots to plots/ora/ directory
for (res_info in list(
  list(res = disease_dgn_up,   name = "disease_dgn_up",   title = "DisGeNET — Upregulated"),
  list(res = disease_dgn_down, name = "disease_dgn_down", title = "DisGeNET — Downregulated"),
  list(res = disease_ncg_up,   name = "disease_ncg_up",   title = "NCG — Upregulated"),
  list(res = disease_ncg_down, name = "disease_ncg_down", title = "NCG — Downregulated")
)) {
  if (!is.null(res_info$res) && nrow(as.data.frame(res_info$res)) > 0) {
    tryCatch({
      n_show <- min(20, nrow(as.data.frame(res_info$res)))
      p <- dotplot(res_info$res, showCategory = n_show) +
        ggtitle(res_info$title) +
        theme_minimal(base_size = 12)
      for (fmt in plot_formats) {
        ggsave(file.path(ora_dir, paste0(res_info$name, "_dotplot.", fmt)),
               plot = p, width = 12, height = 8, dpi = 300)
      }
      cat(sprintf("       Saved %s dotplot.\n", res_info$name))
    }, error = function(e) {
      message(sprintf("       WARNING: %s dotplot failed — %s", res_info$name, conditionMessage(e)))
    })
    all_enrichment_results[[res_info$name]] <- as.data.frame(res_info$res)
  }
}

# =============================================================================
# HPO ENRICHMENT (MSigDB C5 HPO)
# =============================================================================

cat("\n[15/19] Loading HPO gene sets from MSigDB C5 (HPO)...\n")

hpo_t2g <- tryCatch({
  hpo_sets <- msigdbr(species = species_full, category = "C5", subcategory = "HPO")
  hpo_t2g  <- hpo_sets[, c("gs_name", "gene_symbol")]
  cat(sprintf("       HPO gene sets : %d\n", length(unique(hpo_t2g$gs_name))))
  cat(sprintf("       Total mappings: %d\n", nrow(hpo_t2g)))
  hpo_t2g
}, error = function(e) {
  message(sprintf("       WARNING: Failed to load HPO terms — %s", conditionMessage(e)))
  NULL
})

# ── 16. HPO ORA — Upregulated genes ─────────────────────────────────────────
cat("\n[16/19] Running HPO ORA on upregulated genes...\n")

hpo_ora_up <- NULL
if (!is.null(hpo_t2g) && length(deg_splits$up) >= 5L) {
  hpo_ora_up <- tryCatch({
    res <- enricher(
      gene          = deg_splits$up,
      universe      = deg_splits$background,
      TERM2GENE     = hpo_t2g,
      pAdjustMethod = "BH",
      pvalueCutoff  = params$pvalue_cutoff,
      qvalueCutoff  = params$qvalue_cutoff,
      minGSSize     = params$min_gs_size,
      maxGSSize     = params$max_gs_size
    )
    if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
      cat(sprintf("       Significant HPO terms (up): %d\n", nrow(as.data.frame(res))))
      res
    } else {
      cat("       No significant HPO enrichments for upregulated genes.\n")
      NULL
    }
  }, error = function(e) {
    message(sprintf("       WARNING: HPO ORA (up) failed — %s", conditionMessage(e)))
    NULL
  })
} else {
  cat("       Skipping: insufficient upregulated genes or HPO data unavailable.\n")
}

# ── HPO ORA — Downregulated genes ───────────────────────────────────────────
cat("\n[16b/19] Running HPO ORA on downregulated genes...\n")

hpo_ora_down <- NULL
if (!is.null(hpo_t2g) && length(deg_splits$down) >= 5L) {
  hpo_ora_down <- tryCatch({
    res <- enricher(
      gene          = deg_splits$down,
      universe      = deg_splits$background,
      TERM2GENE     = hpo_t2g,
      pAdjustMethod = "BH",
      pvalueCutoff  = params$pvalue_cutoff,
      qvalueCutoff  = params$qvalue_cutoff,
      minGSSize     = params$min_gs_size,
      maxGSSize     = params$max_gs_size
    )
    if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
      cat(sprintf("       Significant HPO terms (down): %d\n", nrow(as.data.frame(res))))
      res
    } else {
      cat("       No significant HPO enrichments for downregulated genes.\n")
      NULL
    }
  }, error = function(e) {
    message(sprintf("       WARNING: HPO ORA (down) failed — %s", conditionMessage(e)))
    NULL
  })
} else {
  cat("       Skipping: insufficient downregulated genes or HPO data unavailable.\n")
}

# ── 17. Save HPO results ────────────────────────────────────────────────────
cat("\n[17/19] Saving HPO enrichment results...\n")

save_enrichment(hpo_ora_up,   tables_dir, "hpo_ora_up",   formats = plot_formats)
save_enrichment(hpo_ora_down, tables_dir, "hpo_ora_down", formats = plot_formats)

for (res_info in list(
  list(res = hpo_ora_up,   name = "hpo_ora_up",   title = "HPO Enrichment — Upregulated"),
  list(res = hpo_ora_down, name = "hpo_ora_down", title = "HPO Enrichment — Downregulated")
)) {
  if (!is.null(res_info$res) && nrow(as.data.frame(res_info$res)) > 0) {
    tryCatch({
      n_show <- min(20, nrow(as.data.frame(res_info$res)))
      p <- dotplot(res_info$res, showCategory = n_show) +
        ggtitle(res_info$title) +
        theme_minimal(base_size = 12)
      for (fmt in plot_formats) {
        ggsave(file.path(ora_dir, paste0(res_info$name, "_dotplot.", fmt)),
               plot = p, width = 12, height = 8, dpi = 300)
      }
      cat(sprintf("       Saved %s dotplot.\n", res_info$name))
    }, error = function(e) {
      message(sprintf("       WARNING: %s dotplot failed — %s", res_info$name, conditionMessage(e)))
    })
    all_enrichment_results[[res_info$name]] <- as.data.frame(res_info$res)
  }
}

# =============================================================================
# NEURO-FLAGGED SUBSET
# =============================================================================

cat("\n[18/19] Scanning enrichment results for neuro-related terms...\n")

neuro_pattern <- "neuro|nerve|motor|muscle|axon|synap|dendrit"

neuro_flagged <- tryCatch({
  flagged_rows <- list()

  for (result_name in names(all_enrichment_results)) {
    df <- all_enrichment_results[[result_name]]

    if (is.null(df) || nrow(df) == 0) next

    # Search in Description or ID columns for neuro-related terms
    desc_col <- if ("Description" %in% colnames(df)) "Description" else if ("ID" %in% colnames(df)) "ID" else NULL

    if (!is.null(desc_col)) {
      neuro_idx <- grep(neuro_pattern, df[[desc_col]], ignore.case = TRUE)
      if (length(neuro_idx) > 0) {
        subset_df <- df[neuro_idx, , drop = FALSE]
        subset_df$source_analysis <- result_name
        subset_df$matched_column  <- desc_col
        flagged_rows[[length(flagged_rows) + 1]] <- subset_df
      }
    }

    # Also search in ID if different from Description
    if ("ID" %in% colnames(df) && !is.null(desc_col) && desc_col != "ID") {
      id_idx <- grep(neuro_pattern, df[["ID"]], ignore.case = TRUE)
      # Exclude rows already captured via Description
      if (!is.null(desc_col)) {
        id_idx <- setdiff(id_idx, grep(neuro_pattern, df[[desc_col]], ignore.case = TRUE))
      }
      if (length(id_idx) > 0) {
        subset_df <- df[id_idx, , drop = FALSE]
        subset_df$source_analysis <- result_name
        subset_df$matched_column  <- "ID"
        flagged_rows[[length(flagged_rows) + 1]] <- subset_df
      }
    }
  }

  if (length(flagged_rows) > 0) {
    neuro_df <- tryCatch(
      dplyr::bind_rows(flagged_rows),
      error = function(e) {
        # Fallback: select common columns only
        common_cols <- Reduce(intersect, lapply(flagged_rows, colnames))
        do.call(rbind, lapply(flagged_rows, function(x) x[, common_cols, drop = FALSE]))
      }
    )
    neuro_df
  } else {
    NULL
  }
}, error = function(e) {
  message(sprintf("       WARNING: Neuro flagging failed — %s", conditionMessage(e)))
  NULL
})

if (!is.null(neuro_flagged) && nrow(neuro_flagged) > 0) {
  neuro_path <- file.path(tables_dir, "neuro_flagged_terms.csv")
  write.csv(neuro_flagged, neuro_path, row.names = FALSE)
  cat(sprintf("       Neuro-flagged terms: %d\n", nrow(neuro_flagged)))
  cat(sprintf("       Saved: %s\n", neuro_path))
} else {
  cat("       No neuro-related terms found across enrichment results.\n")
  # Write empty file for pipeline consistency
  write.csv(
    data.frame(note = "No neuro-related terms identified", stringsAsFactors = FALSE),
    file.path(tables_dir, "neuro_flagged_terms.csv"),
    row.names = FALSE
  )
}

# =============================================================================
# FLAG FILE AND SUMMARY
# =============================================================================

cat("\n[19/19] Writing completion flag...\n")
writeLines(
  sprintf("Step 05 completed: %s", Sys.time()),
  file.path(tables_dir, "tf_done.flag")
)
cat(sprintf("       Flag: %s\n", file.path(tables_dir, "tf_done.flag")))

# ── Summary ──────────────────────────────────────────────────────────────────
cat("\n════════════════════════════════════════════════════════════════════\n")
cat("  Step 05 complete — TF, PPI, Disease & HPO Enrichment\n")
cat("════════════════════════════════════════════════════════════════════\n")
cat(sprintf("  Dataset  : %s\n", args$dataset))
cat(sprintf("  Contrast : %s (%s)\n", args$contrast, contrast_cfg$label))
cat("────────────────────────────────────────────────────────────────────\n")
cat("  TF Target Enrichment:\n")
cat(sprintf("    ORA up       : %s\n",
            if (!is.null(tf_ora_up)) paste(nrow(as.data.frame(tf_ora_up)), "terms") else "skipped/none"))
cat(sprintf("    ORA down     : %s\n",
            if (!is.null(tf_ora_down)) paste(nrow(as.data.frame(tf_ora_down)), "terms") else "skipped/none"))
cat(sprintf("    GSEA         : %s\n",
            if (!is.null(tf_gsea)) paste(nrow(as.data.frame(tf_gsea)), "terms") else "skipped/none"))
cat("  PPI Network:\n")
cat(sprintf("    Up genes     : %d exported\n", length(deg_splits$up)))
cat(sprintf("    Down genes   : %d exported\n", length(deg_splits$down)))
cat(sprintf("    Interactions : %s\n",
            if (!is.null(ppi_network)) paste(nrow(ppi_network), "edges") else "STRINGdb not available"))
cat("  Disease Enrichment:\n")
cat(sprintf("    DisGeNET up  : %s\n",
            if (!is.null(disease_dgn_up)) paste(nrow(as.data.frame(disease_dgn_up)), "terms") else "skipped/none"))
cat(sprintf("    DisGeNET down: %s\n",
            if (!is.null(disease_dgn_down)) paste(nrow(as.data.frame(disease_dgn_down)), "terms") else "skipped/none"))
cat(sprintf("    NCG up       : %s\n",
            if (!is.null(disease_ncg_up)) paste(nrow(as.data.frame(disease_ncg_up)), "terms") else "skipped/none"))
cat(sprintf("    NCG down     : %s\n",
            if (!is.null(disease_ncg_down)) paste(nrow(as.data.frame(disease_ncg_down)), "terms") else "skipped/none"))
cat("  HPO Enrichment:\n")
cat(sprintf("    ORA up       : %s\n",
            if (!is.null(hpo_ora_up)) paste(nrow(as.data.frame(hpo_ora_up)), "terms") else "skipped/none"))
cat(sprintf("    ORA down     : %s\n",
            if (!is.null(hpo_ora_down)) paste(nrow(as.data.frame(hpo_ora_down)), "terms") else "skipped/none"))
cat("  Neuro Flagging:\n")
cat(sprintf("    Flagged terms: %s\n",
            if (!is.null(neuro_flagged) && nrow(neuro_flagged) > 0) nrow(neuro_flagged) else "none"))
cat(sprintf("  Output   : %s\n", out_dir))
cat("════════════════════════════════════════════════════════════════════\n")

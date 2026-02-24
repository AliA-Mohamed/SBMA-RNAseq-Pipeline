# enrichment_utils.R
# ORA and GSEA wrappers for all supported databases in the SBMA RNAseq pipeline.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(ReactomePA)
  library(msigdbr)
  library(enrichplot)
  library(ggplot2)
})

# ---- Constants --------------------------------------------------------------

.ORA_MIN_GENES  <- 5L
.GSEA_MIN_GENES <- 10L


# =============================================================================
# ORA wrappers
# =============================================================================

#' Run ORA with Gene Ontology (enrichGO)
#'
#' @param genes     Character vector of gene symbols.
#' @param universe  Character vector of background gene symbols.
#' @param ont       GO ontology: "BP", "CC", or "MF".
#' @param orgdb     OrgDb annotation package name (e.g. "org.Hs.eg.db").
#' @param params    List with pvalue_cutoff, qvalue_cutoff, simplify_cutoff.
#' @return An enrichResult object, or NULL on failure / insufficient genes.
#' @export
run_ora_go <- function(genes, universe, ont, orgdb, params) {

  if (length(genes) < .ORA_MIN_GENES) {
    message(sprintf("run_ora_go [%s]: fewer than %d genes — skipping.",
                    ont, .ORA_MIN_GENES))
    return(NULL)
  }

  tryCatch({
    res <- enrichGO(
      gene         = genes,
      universe     = universe,
      OrgDb        = orgdb,
      ont          = ont,
      keyType      = "SYMBOL",
      pAdjustMethod = "BH",
      pvalueCutoff = params$pvalue_cutoff,
      qvalueCutoff = params$qvalue_cutoff,
      minGSSize    = params$min_gs_size,
      maxGSSize    = params$max_gs_size,
      readable     = TRUE
    )

    if (!is.null(res) && nrow(res) > 0 && ont == "BP") {
      res <- clusterProfiler::simplify(res, cutoff = params$simplify_cutoff)
    }

    res
  }, error = function(e) {
    message(sprintf("run_ora_go [%s]: %s", ont, conditionMessage(e)))
    NULL
  })
}


#' Run ORA with KEGG (enrichKEGG)
#'
#' @param entrez_ids      Character vector of Entrez gene IDs.
#' @param universe_entrez Character vector of background Entrez IDs.
#' @param organism        KEGG organism code (e.g. "hsa", "mmu").
#' @param params          List with pvalue_cutoff, qvalue_cutoff, etc.
#' @return An enrichResult object, or NULL on failure / insufficient genes.
#' @export
run_ora_kegg <- function(entrez_ids, universe_entrez, organism, params) {

  if (length(entrez_ids) < .ORA_MIN_GENES) {
    message(sprintf("run_ora_kegg: fewer than %d Entrez IDs — skipping.",
                    .ORA_MIN_GENES))
    return(NULL)
  }

  tryCatch({
    enrichKEGG(
      gene         = entrez_ids,
      universe     = universe_entrez,
      organism     = organism,
      pAdjustMethod = "BH",
      pvalueCutoff = params$pvalue_cutoff,
      qvalueCutoff = params$qvalue_cutoff,
      minGSSize    = params$min_gs_size,
      maxGSSize    = params$max_gs_size
    )
  }, error = function(e) {
    message(sprintf("run_ora_kegg: %s", conditionMessage(e)))
    NULL
  })
}


#' Run ORA with Reactome (enrichPathway)
#'
#' @param entrez_ids      Character vector of Entrez gene IDs.
#' @param universe_entrez Character vector of background Entrez IDs.
#' @param organism        Reactome organism label ("human" or "mouse").
#' @param params          List with pvalue_cutoff, qvalue_cutoff, etc.
#' @return An enrichResult object, or NULL on failure / insufficient genes.
#' @export
run_ora_reactome <- function(entrez_ids, universe_entrez, organism, params) {

  if (length(entrez_ids) < .ORA_MIN_GENES) {
    message(sprintf("run_ora_reactome: fewer than %d Entrez IDs — skipping.",
                    .ORA_MIN_GENES))
    return(NULL)
  }

  tryCatch({
    ReactomePA::enrichPathway(
      gene         = entrez_ids,
      universe     = universe_entrez,
      organism     = organism,
      pAdjustMethod = "BH",
      pvalueCutoff = params$pvalue_cutoff,
      qvalueCutoff = params$qvalue_cutoff,
      minGSSize    = params$min_gs_size,
      maxGSSize    = params$max_gs_size,
      readable     = TRUE
    )
  }, error = function(e) {
    message(sprintf("run_ora_reactome: %s", conditionMessage(e)))
    NULL
  })
}


#' Run ORA with MSigDB Hallmark gene sets
#'
#' @param genes    Character vector of gene symbols.
#' @param universe Character vector of background gene symbols.
#' @param species  msigdbr species name ("Homo sapiens" or "Mus musculus").
#' @param params   List with pvalue_cutoff, qvalue_cutoff, etc.
#' @return An enrichResult object, or NULL on failure / insufficient genes.
#' @export
run_ora_hallmarks <- function(genes, universe, species, params) {

  if (length(genes) < .ORA_MIN_GENES) {
    message(sprintf("run_ora_hallmarks: fewer than %d genes — skipping.",
                    .ORA_MIN_GENES))
    return(NULL)
  }

  tryCatch({
    hallmark_sets <- msigdbr(species = species, category = "H")
    hallmark_t2g  <- hallmark_sets[, c("gs_name", "gene_symbol")]

    enricher(
      gene         = genes,
      universe     = universe,
      TERM2GENE    = hallmark_t2g,
      pAdjustMethod = "BH",
      pvalueCutoff = params$pvalue_cutoff,
      qvalueCutoff = params$qvalue_cutoff,
      minGSSize    = params$min_gs_size,
      maxGSSize    = params$max_gs_size
    )
  }, error = function(e) {
    message(sprintf("run_ora_hallmarks: %s", conditionMessage(e)))
    NULL
  })
}


#' Run ORA with WikiPathways gene sets (via msigdbr)
#'
#' @param genes    Character vector of gene symbols.
#' @param universe Character vector of background gene symbols.
#' @param species  msigdbr species name ("Homo sapiens" or "Mus musculus").
#' @param params   List with pvalue_cutoff, qvalue_cutoff, etc.
#' @return An enrichResult object, or NULL on failure / insufficient genes.
#' @export
run_ora_wp <- function(genes, universe, species, params) {

  if (length(genes) < .ORA_MIN_GENES) {
    message(sprintf("run_ora_wp: fewer than %d genes — skipping.",
                    .ORA_MIN_GENES))
    return(NULL)
  }

  tryCatch({
    wp_sets <- msigdbr(species = species, category = "C2",
                       subcategory = "CP:WIKIPATHWAYS")
    wp_t2g  <- wp_sets[, c("gs_name", "gene_symbol")]

    enricher(
      gene         = genes,
      universe     = universe,
      TERM2GENE    = wp_t2g,
      pAdjustMethod = "BH",
      pvalueCutoff = params$pvalue_cutoff,
      qvalueCutoff = params$qvalue_cutoff,
      minGSSize    = params$min_gs_size,
      maxGSSize    = params$max_gs_size
    )
  }, error = function(e) {
    message(sprintf("run_ora_wp: %s", conditionMessage(e)))
    NULL
  })
}


# =============================================================================
# GSEA wrappers
# =============================================================================

#' Run GSEA with Gene Ontology (gseGO)
#'
#' @param ranked Named numeric vector (gene symbols -> ranking metric),
#'               sorted in decreasing order.
#' @param ont    GO ontology: "BP", "CC", or "MF".
#' @param orgdb  OrgDb annotation package name.
#' @param params List with n_perm, seed, min_gs_size, max_gs_size,
#'               pvalue_cutoff.
#' @return A gseaResult object, or NULL on failure / insufficient genes.
#' @export
run_gsea_go <- function(ranked, ont, orgdb, params) {

  if (length(ranked) < .GSEA_MIN_GENES) {
    message(sprintf("run_gsea_go [%s]: fewer than %d genes — skipping.",
                    ont, .GSEA_MIN_GENES))
    return(NULL)
  }

  tryCatch({
    set.seed(params$seed)
    gseGO(
      geneList      = ranked,
      OrgDb         = orgdb,
      ont           = ont,
      keyType       = "SYMBOL",
      pAdjustMethod = "BH",
      pvalueCutoff  = params$pvalue_cutoff,
      minGSSize     = params$min_gs_size,
      maxGSSize     = params$max_gs_size,
      nPermSimple   = params$n_perm
    )
  }, error = function(e) {
    message(sprintf("run_gsea_go [%s]: %s", ont, conditionMessage(e)))
    NULL
  })
}


#' Run GSEA with KEGG (gseKEGG)
#'
#' @param ranked_entrez Named numeric vector (Entrez IDs -> ranking metric),
#'                      sorted in decreasing order.
#' @param organism      KEGG organism code (e.g. "hsa", "mmu").
#' @param params        List with n_perm, seed, min_gs_size, max_gs_size,
#'                      pvalue_cutoff.
#' @return A gseaResult object, or NULL on failure / insufficient genes.
#' @export
run_gsea_kegg <- function(ranked_entrez, organism, params) {

  if (length(ranked_entrez) < .GSEA_MIN_GENES) {
    message(sprintf("run_gsea_kegg: fewer than %d Entrez IDs — skipping.",
                    .GSEA_MIN_GENES))
    return(NULL)
  }

  tryCatch({
    set.seed(params$seed)
    gseKEGG(
      geneList      = ranked_entrez,
      organism      = organism,
      pAdjustMethod = "BH",
      pvalueCutoff  = params$pvalue_cutoff,
      minGSSize     = params$min_gs_size,
      maxGSSize     = params$max_gs_size,
      nPermSimple   = params$n_perm
    )
  }, error = function(e) {
    message(sprintf("run_gsea_kegg: %s", conditionMessage(e)))
    NULL
  })
}


#' Run GSEA with Reactome (gsePathway)
#'
#' @param ranked_entrez Named numeric vector (Entrez IDs -> ranking metric),
#'                      sorted in decreasing order.
#' @param organism      Reactome organism label ("human" or "mouse").
#' @param params        List with n_perm, seed, min_gs_size, max_gs_size,
#'                      pvalue_cutoff.
#' @return A gseaResult object, or NULL on failure / insufficient genes.
#' @export
run_gsea_reactome <- function(ranked_entrez, organism, params) {

  if (length(ranked_entrez) < .GSEA_MIN_GENES) {
    message(sprintf("run_gsea_reactome: fewer than %d Entrez IDs — skipping.",
                    .GSEA_MIN_GENES))
    return(NULL)
  }

  tryCatch({
    set.seed(params$seed)
    ReactomePA::gsePathway(
      geneList      = ranked_entrez,
      organism      = organism,
      pAdjustMethod = "BH",
      pvalueCutoff  = params$pvalue_cutoff,
      minGSSize     = params$min_gs_size,
      maxGSSize     = params$max_gs_size,
      nPermSimple   = params$n_perm
    )
  }, error = function(e) {
    message(sprintf("run_gsea_reactome: %s", conditionMessage(e)))
    NULL
  })
}


#' Run GSEA with MSigDB Hallmark gene sets
#'
#' @param ranked Named numeric vector (gene symbols -> ranking metric),
#'               sorted in decreasing order.
#' @param species msigdbr species name ("Homo sapiens" or "Mus musculus").
#' @param params  List with n_perm, seed, min_gs_size, max_gs_size,
#'                pvalue_cutoff.
#' @return A gseaResult object, or NULL on failure / insufficient genes.
#' @export
run_gsea_hallmarks <- function(ranked, species, params) {

  if (length(ranked) < .GSEA_MIN_GENES) {
    message(sprintf("run_gsea_hallmarks: fewer than %d genes — skipping.",
                    .GSEA_MIN_GENES))
    return(NULL)
  }

  tryCatch({
    hallmark_sets <- msigdbr(species = species, category = "H")
    hallmark_t2g  <- hallmark_sets[, c("gs_name", "gene_symbol")]

    set.seed(params$seed)
    GSEA(
      geneList      = ranked,
      TERM2GENE     = hallmark_t2g,
      pAdjustMethod = "BH",
      pvalueCutoff  = params$pvalue_cutoff,
      minGSSize     = params$min_gs_size,
      maxGSSize     = params$max_gs_size,
      nPermSimple   = params$n_perm
    )
  }, error = function(e) {
    message(sprintf("run_gsea_hallmarks: %s", conditionMessage(e)))
    NULL
  })
}


# =============================================================================
# Orchestrators
# =============================================================================

#' Convert gene symbols to Entrez IDs
#'
#' @param symbols Character vector of gene symbols.
#' @param orgdb   OrgDb annotation package name.
#' @return Character vector of Entrez IDs (NAs removed).
#' @keywords internal
.symbols_to_entrez <- function(symbols, orgdb) {
  tryCatch({
    mapped <- clusterProfiler::bitr(
      symbols,
      fromType = "SYMBOL",
      toType   = "ENTREZID",
      OrgDb    = orgdb
    )
    unique(mapped$ENTREZID)
  }, error = function(e) {
    message(sprintf(".symbols_to_entrez: %s", conditionMessage(e)))
    character(0)
  })
}


#' Resolve msigdbr species name from config species string
#'
#' @param cfg Config list (must contain \code{species}).
#' @return A character string suitable for msigdbr ("Homo sapiens" or
#'         "Mus musculus").
#' @keywords internal
.resolve_msigdbr_species <- function(cfg) {
  switch(tolower(cfg$species),
    human = "Homo sapiens",
    mouse = "Mus musculus",
    stop(sprintf("Unsupported species for msigdbr: '%s'", cfg$species),
         call. = FALSE)
  )
}


#' Run all ORA analyses across configured databases
#'
#' For each database listed in \code{cfg$databases}, runs ORA on up-regulated
#' and down-regulated gene sets separately.  Symbol-to-Entrez conversion is
#' handled internally for KEGG and Reactome.
#'
#' @param deg_splits Named list with elements \code{up} and \code{down}, each
#'                   a character vector of gene symbols.
#' @param universe   Character vector of background gene symbols.
#' @param cfg        Full pipeline config list (needs \code{databases},
#'                   \code{enrichment}, \code{orgdb}, \code{kegg_organism},
#'                   \code{reactome_organism}, \code{species}).
#' @return A named list of enrichResult objects keyed by
#'         \code{"<database>_<direction>"} (e.g. "GO_BP_up", "KEGG_down").
#'         Entries are NULL when the analysis was skipped or failed.
#' @export
run_all_ora <- function(deg_splits, universe, cfg) {

  params    <- cfg$enrichment
  databases <- cfg$databases
  results   <- list()

  # Pre-compute Entrez mappings only if KEGG or Reactome are requested
  needs_entrez <- any(c("KEGG", "Reactome") %in% databases)
  entrez_up    <- character(0)
  entrez_down  <- character(0)
  entrez_univ  <- character(0)

  if (needs_entrez) {
    entrez_up   <- .symbols_to_entrez(deg_splits$up,   cfg$orgdb)
    entrez_down <- .symbols_to_entrez(deg_splits$down,  cfg$orgdb)
    entrez_univ <- .symbols_to_entrez(universe,          cfg$orgdb)
  }

  msigdbr_species <- NULL
  if (any(c("Hallmarks", "WikiPathways") %in% databases)) {
    msigdbr_species <- .resolve_msigdbr_species(cfg)
  }

  for (db in databases) {

    for (direction in c("up", "down")) {

      key   <- paste(db, direction, sep = "_")
      genes <- deg_splits[[direction]]

      res <- switch(db,

        GO_BP = run_ora_go(genes, universe, "BP", cfg$orgdb, params),
        GO_CC = run_ora_go(genes, universe, "CC", cfg$orgdb, params),
        GO_MF = run_ora_go(genes, universe, "MF", cfg$orgdb, params),

        KEGG = {
          ids  <- if (direction == "up") entrez_up else entrez_down
          run_ora_kegg(ids, entrez_univ, cfg$kegg_organism, params)
        },

        Reactome = {
          ids  <- if (direction == "up") entrez_up else entrez_down
          run_ora_reactome(ids, entrez_univ, cfg$reactome_organism, params)
        },

        Hallmarks     = run_ora_hallmarks(genes, universe, msigdbr_species, params),
        WikiPathways  = run_ora_wp(genes, universe, msigdbr_species, params),

        {
          message(sprintf("run_all_ora: unknown database '%s' — skipping.", db))
          NULL
        }
      )

      results[[key]] <- res
    }
  }

  results
}


#' Run all GSEA analyses across configured databases
#'
#' @param ranked_symbol Named numeric vector of gene symbols -> ranking metric,
#'                      sorted in decreasing order.  May be NULL if the contrast
#'                      does not provide a full ranked list.
#' @param ranked_entrez Named numeric vector of Entrez IDs -> ranking metric,
#'                      sorted in decreasing order.  May be NULL.
#' @param cfg           Full pipeline config list.
#' @return A named list of gseaResult objects keyed by database name.
#'         Entries are NULL when the analysis was skipped or failed.
#' @export
run_all_gsea <- function(ranked_symbol, ranked_entrez, cfg) {

  params    <- cfg$enrichment
  databases <- cfg$databases
  results   <- list()

  msigdbr_species <- NULL
  if (any(c("Hallmarks") %in% databases)) {
    msigdbr_species <- .resolve_msigdbr_species(cfg)
  }

  for (db in databases) {

    res <- switch(db,

      GO_BP = {
        if (is.null(ranked_symbol)) NULL
        else run_gsea_go(ranked_symbol, "BP", cfg$orgdb, params)
      },

      GO_CC = {
        if (is.null(ranked_symbol)) NULL
        else run_gsea_go(ranked_symbol, "CC", cfg$orgdb, params)
      },

      GO_MF = {
        if (is.null(ranked_symbol)) NULL
        else run_gsea_go(ranked_symbol, "MF", cfg$orgdb, params)
      },

      KEGG = {
        if (is.null(ranked_entrez)) NULL
        else run_gsea_kegg(ranked_entrez, cfg$kegg_organism, params)
      },

      Reactome = {
        if (is.null(ranked_entrez)) NULL
        else run_gsea_reactome(ranked_entrez, cfg$reactome_organism, params)
      },

      Hallmarks = {
        if (is.null(ranked_symbol)) NULL
        else run_gsea_hallmarks(ranked_symbol, msigdbr_species, params)
      },

      # WikiPathways GSEA not implemented — skip silently
      WikiPathways = NULL,

      {
        message(sprintf("run_all_gsea: unknown database '%s' — skipping.", db))
        NULL
      }
    )

    results[[db]] <- res
  }

  results
}


# =============================================================================
# I/O helpers
# =============================================================================

#' Save an enrichment result to disk
#'
#' Writes the enrichResult/gseaResult as RDS and CSV.  If the result contains
#' at least one significant term, a dotplot is also saved in the requested
#' formats.
#'
#' @param result  An enrichResult or gseaResult object (or NULL, in which case
#'                the function returns invisibly).
#' @param out_dir Character. Directory where files will be written.
#' @param name    Character. Base file name (without extension).
#' @param formats Character vector of plot formats (e.g. c("pdf", "png")).
#' @return NULL (called for side effects).
#' @export
save_enrichment <- function(result, out_dir, name, formats = c("pdf", "png")) {

  if (is.null(result)) return(invisible(NULL))

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # -- RDS --
  rds_path <- file.path(out_dir, paste0(name, ".rds"))
  tryCatch(
    saveRDS(result, rds_path),
    error = function(e) message(sprintf("save_enrichment [RDS]: %s", conditionMessage(e)))
  )

  # -- CSV --
  tryCatch({
    df <- as.data.frame(result)
    if (nrow(df) > 0) {
      csv_path <- file.path(out_dir, paste0(name, ".csv"))
      write.csv(df, csv_path, row.names = FALSE)
    }
  }, error = function(e) {
    message(sprintf("save_enrichment [CSV]: %s", conditionMessage(e)))
  })

  # -- Dotplot --
  tryCatch({
    df <- as.data.frame(result)
    if (nrow(df) > 0) {
      n_show <- min(20, nrow(df))
      p <- dotplot(result, showCategory = n_show) +
        theme_minimal(base_size = 12)

      for (fmt in formats) {
        plot_path <- file.path(out_dir, paste0(name, "_dotplot.", fmt))
        ggsave(plot_path, plot = p, width = 10, height = 8, dpi = 300)
      }
    }
  }, error = function(e) {
    message(sprintf("save_enrichment [dotplot]: %s", conditionMessage(e)))
  })

  invisible(NULL)
}

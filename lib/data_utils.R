# data_utils.R — DEG loading, column standardization, classification, ranked lists
# Part of the SBMA RNAseq Generalized Pipeline
# ══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(clusterProfiler)
})

# ── File I/O helpers ─────────────────────────────────────────────────────────

#' Detect file format from extension and read into a tibble
#' @param path Path to the DEG file
#' @return A tibble of the raw contents
#' @keywords internal
.read_by_extension <- function(path) {
  if (!file.exists(path)) {
    stop("DEG file not found: ", path, call. = FALSE)
  }

  ext <- tolower(tools::file_ext(path))

  tryCatch(
    switch(ext,
      csv  = readr::read_csv(path, show_col_types = FALSE),
      tsv  = readr::read_tsv(path, show_col_types = FALSE),
      txt  = readr::read_tsv(path, show_col_types = FALSE),
      xlsx = {
        if (!requireNamespace("readxl", quietly = TRUE)) {
          stop("Package 'readxl' is required to read .xlsx files. ",
               "Install it with install.packages('readxl').", call. = FALSE)
        }
        readxl::read_excel(path)
      },
      stop("Unsupported file extension: .", ext,
           ". Supported formats: .csv, .tsv, .txt, .xlsx", call. = FALSE)
    ),
    error = function(e) {
      stop("Failed to read '", path, "': ", conditionMessage(e), call. = FALSE)
    }
  )
}

# ── Core functions ───────────────────────────────────────────────────────────

#' Load DEG table and standardize column names
#'
#' Reads a DEG results file (CSV, TSV, TXT, or XLSX) and renames columns
#' according to the mapping provided in the contrast configuration. The
#' returned tibble always has columns: gene_symbol, log2fc, pvalue, padj.
#'
#' @param contrast_cfg A list (typically from config.yaml) with at least:
#'   \describe{
#'     \item{deg_file}{Character. Path to the DEG results file.}
#'     \item{columns}{Named list mapping standardized names to original names,
#'       e.g. \code{list(gene_symbol = "Symbol", log2fc = "logFC",
#'       pvalue = "PValue", padj = "FDR")}.}
#'   }
#' @return A tibble with columns: gene_symbol (character), log2fc (numeric),
#'   pvalue (numeric), padj (numeric). Rows with NA gene_symbol are removed
#'   and duplicates are resolved by keeping the row with the lowest padj.
load_deg_table <- function(contrast_cfg) {
  # --- Validate inputs ---
  required_keys <- c("deg_file", "columns")
  missing_keys  <- setdiff(required_keys, names(contrast_cfg))
  if (length(missing_keys) > 0) {
    stop("contrast_cfg is missing required keys: ",
         paste(missing_keys, collapse = ", "), call. = FALSE)
  }

  col_map        <- contrast_cfg$columns
  standard_names <- c("gene_symbol", "log2fc", "pvalue", "padj")
  missing_cols   <- setdiff(standard_names, names(col_map))
  if (length(missing_cols) > 0) {
    stop("Column mapping is missing entries for: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  # --- Read raw data ---
  raw <- .read_by_extension(contrast_cfg$deg_file)

  # --- Validate that mapped source columns exist in the file ---
  missing_in_file <- setdiff(unlist(col_map[standard_names]), colnames(raw))
  if (length(missing_in_file) > 0) {
    stop("Column(s) not found in '", contrast_cfg$deg_file, "': ",
         paste(missing_in_file, collapse = ", "),
         ". Available columns: ", paste(colnames(raw), collapse = ", "),
         call. = FALSE)
  }

  # --- Build rename vector: new_name = old_name ---
  rename_vec <- setNames(
    vapply(standard_names, function(nm) col_map[[nm]], character(1)),
    standard_names
  )

  # --- Select and rename ---
  df <- raw %>%
    dplyr::select(dplyr::all_of(unname(rename_vec))) %>%
    dplyr::rename(!!!rename_vec)

  # --- Coerce types ---
  df <- df %>%
    dplyr::mutate(
      gene_symbol = as.character(gene_symbol),
      log2fc      = as.numeric(log2fc),
      pvalue      = as.numeric(pvalue),
      padj        = as.numeric(padj)
    )

  # --- Clean: remove NA symbols, deduplicate keeping lowest padj ---
  df <- df %>%
    dplyr::filter(!is.na(gene_symbol), gene_symbol != "") %>%
    dplyr::arrange(padj) %>%
    dplyr::distinct(gene_symbol, .keep_all = TRUE)

  df
}


#' Classify DEGs as Upregulated, Downregulated, or NS (not significant)
#'
#' Adds a \code{diffexpressed} column based on fold-change and adjusted
#' p-value thresholds.
#'
#' @param df A standardized DEG tibble (must contain log2fc and padj columns).
#' @param lfc_thresh Absolute log2 fold-change threshold (default 1).
#' @param padj_thresh Adjusted p-value threshold (default 0.05).
#' @return The input tibble with an added \code{diffexpressed} factor column
#'   with levels: "Downregulated", "NS", "Upregulated".
classify_deg <- function(df, lfc_thresh = 1, padj_thresh = 0.05) {
  stopifnot(
    is.data.frame(df),
    all(c("log2fc", "padj") %in% colnames(df)),
    is.numeric(lfc_thresh), length(lfc_thresh) == 1, lfc_thresh >= 0,
    is.numeric(padj_thresh), length(padj_thresh) == 1,
    padj_thresh > 0, padj_thresh <= 1
  )

  df %>%
    dplyr::mutate(
      diffexpressed = dplyr::case_when(
        padj < padj_thresh & log2fc >  lfc_thresh ~ "Upregulated",
        padj < padj_thresh & log2fc < -lfc_thresh ~ "Downregulated",
        TRUE                                       ~ "NS"
      ),
      diffexpressed = factor(diffexpressed,
                             levels = c("Downregulated", "NS", "Upregulated"))
    )
}


#' Split DEGs into gene-symbol vectors by direction
#'
#' @param df A classified DEG tibble (output of \code{classify_deg}; must
#'   contain gene_symbol and diffexpressed columns).
#' @return A named list with elements:
#'   \describe{
#'     \item{up}{Character vector of upregulated gene symbols.}
#'     \item{down}{Character vector of downregulated gene symbols.}
#'     \item{all_sig}{Character vector of all significant gene symbols.}
#'     \item{background}{Character vector of all gene symbols in the table.}
#'   }
split_deg <- function(df) {
  stopifnot(
    is.data.frame(df),
    all(c("gene_symbol", "diffexpressed") %in% colnames(df))
  )

  list(
    up         = df$gene_symbol[df$diffexpressed == "Upregulated"],
    down       = df$gene_symbol[df$diffexpressed == "Downregulated"],
    all_sig    = df$gene_symbol[df$diffexpressed != "NS"],
    background = df$gene_symbol
  )
}


#' Create a ranked gene list for GSEA (named by gene symbol)
#'
#' The ranking metric is \code{sign(log2fc) * -log10(pvalue)}, which combines
#' direction and statistical significance. P-values are floored at 1e-300 to
#' avoid infinite values.
#'
#' @param df A standardized DEG tibble (must contain gene_symbol, log2fc,
#'   pvalue).
#' @return A named numeric vector sorted in descending order. Names are gene
#'   symbols. Genes with NA metric or duplicated symbols are removed (keeping
#'   the entry with the largest absolute metric).
make_ranked_list <- function(df) {
  stopifnot(
    is.data.frame(df),
    all(c("gene_symbol", "log2fc", "pvalue") %in% colnames(df))
  )

  ranked <- df %>%
    dplyr::filter(!is.na(gene_symbol), gene_symbol != "",
                  !is.na(log2fc), !is.na(pvalue)) %>%
    dplyr::mutate(
      rank_metric = sign(log2fc) * -log10(pmax(pvalue, 1e-300))
    ) %>%
    # Deduplicate: keep the entry with the largest absolute metric per gene
    dplyr::arrange(dplyr::desc(abs(rank_metric))) %>%
    dplyr::distinct(gene_symbol, .keep_all = TRUE)

  rl <- setNames(ranked$rank_metric, ranked$gene_symbol)
  sort(rl, decreasing = TRUE)
}


#' Create a ranked gene list with Entrez IDs for GSEA
#'
#' Converts gene symbols to Entrez IDs using \code{clusterProfiler::bitr},
#' then builds a ranked list with the same metric as \code{make_ranked_list}.
#'
#' @param df A standardized DEG tibble (must contain gene_symbol, log2fc,
#'   pvalue).
#' @param orgdb Character string naming the OrgDb annotation package
#'   (e.g., \code{"org.Hs.eg.db"}).
#' @return A named numeric vector (names = Entrez IDs) sorted descending.
#'   Unmappable symbols are silently dropped. Duplicate Entrez IDs are resolved
#'   by keeping the entry with the largest absolute metric.
make_ranked_entrez <- function(df, orgdb) {
  stopifnot(
    is.data.frame(df),
    all(c("gene_symbol", "log2fc", "pvalue") %in% colnames(df)),
    is.character(orgdb), length(orgdb) == 1
  )

  # Build symbol-level ranked table first
  ranked <- df %>%
    dplyr::filter(!is.na(gene_symbol), gene_symbol != "",
                  !is.na(log2fc), !is.na(pvalue)) %>%
    dplyr::mutate(
      rank_metric = sign(log2fc) * -log10(pmax(pvalue, 1e-300))
    ) %>%
    dplyr::arrange(dplyr::desc(abs(rank_metric))) %>%
    dplyr::distinct(gene_symbol, .keep_all = TRUE)

  # Map symbols to Entrez IDs
  mapping <- symbols_to_entrez(ranked$gene_symbol, orgdb)

  if (length(mapping) == 0) {
    warning("No gene symbols could be mapped to Entrez IDs.", call. = FALSE)
    return(setNames(numeric(0), character(0)))
  }

  # Join mapping onto ranked table
  map_df <- tibble::tibble(
    gene_symbol = names(mapping),
    entrez_id   = unname(mapping)
  )

  merged <- ranked %>%
    dplyr::inner_join(map_df, by = "gene_symbol") %>%
    # Resolve duplicate Entrez IDs (already sorted by descending abs metric)
    dplyr::distinct(entrez_id, .keep_all = TRUE)

  rl <- setNames(merged$rank_metric, merged$entrez_id)
  sort(rl, decreasing = TRUE)
}


#' Convert gene symbols to Entrez IDs
#'
#' Wrapper around \code{clusterProfiler::bitr} that suppresses warnings about
#' unmapped identifiers and returns a convenient named vector.
#'
#' @param symbols Character vector of gene symbols.
#' @param orgdb Character string naming the OrgDb annotation package
#'   (e.g., \code{"org.Hs.eg.db"}).
#' @return A named character vector where names are the original gene symbols
#'   and values are the corresponding Entrez IDs. Symbols that could not be
#'   mapped are silently excluded. If the OrgDb package is not installed, an
#'   informative error is raised.
symbols_to_entrez <- function(symbols, orgdb) {
  stopifnot(
    is.character(symbols),
    is.character(orgdb), length(orgdb) == 1
  )

  if (!requireNamespace(orgdb, quietly = TRUE)) {
    stop("Annotation package '", orgdb, "' is not installed. ",
         "Install it with BiocManager::install('", orgdb, "').",
         call. = FALSE)
  }

  symbols <- unique(symbols[!is.na(symbols) & symbols != ""])
  if (length(symbols) == 0) {
    return(setNames(character(0), character(0)))
  }

  converted <- tryCatch(
    suppressWarnings(
      clusterProfiler::bitr(
        symbols,
        fromType = "SYMBOL",
        toType   = "ENTREZID",
        OrgDb    = orgdb
      )
    ),
    error = function(e) {
      warning("bitr conversion failed: ", conditionMessage(e), call. = FALSE)
      return(data.frame(SYMBOL = character(0), ENTREZID = character(0)))
    }
  )

  if (nrow(converted) == 0) {
    return(setNames(character(0), character(0)))
  }

  # Deduplicate: if one symbol maps to multiple Entrez IDs, keep the first
  converted <- converted[!duplicated(converted$SYMBOL), ]

  setNames(converted$ENTREZID, converted$SYMBOL)
}

# crosstalk_utils.R — Mitochondrial-purinergic crosstalk analysis utilities
# Part of the SBMA RNAseq Generalized Pipeline
# ══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(glue)
})

# ── Overlap between purinergic and mitochondrial gene sets ─────────────────

#' Compute gene overlap between purinergic and mitochondrial DEG tables
#'
#' Identifies genes present in both the purinergic and mitochondrial merged
#' DEG tables (output of \code{merge_with_deg}), returning a combined tibble
#' that preserves the category membership from each source.
#'
#' @param puri_df Tibble from \code{merge_with_deg} for purinergic genes.
#'   Must contain: gene_symbol, log2fc, padj, sig, direction.
#' @param mito_df Tibble from \code{merge_with_deg} for mitochondrial genes.
#'   Must contain: gene_symbol, log2fc, padj, sig, direction.
#' @return A tibble with columns: gene_symbol, puri_category, mito_category,
#'   log2fc, padj, sig, direction. Returns an empty tibble with the correct
#'   schema if either input is empty or there is no overlap.
compute_overlap <- function(puri_df, mito_df) {
  stopifnot(is.data.frame(puri_df), is.data.frame(mito_df))

  .required_cols <- c("gene_symbol", "log2fc", "padj", "sig", "direction")

  missing_puri <- setdiff(.required_cols, colnames(puri_df))
  if (length(missing_puri) > 0) {
    stop("puri_df is missing columns: ",
         paste(missing_puri, collapse = ", "), call. = FALSE)
  }
  missing_mito <- setdiff(.required_cols, colnames(mito_df))
  if (length(missing_mito) > 0) {
    stop("mito_df is missing columns: ",
         paste(missing_mito, collapse = ", "), call. = FALSE)
  }

  # Empty-input guard: return correctly-typed empty tibble
  .empty_result <- tibble::tibble(
    gene_symbol   = character(0),
    puri_category = character(0),
    mito_category = character(0),
    log2fc        = numeric(0),
    padj          = numeric(0),
    sig           = logical(0),
    direction     = character(0)
  )

  if (nrow(puri_df) == 0 || nrow(mito_df) == 0) {
    return(.empty_result)
  }

  # Determine category column names (first column that is not in .required_cols

  # and looks like a category). Convention: merge_with_deg adds a "category"
  # column; fall back to NA if absent.
  puri_cat_col <- intersect(c("category", "subcategory", "family"), colnames(puri_df))
  mito_cat_col <- intersect(c("category", "subcategory", "family"), colnames(mito_df))

  puri_prep <- puri_df %>%
    dplyr::select(gene_symbol,
                  puri_category = dplyr::any_of(puri_cat_col[1]),
                  log2fc, padj, sig, direction) %>%
    dplyr::distinct(gene_symbol, .keep_all = TRUE)

  mito_prep <- mito_df %>%
    dplyr::select(gene_symbol,
                  mito_category = dplyr::any_of(mito_cat_col[1])) %>%
    dplyr::distinct(gene_symbol, .keep_all = TRUE)

  # If no category column was found, add NA placeholders

  if (!"puri_category" %in% colnames(puri_prep)) {
    puri_prep <- dplyr::mutate(puri_prep, puri_category = NA_character_,
                               .after = gene_symbol)
  }
  if (!"mito_category" %in% colnames(mito_prep)) {
    mito_prep <- dplyr::mutate(mito_prep, mito_category = NA_character_,
                               .after = gene_symbol)
  }

  overlap <- dplyr::inner_join(puri_prep, mito_prep, by = "gene_symbol")

  if (nrow(overlap) == 0) {
    return(.empty_result)
  }

  overlap %>%
    dplyr::select(gene_symbol, puri_category, mito_category,
                  log2fc, padj, sig, direction)
}


# ── Build crosstalk dataframe ──────────────────────────────────────────────

#' Build a crosstalk analysis dataframe
#'
#' Left-joins a predefined crosstalk gene list (from \code{get_crosstalk_genes})
#' with a standardized DEG table, adding significance and direction annotations.
#'
#' @param deg_df Standardized DEG tibble with at least: gene_symbol, log2fc,
#'   pvalue, padj.
#' @param crosstalk_genes_df Tibble with columns: gene_symbol, function_group.
#'   Defines the crosstalk gene universe and their functional grouping.
#' @return A tibble with columns: gene_symbol, function_group, log2fc, pvalue,
#'   padj, sig, direction. Genes in the crosstalk list that are absent from the
#'   DEG table will have NA for expression values and sig = FALSE.
build_crosstalk_df <- function(deg_df, crosstalk_genes_df) {
  stopifnot(is.data.frame(deg_df), is.data.frame(crosstalk_genes_df))

  if (!all(c("gene_symbol", "log2fc", "pvalue", "padj") %in% colnames(deg_df))) {
    stop("deg_df must contain columns: gene_symbol, log2fc, pvalue, padj",
         call. = FALSE)
  }
  if (!all(c("gene_symbol", "function_group") %in% colnames(crosstalk_genes_df))) {
    stop("crosstalk_genes_df must contain columns: gene_symbol, function_group",
         call. = FALSE)
  }

  # Empty-input guard
  if (nrow(crosstalk_genes_df) == 0) {
    return(tibble::tibble(
      gene_symbol    = character(0),
      function_group = character(0),
      log2fc         = numeric(0),
      pvalue         = numeric(0),
      padj           = numeric(0),
      sig            = logical(0),
      direction      = character(0)
    ))
  }

  # Deduplicate DEG table to one row per gene (lowest padj wins)
  deg_dedup <- deg_df %>%
    dplyr::filter(!is.na(gene_symbol), gene_symbol != "") %>%
    dplyr::arrange(padj) %>%
    dplyr::distinct(gene_symbol, .keep_all = TRUE) %>%
    dplyr::select(gene_symbol, log2fc, pvalue, padj)

  result <- crosstalk_genes_df %>%
    dplyr::select(gene_symbol, function_group) %>%
    dplyr::left_join(deg_dedup, by = "gene_symbol") %>%
    dplyr::mutate(
      sig = dplyr::if_else(!is.na(padj) & padj < 0.05, TRUE, FALSE),
      direction = dplyr::case_when(
        is.na(log2fc)     ~ NA_character_,
        sig & log2fc > 0  ~ "up",
        sig & log2fc < 0  ~ "down",
        TRUE              ~ "ns"
      )
    )

  result
}


# ── Mechanism summary generation ───────────────────────────────────────────

#' Auto-generate a mechanistic narrative from crosstalk data
#'
#' Inspects each functional group in the crosstalk dataframe for significant
#' genes and their direction, then composes a plain-text summary with
#' data-driven mechanistic interpretation.
#'
#' @param crosstalk_df Tibble as returned by \code{build_crosstalk_df}. Must
#'   contain: function_group, gene_symbol, log2fc, sig, direction.
#' @param contrast_label Character string identifying the contrast
#'   (e.g., "SBMA vs Control").
#' @return A single character string containing the full mechanism summary.
#'   Returns a brief "no data" message if the input is empty.
generate_mechanism_summary <- function(crosstalk_df, contrast_label) {
  stopifnot(is.data.frame(crosstalk_df), is.character(contrast_label),
            length(contrast_label) == 1)

  if (nrow(crosstalk_df) == 0) {
    return(glue::glue("In {contrast_label}: No crosstalk genes available for analysis."))
  }

  required <- c("function_group", "gene_symbol", "log2fc", "sig", "direction")
  missing  <- setdiff(required, colnames(crosstalk_df))
  if (length(missing) > 0) {
    stop("crosstalk_df is missing columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  # ---- Per-group summary lines ----
  group_stats <- crosstalk_df %>%
    dplyr::group_by(function_group) %>%
    dplyr::summarise(
      total  = dplyr::n(),
      n_sig  = sum(sig, na.rm = TRUE),
      n_up   = sum(direction == "up", na.rm = TRUE),
      n_down = sum(direction == "down", na.rm = TRUE),
      mean_lfc = mean(log2fc, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(function_group)

  group_lines <- purrr::map_chr(seq_len(nrow(group_stats)), function(i) {
    row <- group_stats[i, ]
    lfc_str <- if (is.nan(row$mean_lfc)) "NA" else sprintf("%.2f", row$mean_lfc)
    glue::glue("- {row$function_group}: {row$n_sig} of {row$total} genes ",
               "significant ({row$n_up} up, {row$n_down} down). ",
               "Mean LFC = {lfc_str}.")
  })

  # ---- Mechanistic interpretation ----
  interpretations <- character(0)

  # Helper: check if a group has predominantly a given direction among sig genes
  .group_direction <- function(group_pattern, dir) {
    grp <- group_stats %>%
      dplyr::filter(grepl(group_pattern, function_group, ignore.case = TRUE))
    if (nrow(grp) == 0) return(FALSE)
    grp <- dplyr::summarise(grp,
                            n_dir = sum(get(paste0("n_", dir)), na.rm = TRUE),
                            n_sig = sum(n_sig, na.rm = TRUE))
    grp$n_sig > 0 && grp$n_dir > 0
  }

  .group_has_sig <- function(group_pattern) {
    grp <- group_stats %>%
      dplyr::filter(grepl(group_pattern, function_group, ignore.case = TRUE))
    if (nrow(grp) == 0) return(FALSE)
    sum(grp$n_sig, na.rm = TRUE) > 0
  }

  # Pattern 1: ATP production down AND purinergic receptors up
  if (.group_direction("ATP.production|oxidative.phosphorylation|OXPHOS|electron.transport",
                       "down") &&
      .group_direction("purinergic.receptor|P2[XY]|P1", "up")) {
    interpretations <- c(interpretations,
      "Pattern suggests compensatory purinergic signaling in response to mitochondrial ATP deficit.")
  }

  # Pattern 2: Inflammasome up AND ROS defense down
  if (.group_direction("inflammasome|NLRP|IL.1|caspase", "up") &&
      .group_direction("ROS.defense|antioxidant|SOD|catalase|glutathione", "down")) {
    interpretations <- c(interpretations,
      "Pattern suggests inflammasome activation potentially driven by oxidative stress.")
  }

  # Pattern 3: ATP release up
  if (.group_direction("ATP.release|pannexin|connexin|VDAC|channel", "up")) {
    interpretations <- c(interpretations,
      "Elevated ATP release channels suggest increased extracellular ATP signaling.")
  }

  if (length(interpretations) == 0) {
    interpretations <- "No canonical crosstalk patterns detected at current significance thresholds."
  }

  # ---- Assemble full summary ----
  header <- glue::glue("In {contrast_label}:")
  body   <- paste(group_lines, collapse = "\n")
  interp <- paste(interpretations, collapse = "\n")

  paste(header, body, "", "Mechanistic interpretation:", interp, sep = "\n")
}


# ── Group-level summary ────────────────────────────────────────────────────

#' Summarize crosstalk genes by functional group
#'
#' Produces a one-row-per-group summary with gene counts, significance
#' breakdown, mean fold-change, and a comma-separated list of significant
#' gene symbols.
#'
#' @param crosstalk_df Tibble as returned by \code{build_crosstalk_df}. Must
#'   contain: function_group, gene_symbol, log2fc, sig, direction.
#' @return A summary tibble with columns: function_group, n_genes, n_sig,
#'   n_up, n_down, mean_lfc, genes_sig. Returns an empty tibble with the
#'   correct schema if the input is empty.
summarize_crosstalk_by_group <- function(crosstalk_df) {
  stopifnot(is.data.frame(crosstalk_df))

  .empty_summary <- tibble::tibble(
    function_group = character(0),
    n_genes        = integer(0),
    n_sig          = integer(0),
    n_up           = integer(0),
    n_down         = integer(0),
    mean_lfc       = numeric(0),
    genes_sig      = character(0)
  )

  if (nrow(crosstalk_df) == 0) {
    return(.empty_summary)
  }

  required <- c("function_group", "gene_symbol", "log2fc", "sig", "direction")
  missing  <- setdiff(required, colnames(crosstalk_df))
  if (length(missing) > 0) {
    stop("crosstalk_df is missing columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  crosstalk_df %>%
    dplyr::group_by(function_group) %>%
    dplyr::summarise(
      n_genes   = dplyr::n(),
      n_sig     = sum(sig, na.rm = TRUE),
      n_up      = sum(direction == "up", na.rm = TRUE),
      n_down    = sum(direction == "down", na.rm = TRUE),
      mean_lfc  = mean(log2fc, na.rm = TRUE),
      genes_sig = paste(gene_symbol[which(sig)], collapse = ", "),
      .groups   = "drop"
    ) %>%
    dplyr::mutate(
      genes_sig = dplyr::if_else(genes_sig == "", NA_character_, genes_sig)
    )
}

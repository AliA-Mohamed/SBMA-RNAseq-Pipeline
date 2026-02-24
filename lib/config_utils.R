# config_utils.R
# Config loading, validation, and path helpers for the SBMA RNAseq pipeline.
# ---------------------------------------------------------------------------

library(yaml)

# ---- Required top-level fields and their expected types ----
.REQUIRED_TOP_FIELDS <- c("species", "orgdb", "kegg_organism",
                          "reactome_organism", "thresholds", "datasets")

.REQUIRED_THRESHOLD_FIELDS <- c("lfc", "padj")

.REQUIRED_ENRICHMENT_FIELDS <- c("pvalue_cutoff", "qvalue_cutoff",
                                 "min_gs_size", "max_gs_size",
                                 "simplify_cutoff", "n_perm", "seed")

.REQUIRED_CONTRAST_FIELDS <- c("label", "deg_file", "full_ranked",
                                "androgen_treatment", "columns")

.REQUIRED_COLUMN_FIELDS <- c("gene_symbol", "log2fc", "pvalue", "padj")

.OUTPUT_SUBDIRS <- c("data",
                     "tables",
                     "plots/qc",
                     "plots/volcano",
                     "plots/ora",
                     "plots/gsea",
                     "plots/purinergic",
                     "plots/mitochondrial",
                     "reports")


#' Load and validate pipeline configuration
#'
#' Reads a YAML configuration file and checks that all required sections and
#' fields are present.  Stops with an informative error when anything is
#' missing.
#'
#' @param config_path Character. Path to \code{config.yaml}.
#' @return A validated configuration list (invisibly).
#' @export
load_config <- function(config_path) {


  if (missing(config_path) || is.null(config_path)) {
    stop("load_config: 'config_path' must be provided.", call. = FALSE)
  }

  if (!file.exists(config_path)) {
    stop(sprintf("load_config: config file not found: %s", config_path),
         call. = FALSE)
  }

  cfg <- tryCatch(
    yaml::read_yaml(config_path),
    error = function(e) {
      stop(sprintf("load_config: failed to parse YAML — %s", conditionMessage(e)),
           call. = FALSE)
    }
  )

  # -- Top-level required fields --
  missing_top <- setdiff(.REQUIRED_TOP_FIELDS, names(cfg))
  if (length(missing_top) > 0L) {
    stop(sprintf("load_config: missing required top-level fields: %s",
                 paste(missing_top, collapse = ", ")),
         call. = FALSE)
  }

  # -- Thresholds --
  .validate_section(cfg$thresholds, .REQUIRED_THRESHOLD_FIELDS,
                    section_label = "thresholds")

  # -- Enrichment (optional section, but if present validate fully) --
  if (!is.null(cfg$enrichment)) {
    .validate_section(cfg$enrichment, .REQUIRED_ENRICHMENT_FIELDS,
                      section_label = "enrichment")
  }

  # -- Datasets and contrasts --
  if (!is.list(cfg$datasets) || length(cfg$datasets) == 0L) {
    stop("load_config: 'datasets' must be a non-empty named list.", call. = FALSE)
  }

  for (ds_id in names(cfg$datasets)) {
    ds <- cfg$datasets[[ds_id]]

    if (is.null(ds$tissue_type)) {
      stop(sprintf("load_config: dataset '%s' is missing 'tissue_type'.", ds_id),
           call. = FALSE)
    }

    if (!is.list(ds$contrasts) || length(ds$contrasts) == 0L) {
      stop(sprintf("load_config: dataset '%s' must have at least one contrast.",
                   ds_id),
           call. = FALSE)
    }

    for (ct_id in names(ds$contrasts)) {
      ct <- ds$contrasts[[ct_id]]
      .validate_section(ct, .REQUIRED_CONTRAST_FIELDS,
                        section_label = sprintf("datasets/%s/contrasts/%s",
                                                ds_id, ct_id))
      .validate_section(ct$columns, .REQUIRED_COLUMN_FIELDS,
                        section_label = sprintf("datasets/%s/contrasts/%s/columns",
                                                ds_id, ct_id))
    }
  }

  invisible(cfg)
}


#' Extract a single contrast configuration
#'
#' Retrieves the contrast sub-list for the given dataset and contrast
#' identifiers, attaching useful context fields (\code{dataset_id},
#' \code{contrast_id}, \code{tissue_type}).
#'
#' @param cfg   Config list (as returned by \code{load_config}).
#' @param dataset_id  Character. Name of the dataset.
#' @param contrast_id Character. Name of the contrast within the dataset.
#' @return A list combining the contrast config with \code{dataset_id},
#'   \code{contrast_id}, and dataset-level fields such as \code{tissue_type}.
#' @export
get_contrast <- function(cfg, dataset_id, contrast_id) {

  if (is.null(cfg$datasets[[dataset_id]])) {
    stop(sprintf("get_contrast: dataset '%s' not found in config. Available: %s",
                 dataset_id, paste(names(cfg$datasets), collapse = ", ")),
         call. = FALSE)
  }

  ds <- cfg$datasets[[dataset_id]]

  if (is.null(ds$contrasts[[contrast_id]])) {
    stop(sprintf("get_contrast: contrast '%s' not found in dataset '%s'. Available: %s",
                 contrast_id, dataset_id,
                 paste(names(ds$contrasts), collapse = ", ")),
         call. = FALSE)
  }

  contrast <- ds$contrasts[[contrast_id]]


  # Attach identifiers and dataset-level context
  contrast$dataset_id  <- dataset_id
  contrast$contrast_id <- contrast_id
  contrast$tissue_type <- ds$tissue_type

  contrast
}


#' Build output directory path and create subdirectories
#'
#' Constructs \code{results/<dataset_id>/<contrast_id>} and ensures all
#' standard subdirectories exist.
#'
#' @param cfg         Config list.
#' @param dataset_id  Character. Dataset name.
#' @param contrast_id Character. Contrast name.
#' @return Character string with the base output directory path (invisibly).
#' @export
get_output_dir <- function(cfg, dataset_id, contrast_id) {

  base_dir <- file.path("results", dataset_id, contrast_id)

  for (sub in .OUTPUT_SUBDIRS) {
    dir.create(file.path(base_dir, sub), recursive = TRUE, showWarnings = FALSE)
  }

  invisible(base_dir)
}


#' Get species-specific parameters
#'
#' Extracts species annotation parameters from the config for use in
#' enrichment and annotation steps.
#'
#' @param cfg Config list.
#' @return A named list with elements \code{orgdb}, \code{kegg_organism},
#'   and \code{reactome_organism}.
#' @export
get_species_params <- function(cfg) {

  required <- c("orgdb", "kegg_organism", "reactome_organism")
  missing  <- setdiff(required, names(cfg))

  if (length(missing) > 0L) {
    stop(sprintf("get_species_params: missing species fields in config: %s",
                 paste(missing, collapse = ", ")),
         call. = FALSE)
  }

  list(
    orgdb            = cfg$orgdb,
    kegg_organism    = cfg$kegg_organism,
    reactome_organism = cfg$reactome_organism
  )
}


#' Parse command-line arguments for pipeline scripts
#'
#' Uses \pkg{optparse} to parse \code{--config}, \code{--dataset}, and
#' \code{--contrast} from the command line.
#'
#' @return A named list with elements \code{config}, \code{dataset}, and
#'   \code{contrast}.
#' @export
parse_pipeline_args <- function() {

  if (!requireNamespace("optparse", quietly = TRUE)) {
    stop("parse_pipeline_args: the 'optparse' package is required.", call. = FALSE)
  }

  option_list <- list(
    optparse::make_option(
      c("-c", "--config"),
      type    = "character",
      default = "config/config.yaml",
      help    = "Path to the pipeline config.yaml [default: %default]"
    ),
    optparse::make_option(
      c("-d", "--dataset"),
      type = "character",
      default = NULL,
      help = "Dataset identifier (must match a key under 'datasets' in config)"
    ),
    optparse::make_option(
      c("-t", "--contrast"),
      type = "character",
      default = NULL,
      help = "Contrast identifier (must match a key under the chosen dataset)"
    )
  )

  parser <- optparse::OptionParser(
    usage       = "usage: %prog [options]",
    option_list = option_list,
    description = "SBMA RNAseq Pipeline — run a single contrast analysis"
  )

  args <- optparse::parse_args(parser)

  if (is.null(args$dataset)) {
    stop("parse_pipeline_args: --dataset is required.", call. = FALSE)
  }
  if (is.null(args$contrast)) {
    stop("parse_pipeline_args: --contrast is required.", call. = FALSE)
  }

  list(
    config   = args$config,
    dataset  = args$dataset,
    contrast = args$contrast
  )
}


# ---- Internal helpers ------------------------------------------------------

#' Validate that all expected fields are present in a config section
#' @param section       List to check.
#' @param required      Character vector of required field names.
#' @param section_label Human-readable label for error messages.
#' @return NULL (called for side effects).
#' @keywords internal
.validate_section <- function(section, required, section_label) {
  if (!is.list(section)) {
    stop(sprintf("load_config: '%s' must be a list/mapping.", section_label),
         call. = FALSE)
  }
  missing <- setdiff(required, names(section))
  if (length(missing) > 0L) {
    stop(sprintf("load_config: '%s' is missing fields: %s",
                 section_label, paste(missing, collapse = ", ")),
         call. = FALSE)
  }
  invisible(NULL)
}

#!/usr/bin/env Rscript
# ==============================================================================
# 00_deseq2_from_counts.R — Compute DEGs from raw count matrix using DESeq2
# ==============================================================================
#
# Optional pre-processing step. When raw counts are available, this script
# runs DESeq2 to produce a full genome-wide DEG table with proper statistical
# modeling (negative binomial, Wald test, BH-adjusted p-values, shrunken LFC).
#
# The output replaces the need for pre-computed DEG files and TPM-based
# workarounds, providing statistically rigorous input for the pipeline.
#
# Usage:
#   Rscript scripts/00_deseq2_from_counts.R \
#     --counts data/MyDataset/count_matrix.csv \
#     --samples data/MyDataset/sample_info.csv \
#     --design "~ condition" \
#     --contrast "condition,disease,control" \
#     --output data/MyDataset/deseq2_degs.csv \
#     --min-count 10
#
# Inputs:
#   --counts   : CSV/TSV with gene IDs as first column, samples as remaining
#                columns. Values must be raw integer counts (not TPM/RPKM).
#   --samples  : CSV/TSV with sample metadata. First column = sample IDs
#                (must match count matrix column names). Must contain the
#                variable(s) used in --design.
#   --design   : DESeq2 design formula (e.g., "~ condition" or
#                "~ batch + condition").
#   --contrast : Comma-separated: "variable,numerator,denominator"
#                (e.g., "condition,SBMA,Control").
#   --output   : Path for the output DEG CSV file.
#   --min-count: Minimum total count across all samples to keep a gene
#                (default: 10). Genes below this threshold are filtered.
#
# Outputs:
#   A CSV file with columns: gene_symbol, log2fc, pvalue, padj,
#   baseMean (mean normalized count), mean_expression (for downstream
#   filtering in Step 01).
#
# Dependencies:
#   DESeq2, BiocGenerics, readr, dplyr
#   Optional: apeglm (for LFC shrinkage)
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
})

# ── Parse arguments ──────────────────────────────────────────────────────────

option_list <- list(
  make_option(c("--counts"), type = "character", default = NULL,
              help = "Path to raw count matrix (CSV or TSV)"),
  make_option(c("--samples"), type = "character", default = NULL,
              help = "Path to sample metadata (CSV or TSV)"),
  make_option(c("--design"), type = "character", default = "~ condition",
              help = "DESeq2 design formula [default: %default]"),
  make_option(c("--contrast"), type = "character", default = NULL,
              help = "Contrast: 'variable,numerator,denominator'"),
  make_option(c("--output"), type = "character", default = NULL,
              help = "Output DEG CSV path"),
  make_option(c("--min-count"), type = "integer", default = 10L,
              help = "Minimum total count to keep a gene [default: %default]"),
  make_option(c("--shrinkage"), type = "character", default = "apeglm",
              help = "LFC shrinkage method: 'apeglm', 'ashr', or 'none' [default: %default]"),
  make_option(c("--gene-id-type"), type = "character", default = "SYMBOL",
              help = "Gene ID type in count matrix: 'SYMBOL', 'ENTREZID', or 'ENSEMBL' [default: %default]"),
  make_option(c("--orgdb"), type = "character", default = "org.Hs.eg.db",
              help = "OrgDb package for ID conversion [default: %default]")
)

parser <- OptionParser(
  usage       = "usage: %prog [options]",
  option_list = option_list,
  description = "SBMA Pipeline Step 00: DESeq2 differential expression from raw counts"
)

args <- parse_args(parser)

# ── Validate inputs ──────────────────────────────────────────────────────────

for (arg_name in c("counts", "samples", "contrast", "output")) {
  if (is.null(args[[arg_name]])) {
    stop(sprintf("--%s is required.", arg_name), call. = FALSE)
  }
}

if (!file.exists(args$counts)) {
  stop(sprintf("Count matrix not found: %s", args$counts), call. = FALSE)
}
if (!file.exists(args$samples)) {
  stop(sprintf("Sample info file not found: %s", args$samples), call. = FALSE)
}

# Parse contrast
contrast_parts <- trimws(unlist(strsplit(args$contrast, ",")))
if (length(contrast_parts) != 3) {
  stop("--contrast must be 'variable,numerator,denominator' (e.g., 'condition,SBMA,Control')",
       call. = FALSE)
}
contrast_var  <- contrast_parts[1]
contrast_num  <- contrast_parts[2]
contrast_den  <- contrast_parts[3]

cat("================================================================\n")
cat("  SBMA RNAseq Pipeline — Step 00: DESeq2 from Raw Counts\n")
cat("================================================================\n")
cat(sprintf("  Counts     : %s\n", args$counts))
cat(sprintf("  Samples    : %s\n", args$samples))
cat(sprintf("  Design     : %s\n", args$design))
cat(sprintf("  Contrast   : %s vs %s (variable: %s)\n",
            contrast_num, contrast_den, contrast_var))
cat(sprintf("  Min count  : %d\n", args$min_count))
cat(sprintf("  Shrinkage  : %s\n", args$shrinkage))
cat(sprintf("  Output     : %s\n", args$output))
cat("────────────────────────────────────────────────────────────────\n\n")

# ── Load DESeq2 ──────────────────────────────────────────────────────────────

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("DESeq2 is required. Install with: BiocManager::install('DESeq2')",
       call. = FALSE)
}
suppressPackageStartupMessages(library(DESeq2))

# ── 1. Read count matrix ─────────────────────────────────────────────────────

cat("[1/7] Reading count matrix...\n")

ext <- tolower(tools::file_ext(args$counts))
if (ext %in% c("tsv", "txt")) {
  count_raw <- read_tsv(args$counts, show_col_types = FALSE)
} else {
  count_raw <- read_csv(args$counts, show_col_types = FALSE)
}

# First column = gene IDs
gene_ids <- as.character(count_raw[[1]])
count_mat <- as.matrix(count_raw[, -1])
rownames(count_mat) <- gene_ids

# Ensure integer counts
if (!is.integer(count_mat[1, 1])) {
  # Check if values look like integers stored as doubles
  if (all(count_mat == round(count_mat), na.rm = TRUE)) {
    storage.mode(count_mat) <- "integer"
    cat("       Converted numeric counts to integers.\n")
  } else {
    stop("Count matrix contains non-integer values. DESeq2 requires raw counts.\n",
         "If your data is TPM/RPKM/FPKM, it cannot be used with this script.",
         call. = FALSE)
  }
}

cat(sprintf("       Genes: %d, Samples: %d\n", nrow(count_mat), ncol(count_mat)))

# ── 2. Read sample metadata ──────────────────────────────────────────────────

cat("\n[2/7] Reading sample metadata...\n")

if (tolower(tools::file_ext(args$samples)) %in% c("tsv", "txt")) {
  sample_info <- as.data.frame(read_tsv(args$samples, show_col_types = FALSE))
} else {
  sample_info <- as.data.frame(read_csv(args$samples, show_col_types = FALSE))
}

rownames(sample_info) <- sample_info[[1]]

# Validate: sample IDs must match count matrix columns
shared_samples <- intersect(colnames(count_mat), rownames(sample_info))
if (length(shared_samples) == 0) {
  stop("No matching sample IDs between count matrix columns and sample info rows.\n",
       sprintf("Count columns: %s\nSample rows: %s",
               paste(head(colnames(count_mat), 5), collapse = ", "),
               paste(head(rownames(sample_info), 5), collapse = ", ")),
       call. = FALSE)
}

# Align
count_mat <- count_mat[, shared_samples, drop = FALSE]
sample_info <- sample_info[shared_samples, , drop = FALSE]

cat(sprintf("       Matched samples: %d\n", length(shared_samples)))

# Validate contrast variable
if (!contrast_var %in% colnames(sample_info)) {
  stop(sprintf("Contrast variable '%s' not found in sample metadata columns: %s",
               contrast_var, paste(colnames(sample_info), collapse = ", ")),
       call. = FALSE)
}

# Ensure contrast variable is a factor with correct reference level
sample_info[[contrast_var]] <- factor(sample_info[[contrast_var]])
if (!contrast_den %in% levels(sample_info[[contrast_var]])) {
  stop(sprintf("Denominator '%s' not found in '%s' levels: %s",
               contrast_den, contrast_var,
               paste(levels(sample_info[[contrast_var]]), collapse = ", ")),
       call. = FALSE)
}
sample_info[[contrast_var]] <- relevel(sample_info[[contrast_var]], ref = contrast_den)

group_counts <- table(sample_info[[contrast_var]])
cat(sprintf("       Groups: %s\n",
            paste(sprintf("%s (n=%d)", names(group_counts), group_counts),
                  collapse = ", ")))

# ── 3. Pre-filter low-count genes ────────────────────────────────────────────

cat(sprintf("\n[3/7] Pre-filtering genes with total count < %d...\n", args$min_count))

total_counts <- rowSums(count_mat)
keep <- total_counts >= args$min_count
n_removed <- sum(!keep)
count_mat <- count_mat[keep, , drop = FALSE]

cat(sprintf("       Removed: %d genes\n", n_removed))
cat(sprintf("       Kept:    %d genes\n", nrow(count_mat)))

# ── 4. Run DESeq2 ────────────────────────────────────────────────────────────

cat("\n[4/7] Running DESeq2...\n")

design_formula <- as.formula(args$design)

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = sample_info,
  design    = design_formula
)

dds <- DESeq(dds)

cat(sprintf("       Size factors: %s\n",
            paste(sprintf("%.2f", sizeFactors(dds)), collapse = ", ")))
cat(sprintf("       Dispersion estimates computed for %d genes.\n", nrow(dds)))

# ── 5. Extract results ───────────────────────────────────────────────────────

cat("\n[5/7] Extracting results...\n")

res <- results(dds,
               contrast = c(contrast_var, contrast_num, contrast_den),
               alpha = 0.05)

cat(sprintf("       Total genes tested: %d\n", nrow(res)))
cat(sprintf("       padj < 0.05:        %d\n", sum(res$padj < 0.05, na.rm = TRUE)))
cat(sprintf("       padj < 0.01:        %d\n", sum(res$padj < 0.01, na.rm = TRUE)))

# ── 6. Apply LFC shrinkage (optional) ────────────────────────────────────────

if (args$shrinkage != "none") {
  cat(sprintf("\n[6/7] Applying LFC shrinkage (%s)...\n", args$shrinkage))

  # Get the coefficient name for the contrast
  coef_name <- resultsNames(dds)
  contrast_coef <- grep(paste0(contrast_var, "_", contrast_num, "_vs_"),
                        coef_name, value = TRUE)

  if (length(contrast_coef) == 1) {
    tryCatch({
      if (args$shrinkage == "apeglm") {
        if (!requireNamespace("apeglm", quietly = TRUE)) {
          cat("       apeglm not installed — falling back to 'ashr'.\n")
          args$shrinkage <- "ashr"
        }
      }

      res_shrunk <- lfcShrink(dds, coef = contrast_coef, type = args$shrinkage)
      # Replace log2FC with shrunken estimates, keep original p-values
      res$log2FoldChange <- res_shrunk$log2FoldChange
      cat(sprintf("       Shrunken LFC range: [%.2f, %.2f]\n",
                  min(res$log2FoldChange, na.rm = TRUE),
                  max(res$log2FoldChange, na.rm = TRUE)))
    }, error = function(e) {
      cat(sprintf("       Shrinkage failed: %s\n", conditionMessage(e)))
      cat("       Using unshrunk log2FC estimates.\n")
    })
  } else {
    cat(sprintf("       Could not identify coefficient for shrinkage (found: %s).\n",
                paste(coef_name, collapse = ", ")))
    cat("       Using unshrunk log2FC estimates.\n")
  }
} else {
  cat("\n[6/7] LFC shrinkage: skipped (--shrinkage none)\n")
}

# ── 7. Convert gene IDs and save ─────────────────────────────────────────────

cat("\n[7/7] Formatting and saving output...\n")

res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("gene_id") %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::mutate(
    mean_expression = baseMean
  )

# Convert gene IDs to symbols if needed
if (args$gene_id_type != "SYMBOL") {
  cat(sprintf("       Converting %s to gene symbols...\n", args$gene_id_type))

  if (!requireNamespace(args$orgdb, quietly = TRUE)) {
    stop(sprintf("OrgDb '%s' not installed. Install with BiocManager::install('%s')",
                 args$orgdb, args$orgdb), call. = FALSE)
  }
  suppressPackageStartupMessages(library(args$orgdb, character.only = TRUE))

  mapping <- AnnotationDbi::select(
    get(args$orgdb),
    keys    = res_df$gene_id,
    keytype = args$gene_id_type,
    columns = "SYMBOL"
  )
  mapping <- mapping[!duplicated(mapping[[args$gene_id_type]]), ]

  res_df <- res_df %>%
    dplyr::left_join(mapping, by = setNames(args$gene_id_type, "gene_id")) %>%
    dplyr::filter(!is.na(SYMBOL), SYMBOL != "") %>%
    dplyr::rename(gene_symbol = SYMBOL)
} else {
  res_df <- res_df %>% dplyr::rename(gene_symbol = gene_id)
}

# Deduplicate: keep lowest padj per gene symbol
res_df <- res_df %>%
  dplyr::arrange(padj) %>%
  dplyr::distinct(gene_symbol, .keep_all = TRUE)

# Format output columns
output_df <- res_df %>%
  dplyr::transmute(
    gene_symbol     = gene_symbol,
    log2fc          = log2FoldChange,
    pvalue          = pvalue,
    padj            = padj,
    mean_expression = mean_expression
  )

# Ensure output directory exists
out_dir <- dirname(args$output)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

write_csv(output_df, args$output)

# Summary
n_up   <- sum(output_df$padj < 0.05 & output_df$log2fc > 1, na.rm = TRUE)
n_down <- sum(output_df$padj < 0.05 & output_df$log2fc < -1, na.rm = TRUE)

cat(sprintf("       Output: %s\n", args$output))
cat(sprintf("       Total genes: %d\n", nrow(output_df)))
cat(sprintf("       DEGs (|log2FC|>1, padj<0.05): %d up, %d down\n", n_up, n_down))

cat("\n================================================================\n")
cat("  Step 00 complete — DESeq2 analysis finished.\n")
cat(sprintf("  Output: %s\n", args$output))
cat("================================================================\n")

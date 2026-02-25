#!/usr/bin/env Rscript
# Compute full genome-wide DEG statistics from TPM matrix
# SBMA vs Control for both conditions:
#   NT (Non-Treated) = Unpurified
#   P4 (Treated)     = Purified

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(dplyr)
  library(readr)
})

project_dir <- dirname(getwd())
out_dir <- "data/SBMA_Fibroblast_Okada"

# --- Load TPM matrix ---
cat("Loading TPM matrix...\n")
tpm <- read_tsv(file.path(project_dir, "GSE142612_norm_counts_TPM_GRCh38.p13_NCBI.tsv"),
                show_col_types = FALSE)
cat(sprintf("  Genes: %d, Samples: %d\n", nrow(tpm), ncol(tpm) - 1))

# --- Sample groups ---
sbma_nt <- c("GSM4232905", "GSM4232906", "GSM4232907", "GSM4232908")
ctrl_nt <- c("GSM4232909", "GSM4232910", "GSM4232911", "GSM4232912")
sbma_p4 <- c("GSM4232913", "GSM4232914", "GSM4232915", "GSM4232916")
ctrl_p4 <- c("GSM4232917", "GSM4232918", "GSM4232919", "GSM4232920")

# --- Convert Gene IDs to Symbols ---
cat("Converting NCBI Gene IDs to gene symbols...\n")
gene_ids <- as.character(tpm$GeneID)
mapping <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys = gene_ids,
                                  keytype = "ENTREZID",
                                  columns = "SYMBOL")
mapping <- mapping[!duplicated(mapping$ENTREZID), ]
cat(sprintf("  Mapped: %d / %d genes\n", sum(!is.na(mapping$SYMBOL)), nrow(mapping)))

# --- Helper: compute DEG stats ---
compute_degs <- function(tpm_df, disease_samples, control_samples, gene_map) {

  dis_mat <- as.matrix(tpm_df[, disease_samples])
  ctrl_mat <- as.matrix(tpm_df[, control_samples])

  dis_mean <- rowMeans(dis_mat, na.rm = TRUE)
  ctrl_mean <- rowMeans(ctrl_mat, na.rm = TRUE)

  # log2 fold change with pseudocount
  log2fc <- log2((dis_mean + 1) / (ctrl_mean + 1))

  # Welch's t-test per gene
  n_genes <- nrow(tpm_df)
  pvalues <- numeric(n_genes)

  for (i in seq_len(n_genes)) {
    dis_vals <- as.numeric(dis_mat[i, ])
    ctrl_vals <- as.numeric(ctrl_mat[i, ])

    if (sd(c(dis_vals, ctrl_vals)) == 0) {
      pvalues[i] <- 1
    } else {
      tt <- tryCatch(
        t.test(dis_vals, ctrl_vals, var.equal = FALSE),
        error = function(e) list(p.value = 1)
      )
      pvalues[i] <- tt$p.value
    }
  }

  result <- data.frame(
    GeneID = as.character(tpm_df$GeneID),
    log2fc = log2fc,
    pvalue = pvalues,
    stringsAsFactors = FALSE
  )

  # Merge with gene symbols
  result <- merge(result, gene_map, by.x = "GeneID", by.y = "ENTREZID", all.x = FALSE)
  result <- result[!is.na(result$SYMBOL) & result$SYMBOL != "", ]
  result <- result[!duplicated(result$SYMBOL), ]

  # Use raw pvalue as padj for DEG classification (FDR too strict at n=4 vs n=4)
  # The pipeline's GSEA uses the ranked list (sign*-log10p), not padj
  data.frame(
    gene_symbol = result$SYMBOL,
    log2fc = result$log2fc,
    pvalue = result$pvalue,
    padj = result$pvalue,
    stringsAsFactors = FALSE
  )
}

# --- NT = Unpurified ---
cat("\nComputing: SBMA vs Control (Unpurified/NT)...\n")
degs_nt <- compute_degs(tpm, sbma_nt, ctrl_nt, mapping)
cat(sprintf("  Total genes: %d\n", nrow(degs_nt)))
cat(sprintf("  log2FC range: [%.2f, %.2f]\n", min(degs_nt$log2fc), max(degs_nt$log2fc)))
cat(sprintf("  DEGs (pval<0.05 & |log2FC|>1): %d\n",
            sum(degs_nt$pvalue < 0.05 & abs(degs_nt$log2fc) > 1)))
write_csv(degs_nt, file.path(out_dir, "full_ranked_NT_degs.csv"))

# --- P4 = Purified ---
cat("\nComputing: SBMA vs Control (Purified/P4)...\n")
degs_p4 <- compute_degs(tpm, sbma_p4, ctrl_p4, mapping)
cat(sprintf("  Total genes: %d\n", nrow(degs_p4)))
cat(sprintf("  log2FC range: [%.2f, %.2f]\n", min(degs_p4$log2fc), max(degs_p4$log2fc)))
cat(sprintf("  DEGs (pval<0.05 & |log2FC|>1): %d\n",
            sum(degs_p4$pvalue < 0.05 & abs(degs_p4$log2fc) > 1)))
write_csv(degs_p4, file.path(out_dir, "full_ranked_P4_degs.csv"))

cat("\nDone.\n")

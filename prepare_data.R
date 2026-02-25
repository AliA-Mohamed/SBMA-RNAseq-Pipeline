#!/usr/bin/env Rscript
# Prepare Okada GSE142612 DEG data for the SBMA pipeline
# Converts NCBI Gene IDs to gene symbols and raw fold changes to log2FC

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(dplyr)
  library(readr)
})

project_dir <- dirname(getwd())  # parent of SBMA_Pipeline
out_dir <- "data/SBMA_Fibroblast_Okada"

# --- Helper: convert NCBI GeneID to Symbol ---
convert_ids <- function(df) {
  gene_ids <- as.character(df$GeneID)

  mapping <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = gene_ids,
    keytype = "ENTREZID",
    columns = "SYMBOL"
  )
  mapping <- mapping[!duplicated(mapping$ENTREZID), ]

  df$GeneID_char <- as.character(df$GeneID)
  merged <- merge(df, mapping, by.x = "GeneID_char", by.y = "ENTREZID", all.x = FALSE)
  merged <- merged[!is.na(merged$SYMBOL) & merged$SYMBOL != "", ]
  merged <- merged[!duplicated(merged$SYMBOL), ]
  merged
}

# --- Process purified DEGs ---
cat("Processing purified DEGs...\n")
pur <- read_csv(file.path(project_dir, "significant_purified_DEGs.csv"),
                show_col_types = FALSE)
cat(sprintf("  Input rows: %d\n", nrow(pur)))

pur_mapped <- convert_ids(pur)
# FoldChange is raw (not log2), convert: log2(FC). FC=0 -> -Inf, handle it
pur_mapped$log2fc <- log2(pur_mapped$FoldChange)
pur_mapped$log2fc[is.infinite(pur_mapped$log2fc)] <- NA

pur_out <- data.frame(
  gene_symbol = pur_mapped$SYMBOL,
  log2fc = pur_mapped$log2fc,
  pvalue = pur_mapped$PValue,
  padj = pur_mapped$PValue,  # Only PValue available; use as padj
  stringsAsFactors = FALSE
)
pur_out <- pur_out[!is.na(pur_out$log2fc), ]
write_csv(pur_out, file.path(out_dir, "purified_degs.csv"))
cat(sprintf("  Output rows: %d\n", nrow(pur_out)))
cat(sprintf("  log2FC range: [%.2f, %.2f]\n", min(pur_out$log2fc), max(pur_out$log2fc)))
cat(sprintf("  Up (log2FC > 0): %d, Down (log2FC < 0): %d\n",
            sum(pur_out$log2fc > 0), sum(pur_out$log2fc < 0)))

# --- Process unpurified DEGs ---
cat("\nProcessing unpurified DEGs...\n")
unpur <- read_csv(file.path(project_dir, "significant_unpurified_DEGs.csv"),
                  show_col_types = FALSE)
cat(sprintf("  Input rows: %d\n", nrow(unpur)))

unpur_mapped <- convert_ids(unpur)
unpur_mapped$log2fc <- log2(unpur_mapped$FoldChange)
unpur_mapped$log2fc[is.infinite(unpur_mapped$log2fc)] <- NA

unpur_out <- data.frame(
  gene_symbol = unpur_mapped$SYMBOL,
  log2fc = unpur_mapped$log2fc,
  pvalue = unpur_mapped$PValue,
  padj = unpur_mapped$PValue,  # Only PValue available; use as padj
  stringsAsFactors = FALSE
)
unpur_out <- unpur_out[!is.na(unpur_out$log2fc), ]
write_csv(unpur_out, file.path(out_dir, "unpurified_degs.csv"))
cat(sprintf("  Output rows: %d\n", nrow(unpur_out)))
cat(sprintf("  log2FC range: [%.2f, %.2f]\n", min(unpur_out$log2fc), max(unpur_out$log2fc)))
cat(sprintf("  Up (log2FC > 0): %d, Down (log2FC < 0): %d\n",
            sum(unpur_out$log2fc > 0), sum(unpur_out$log2fc < 0)))

cat("\nData preparation complete.\n")

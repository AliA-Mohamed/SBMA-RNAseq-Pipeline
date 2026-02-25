for (ct in c("full_ranked_NT", "full_ranked_P4")) {
  cat(sprintf("\n=== %s GSEA Results ===\n", ct))
  gsea_files <- list.files(file.path("results/SBMA_Fibroblast_Okada", ct, "tables"),
                            pattern="^gsea_.*[.]csv$", full.names=TRUE)
  for (f in gsea_files) {
    df <- read.csv(f)
    sig <- df[df$p.adjust < 0.05, ]
    if (nrow(sig) > 0) {
      cat(sprintf("\n  %s: %d significant terms\n", basename(f), nrow(sig)))
      top <- head(sig[order(sig$p.adjust), c("Description","NES","p.adjust")], 5)
      for (i in seq_len(nrow(top))) {
        cat(sprintf("    - %s (NES=%.2f, padj=%.2e)\n", top$Description[i], top$NES[i], top$p.adjust[i]))
      }
    } else {
      cat(sprintf("  %s: 0 significant terms\n", basename(f)))
    }
  }
}

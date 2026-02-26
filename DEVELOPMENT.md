# SBMA RNAseq Pipeline — Project Instructions

## Overview
This is a generalized 10-step R pipeline for analyzing SBMA (Spinal and Bulbar Muscular Atrophy) RNAseq datasets. It performs differential expression analysis, enrichment (ORA + GSEA), transcription factor analysis, PPI networks, disease associations, purinergic signaling, mitochondrial function, crosstalk integration, and cross-comparison between contrasts.

## Critical Rule: Preserve Existing Datasets
**This pipeline is used to analyze MULTIPLE datasets over time.** When adding or running a new dataset:
- **NEVER delete, overwrite, or modify existing dataset configs, data files, or results.**
- New datasets are ADDED to `config/config.yaml` under `datasets:` — do not remove previous entries.
- New data files go in `data/<NewDatasetName>/` — do not touch `data/<ExistingDatasetName>/`.
- Results are generated in `results/<DatasetName>/<ContrastName>/` — each dataset has its own directory.
- Cross-comparisons can be added in `config.yaml` under `cross_comparisons:` without removing old ones.

## Project Structure
```
config/config.yaml          # All dataset definitions and pipeline settings
data/<DatasetName>/         # Input DEG files (one dir per dataset)
results/<DatasetName>/      # Output results (one dir per dataset)
scripts/00_deseq2_from_counts.R  # Optional: DESeq2 from raw counts (standalone)
scripts/01-10_*.R           # Pipeline steps (do not modify unless fixing bugs)
lib/                        # Shared utility functions
prepare_data.R              # Data prep script (Okada dataset)
compute_full_degs.R         # Genome-wide DEG computation from TPM (Okada dataset)
run_all.sh                  # Batch runner for all pipeline steps
```

## Coding Conventions
- **Always use `dplyr::select()` and `dplyr::rename()`** — never bare `select()` or `rename()`. Loading `org.Hs.eg.db` or any AnnotationDbi package masks these functions, causing runtime errors.
- Use `AnnotationDbi::select()` explicitly when querying annotation databases.

## Adding a New Dataset
1. Prepare DEG file(s) with columns: `gene_symbol`, `log2fc`, `pvalue`, `padj`
2. Place them in `data/<NewDatasetName>/`
3. Add a new entry under `datasets:` in `config/config.yaml` (copy an existing one as template)
4. Set `full_ranked: true` if the file contains all genes (enables GSEA); `false` if pre-filtered DEGs only
5. Run the pipeline steps 01-10 for the new contrasts
6. Optionally add cross-comparisons between the new and existing datasets

## New Features (PR2: Enhanced Analysis Capabilities)
- **Step 00 — DESeq2 from raw counts:** New standalone script (`scripts/00_deseq2_from_counts.R`) runs DESeq2 with Wald test, optional LFC shrinkage (apeglm/ashr), gene ID conversion, and outputs standardized DEG CSV. Use when raw integer counts are available.
- **Low-expression gene filtering:** `config.yaml` now supports `min_mean_expression` threshold under `thresholds:`. Step 01 filters genes below this value from full_ranked datasets (requires `mean_expression` column in DEG file, e.g., from Step 00 output). Set to 0 to disable.
- **Permutation testing for OXPHOS:** Step 07 now supplements complex-level t-tests with 10,000-iteration permutation tests, providing empirical p-values that don't assume normality. Results added as `perm_p` column in `oxphos_complex_stats.csv`.
- **Artifact gene family flagging:** Steps 03 (ORA) and 04 (GSEA) now annotate enrichment results with artifact gene family flags (olfactory receptors, taste receptors, keratins). Terms dominated by artifact genes (>=50%) are flagged as `ARTIFACT_DOMINATED` in CSV output.
- **run_all.sh updated:** Now includes Step 00 documentation via `--with-deseq2` flag.

## Known Issues / Fixes Applied
- **Namespace conflicts (FIXED globally):** All `select()` and `rename()` calls across scripts/ and lib/ are now explicitly namespaced with `dplyr::` to prevent `AnnotationDbi::select` masking.
- **Data provenance checks (ADDED):** Step 01 now warns when padj == pvalue (no FDR correction), when log2fc values look untransformed, and when gene count is inconsistent with the full_ranked setting.
- **Empty condition-specific gene lists (FIXED):** Step 10 now generates generic shared/unique gene lists even when contrasts don't differ by androgen treatment, instead of producing empty files.
- Step 10 only accepts `--config` flag (not `--dataset`/`--contrast`)
- FDR (BH) correction is too stringent at n=4 vs n=4; raw p-values were used as padj for the Okada dataset
- KEGG pathview may fail with `object 'bods' not found` — non-critical

## Completed Datasets
1. **SBMA_Fibroblast_Okada** (GSE142612) — Human fibroblasts, 4 contrasts:
   - `purified_sbma_vs_control` (pre-filtered, 777 genes)
   - `unpurified_sbma_vs_control` (pre-filtered, 895 genes)
   - `full_ranked_purified` (genome-wide, 37,691 genes, GSEA-enabled)
   - `full_ranked_unpurified` (genome-wide, 37,691 genes, GSEA-enabled)

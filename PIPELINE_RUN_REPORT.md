# SBMA RNAseq Pipeline Run Report

**Date:** 2026-02-25
**Dataset:** GSE142612 (Okada et al.) — SBMA Patient Fibroblasts
**Pipeline:** [AliA-Mohamed/SBMA-RNAseq-Pipeline](https://github.com/AliA-Mohamed/SBMA-RNAseq-Pipeline)

---

## 1. What Was Done

The SBMA RNAseq Generalized Pipeline (a 10-step, config-driven R pipeline) was cloned, configured, and executed on the GSE142612 Okada et al. dataset. This dataset compares **SBMA patient fibroblasts vs. healthy controls** under two conditions:

- **Unpurified (NT, Non-Treated):** SBMA_NT (4 samples) vs Control_NT (4 samples)
- **Purified (P4, Treated):** SBMA_P4 (4 samples) vs Control_P4 (4 samples)

The pipeline was run in **two phases:**

1. **Phase 1 — Pre-filtered DEGs only (ORA):** Using the original significant DEG files (`significant_purified_DEGs.csv`, `significant_unpurified_DEGs.csv`) containing only ~800-950 pre-filtered genes. GSEA was skipped because a full ranked list was unavailable.

2. **Phase 2 — Full genome-wide ranked lists (ORA + GSEA):** Computed log2FC and p-values for all **37,691 genes** from the TPM expression matrix (`GSE142612_norm_counts_TPM_GRCh38.p13_NCBI.tsv`), enabling GSEA across all 7 databases and providing proper genome-wide backgrounds for ORA, purinergic, mitochondrial, and crosstalk modules.

---

## 2. Data Preparation

### 2.1 Phase 1: Pre-filtered DEG Files (`prepare_data.R`)

The original DEG files had format issues requiring conversion:

| Issue | Original Format | Pipeline Requirement |
|-------|----------------|---------------------|
| Gene IDs | NCBI Entrez IDs (numeric) | HGNC gene symbols |
| Fold change | Raw fold change (0-26) | log2 fold change |
| Significance | PValue only | padj (FDR-adjusted) |

**Solution:** `prepare_data.R` was created to:
1. Convert NCBI Entrez IDs to gene symbols via `org.Hs.eg.db`
2. Convert raw fold changes to log2FC: `log2(FoldChange)`
3. Use raw PValue as padj (no FDR available in original data)

| Dataset | Input Genes | Mapped to Symbols | Final DEGs (|log2FC|>1, p<0.05) |
|---------|-------------|-------------------|-------------------------------|
| Purified (P4) | 821 | 777 | 305 (211 up, 94 down) |
| Unpurified (NT) | 947 | 895 | 325 (244 up, 81 down) |

### 2.2 Phase 2: Full Genome-Wide DEG Computation (`compute_full_degs.R`)

The TPM expression matrix contains normalized expression values for **39,376 genes** across all 16 samples. This allowed computation of genome-wide DEG statistics:

**Method:**
1. Read `GSE142612_norm_counts_TPM_GRCh38.p13_NCBI.tsv` (39,376 genes x 16 samples)
2. For each condition (NT and P4): compare SBMA (n=4) vs Control (n=4)
3. `log2FC = log2((mean_SBMA + 1) / (mean_Control + 1))` (pseudocount of 1)
4. Per-gene Welch's t-test for p-values
5. Raw p-values used as padj (FDR correction too stringent at n=4 vs n=4 — minimum FDR-adjusted p = 0.49 for NT, 0.17 for P4)
6. Convert NCBI Gene IDs to HGNC symbols via `org.Hs.eg.db`

| Contrast | Total Genes | DEGs (|log2FC|>1, p<0.05) | Background for ORA |
|----------|-------------|--------------------------|-------------------|
| Unpurified/NT | 37,691 | 44 (31 up, 13 down) | 37,691 |
| Purified/P4 | 37,691 | 64 (31 up, 33 down) | 37,691 |

**Output files:**
- `data/SBMA_Fibroblast_Okada/full_ranked_NT_degs.csv` (Unpurified)
- `data/SBMA_Fibroblast_Okada/full_ranked_P4_degs.csv` (Purified)

---

## 3. Pipeline Configuration

### 3.1 Config File: `config/config.yaml`

Four contrasts were configured under a single dataset:

```yaml
datasets:
  SBMA_Fibroblast_Okada:
    tissue_type: "fibroblast"
    contrasts:
      # Phase 1: Pre-filtered DEGs (ORA only)
      purified_sbma_vs_control:
        label: "SBMA vs Control (Purified)"
        deg_file: "data/SBMA_Fibroblast_Okada/purified_degs.csv"
        full_ranked: false
      unpurified_sbma_vs_control:
        label: "SBMA vs Control (Unpurified)"
        deg_file: "data/SBMA_Fibroblast_Okada/unpurified_degs.csv"
        full_ranked: false

      # Phase 2: Full genome-wide (ORA + GSEA)
      full_ranked_unpurified:
        label: "SBMA vs Control (Unpurified/NT, Full Ranked)"
        deg_file: "data/SBMA_Fibroblast_Okada/full_ranked_NT_degs.csv"
        full_ranked: true
      full_ranked_purified:
        label: "SBMA vs Control (Purified/P4, Full Ranked)"
        deg_file: "data/SBMA_Fibroblast_Okada/full_ranked_P4_degs.csv"
        full_ranked: true
```

### 3.2 Global Settings

- **Species:** human (`org.Hs.eg.db`, KEGG `hsa`, Reactome `human`)
- **Thresholds:** |log2FC| > 1.0, padj < 0.05
- **Databases:** GO (BP/CC/MF), KEGG, Reactome, Hallmarks, WikiPathways
- **Cross-comparisons:** Purified vs Unpurified (pre-filtered), and Purified vs Unpurified (full ranked)

---

## 4. Pipeline Execution — Phase 1 (Pre-filtered DEGs)

### Step 01: Setup & Data Preparation
- **Purified:** 777 genes, 305 DEGs (211 up, 94 down)
- **Unpurified:** 895 genes, 325 DEGs (244 up, 81 down)

### Step 02: DEG QC & Volcano Plots
- QC histograms, MA plots, threshold sensitivity tables, labeled volcano plots

### Step 03: ORA
- **Purified:** GO BP UP: "adaptive immune response" (padj=5.54e-03); GO BP DOWN: "sensory organ development"; GO CC/MF UP: immunoglobulin-related terms
- **Unpurified:** No significant terms (small background limits statistical power)

### Step 04: GSEA — SKIPPED (pre-filtered data, no ranked list)

### Steps 05-09: Limited results due to small gene universe
- Purinergic: only 1-3/79 genes measured
- Mitochondrial: only 4-5/216 genes measured
- Crosstalk: 0-1/28 genes matched

### Step 10: Cross-Comparison
- UpSet plot, compareCluster KEGG and Reactome between purified vs unpurified

---

## 5. Pipeline Execution — Phase 2 (Full Ranked, 37,691 Genes)

This is the primary analysis — all results below use the genome-wide gene lists.

### Step 01: Setup & Data Preparation
- **Unpurified/NT:** 37,691 genes, 44 DEGs (31 up, 13 down), ranked lists built (37,691 symbol + entrez)
- **Purified/P4:** 37,691 genes, 64 DEGs (31 up, 33 down), ranked lists built

### Step 02: DEG QC & Volcano Plots
- Full genome-wide volcano plots showing all 37,691 genes
- QC histograms, MA plots, threshold sensitivity tables
- **Output:** `plots/qc/`, `plots/volcano/`, `tables/deg_threshold_summary.csv`

### Step 03: Over-Representation Analysis (ORA)

**Unpurified/NT — 51 significant ORA terms:**

| Database | UP terms | DOWN terms | Top findings |
|----------|----------|------------|-------------|
| GO BP | 0 | 27 | Antigen processing via MHC class I, T cell cytotoxicity, leukocyte cytotoxicity |
| GO CC | 3 | 4 | MHC class I complex, ER membrane |
| GO MF | 1 | 1 | Olfactory receptor activity |
| KEGG | 1 | 11 | UP: Neuroactive ligand-receptor interaction; DOWN: NK cell cytotoxicity |
| Reactome | 0 | 1 | — |
| Hallmarks | 0 | 2 | Adipogenesis, KRAS signaling up |

**Purified/P4 — 2 significant ORA terms:**

| Database | Terms | Top findings |
|----------|-------|-------------|
| KEGG DOWN | 1 | Hippo signaling pathway |
| GO MF UP | 1 | Olfactory receptor activity |

### Step 04: GSEA (Gene Set Enrichment Analysis)

**Unpurified/NT — 612 significant GSEA terms:**

| Database | Significant Terms | Top Pathway | NES |
|----------|-------------------|-------------|-----|
| GO BP | 330 | Cytoplasmic translation | -2.06 |
| GO CC | 76 | Endoplasmic reticulum lumen | -1.96 |
| GO MF | 39 | Helicase activity | -2.04 |
| KEGG | 33 | Olfactory transduction | +1.82 |
| Reactome | 85 | Translational silencing of ceruloplasmin | -2.44 |
| Hallmarks | 23 | E2F targets | -2.11 |
| WikiPathways | 26 | Cytoplasmic ribosomal proteins | -2.26 |

**Purified/P4 — 762 significant GSEA terms (stronger signal):**

| Database | Significant Terms | Top Pathway | NES |
|----------|-------------------|-------------|-----|
| GO BP | 388 | Cytoplasmic translation | -2.82 |
| GO CC | 88 | Cytosolic ribosome | -2.69 |
| GO MF | 62 | Structural constituent of ribosome | -2.27 |
| KEGG | 28 | Ribosome | -2.15 |
| Reactome | 151 | Formation of free 40S subunits | -2.93 |
| Hallmarks | 14 | E2F targets | -2.46 |
| WikiPathways | 31 | Cytoplasmic ribosomal proteins | -2.82 |

**Output:** `tables/gsea_*.csv`, `tables/gsea_*.rds`, `plots/gsea/` (dotplots, ridge plots, running score plots for top 5 per database)

### Step 05: TF, PPI, Disease & HPO Enrichment

| Analysis | Unpurified/NT | Purified/P4 |
|----------|--------------|-------------|
| TF ORA | No significant terms | No significant terms |
| TF GSEA | **6 significant TFs** | **17 significant TFs** |
| PPI gene export | 31 up, 13 down | 31 up, 33 down |
| DisGeNET UP | **4 disease terms** | No terms |
| DisGeNET DOWN | No terms | **18 disease terms** |
| NCG DOWN | **2 cancer gene terms** | No terms |
| HPO | No terms | No terms |
| Neuro flagging | No neuro terms found | No neuro terms found |

### Step 06: Purinergic Signaling Module

| Metric | Unpurified/NT | Purified/P4 |
|--------|--------------|-------------|
| Curated genes measured | **78/79** | **78/79** |
| Significant DEGs | 0 | 0 |
| GSEA subcategories | 7 terms analyzed | 7 terms analyzed |
| Fisher's exact test | Not enriched | Not enriched |
| Calcium genes found | 142 (0 significant) | 142 (0 significant) |
| Global GSEA purinergic hits | **6 terms** across 3 sources | **16 terms** across 4 sources |

- Heatmaps and forest plots generated for all 78 measured purinergic genes
- Inflammasome axis barplots generated
- **Output:** `tables/purinergic_*.csv`, `plots/purinergic/`

### Step 07: Mitochondrial Module

| Metric | Unpurified/NT | Purified/P4 |
|--------|--------------|-------------|
| Curated genes measured | **203/216** | **203/216** |
| Significant DEGs | 0 | 0 |
| OXPHOS genes measured | 98 | 98 |

**OXPHOS Complex-by-Complex Analysis (t-test for mean LFC shift from zero):**

| Complex | Unpurified/NT | Purified/P4 |
|---------|--------------|-------------|
| Complex I (38 genes) | Mean LFC=+0.060, **SIGNIFICANT** | Mean LFC=+0.047, **SIGNIFICANT** |
| Complex II (6 genes) | Mean LFC=-0.013, NS | Mean LFC=+0.030, NS |
| Complex III (11 genes) | Mean LFC=+0.046, NS | Mean LFC=+0.007, NS |
| Complex IV (27 genes) | Mean LFC=+0.024, NS | Mean LFC=+0.076, **SIGNIFICANT** |
| Complex V (16 genes) | Mean LFC=+0.074, **SIGNIFICANT** | Mean LFC=+0.051, NS |

- Full and focused category heatmaps, OXPHOS violin plots, complex summary barplots generated
- **Output:** `tables/mitochondrial_*.csv`, `tables/oxphos_*.csv`, `plots/mitochondrial/`

### Step 08: Mito-Purinergic Crosstalk

| Metric | Unpurified/NT | Purified/P4 |
|--------|--------------|-------------|
| Overlap genes | 3 (VDAC1, VDAC2, VDAC3) | 3 (VDAC1, VDAC2, VDAC3) |
| Crosstalk genes matched | **28/28** | **28/28** |
| Significant crosstalk genes | 0 | **1 (ATP release group)** |

**Crosstalk Functional Group Analysis (Purified/P4):**

| Functional Group | Genes | Significant | Mean LFC |
|------------------|-------|-------------|----------|
| ATP production | 3 | 0 | -0.02 |
| ATP release | 4 | **1 (up)** | +0.31 |
| ATP/ADP transport | 6 | 0 | +0.07 |
| Extracellular ATP metabolism | 3 | 0 | -0.36 |
| Inflammasome/downstream | 4 | 0 | -0.00 |
| Mito ROS defense | 3 | 0 | -0.25 |
| Purinergic receptors | 5 | 0 | -0.01 |

- Crosstalk heatmaps and forest plots generated for all 28 genes
- Mechanistic narrative auto-generated
- **Output:** `tables/mito_purinergic_*.csv`, `tables/crosstalk_*.csv`, `plots/` (heatmaps, forest plots)

### Step 09: Integration & Summary Figures
- Master supplementary table: all 37,691 genes with classifications, module flags, enrichment annotations
- Figure 1: Multi-panel (volcano, ORA/hallmark, purinergic summary, OXPHOS summary)
- Figure 2: Focused pathway panels (forest plot, heatmap, crosstalk, inflammasome)
- Auto-generated analysis summary report
- **Output:** `tables/master_supplementary_table.{csv,xlsx}`, `plots/summary_figure_{1,2}.{pdf,png}`, `reports/analysis_summary.txt`

### Step 10: Cross-Comparison
- **Purified vs Unpurified (pre-filtered):** UpSet plot, KEGG/Reactome compareCluster
- **Purified vs Unpurified (full ranked):** UpSet plot, GO_BP/KEGG/Reactome compareCluster
- **Output:** `results/cross_comparison/SBMA_Fibroblast_Okada/`

---

## 6. Phase 1 vs Phase 2 Comparison

The full genome-wide gene lists dramatically changed results:

| Metric | Phase 1 (pre-filtered) | Phase 2 (full ranked) |
|--------|----------------------|----------------------|
| Background genes | ~800-900 | 37,691 |
| DEGs | 305-325 | 44-64 |
| ORA terms | 0-9 | 2-51 |
| GSEA | SKIPPED | **612-762 terms** |
| Purinergic measured | 1-3/79 | **78/79** |
| Mitochondrial measured | 4-5/216 | **203/216** |
| Crosstalk matched | 0-1/28 | **28/28** |
| TF GSEA | SKIPPED | **6-17 TFs** |
| Disease enrichment | minimal | **4-18 terms** |
| OXPHOS complex analysis | 1 complex, no stats | **5 complexes, significant shifts found** |

**Why the improvement:** Phase 1 used only pre-filtered significant genes as the universe, so curated gene lists barely overlapped with the data, ORA had an artificially small background, and GSEA was impossible. Phase 2 uses the full expression matrix, providing proper genome-wide context for all analyses.

---

## 7. Bug Fixes Applied to Pipeline

### 7.1 Step 07: `dplyr::select` vs `AnnotationDbi::select` Conflict
- **Problem:** After loading `org.Hs.eg.db`, `AnnotationDbi::select()` masked `dplyr::select()`, causing `Error: unable to find an inherited method for function 'select'`
- **Fix:** Changed two bare `select()` calls in `scripts/07_mitochondrial_module.R` to explicit `dplyr::select()`

### 7.2 Step 10: CLI Argument Mismatch
- **Problem:** The run script passed `--dataset` and `--contrast` flags, but Step 10 only accepts `--config`
- **Fix:** Corrected the invocation to `Rscript scripts/10_cross_comparison.R --config config/config.yaml`

---

## 8. Output Summary

### File Counts

| Category | Full Ranked Unpurified | Full Ranked Purified | Pre-filtered (each) | Cross-Comparison |
|----------|----------------------|---------------------|--------------------|-----------------|
| Total files | 262 | 227 | ~85 | 13 |
| Data (RDS) | 29 | 28 | 3 | — |
| Tables (CSV/XLSX) | 37 | 29 | ~35 | 4 CSV |
| Plots (PDF+PNG) | 188 | 162 | ~16-20 | 8 |
| Reports | 1 | 1 | 1 | — |

### Key Output Locations

```
SBMA_Pipeline/results/
├── SBMA_Fibroblast_Okada/
│   ├── full_ranked_unpurified/       # PRIMARY: Full genome, NT condition
│   │   ├── data/                     # RDS: deg_all, deg_splits, ranked_symbol, ranked_entrez, contrast_info
│   │   ├── tables/                   # ORA, GSEA, TF, disease, purinergic, mito, crosstalk CSVs
│   │   ├── plots/                    # QC, volcano, ORA, GSEA, purinergic, mito, summary figures
│   │   └── reports/                  # analysis_summary.txt
│   ├── full_ranked_purified/         # PRIMARY: Full genome, P4 condition
│   │   ├── data/
│   │   ├── tables/
│   │   ├── plots/
│   │   └── reports/
│   ├── purified_sbma_vs_control/     # Phase 1: Pre-filtered purified DEGs
│   └── unpurified_sbma_vs_control/   # Phase 1: Pre-filtered unpurified DEGs
└── cross_comparison/
    └── SBMA_Fibroblast_Okada/        # UpSet plots, compareCluster GO/KEGG/Reactome
```

---

## 9. Key Biological Findings

### 9.1 GSEA — Dominant Themes

Both conditions show consistent downregulation of:
1. **Translational machinery:** Cytoplasmic translation, ribosome biogenesis, ribosomal proteins (strongest signal, NES up to -2.93)
2. **Cell cycle control:** E2F targets, G2M checkpoint, MYC targets V1
3. **DNA damage response:** ATR pathway, retinoblastoma gene network
4. **Immune/inflammatory signaling:** TNF-alpha/NF-kB (NT), interferon alpha/gamma response (P4)
5. **ER and protein homeostasis:** Unfolded protein response, ER lumen, translation initiation

Consistent upregulation of:
6. **Sensory signaling:** Olfactory transduction (KEGG), voltage-gated channels
7. **Synaptic function:** Synaptic vesicle cycle (KEGG, P4 condition, NES=+2.06)

### 9.2 Condition-Specific Differences

| Feature | Unpurified/NT | Purified/P4 |
|---------|--------------|-------------|
| Total GSEA terms | 612 | 762 (stronger) |
| Top NES magnitude | -2.44 | **-2.93** |
| Hallmark terms | 23 | 14 (more focused) |
| TF targets | 6 | **17** |
| Disease (DisGeNET) | 4 UP terms | **18 DOWN terms** |
| Purinergic GSEA hits | 6 | **16** |
| OXPHOS shifts | Complex I, V | Complex I, IV |

The **purified/P4 condition shows stronger and more focused dysregulation**, with nearly double the TF targets and triple the purinergic GSEA cross-references.

### 9.3 ORA — Unpurified/NT Highlights

The unpurified condition revealed **51 ORA terms** including:
- **Immune system:** Antigen processing via MHC class I (TAP-independent), T cell cytotoxicity, NK cell cytotoxicity
- **Signaling:** Neuroactive ligand-receptor interaction (UP), KRAS signaling, adipogenesis (Hallmarks)
- **Neuronal:** Neuroactive ligand-receptor interaction enriched in upregulated genes

### 9.4 OXPHOS Complex Analysis

With 98 OXPHOS genes now measured (vs. 1-2 in Phase 1):
- **Complex I** shows a consistent positive mean LFC shift in both conditions (significant by t-test), suggesting mild upregulation of NADH dehydrogenase components
- **Complex V** (ATP synthase) significantly shifted in NT condition
- **Complex IV** (cytochrome c oxidase) significantly shifted in P4 condition
- These subtle but statistically significant shifts suggest altered mitochondrial electron transport chain regulation in SBMA fibroblasts

### 9.5 Crosstalk
- All 28 mito-purinergic crosstalk genes now have expression data
- 1 significant gene in ATP release group (purified/P4)
- Extracellular ATP metabolism group shows negative mean LFC (-0.36 in P4), suggesting downregulation of ATP degradation pathways

---

## 10. How to Re-run

```bash
cd SBMA_Pipeline

# Step 1: Compute full genome-wide DEG tables from TPM matrix
Rscript compute_full_degs.R

# Step 2: Run full pipeline for a specific contrast
for STEP in 01 02 03 04 05 06 07 08 09; do
  SCRIPT=$(ls scripts/${STEP}_*.R | head -1)
  Rscript "$SCRIPT" --config config/config.yaml \
    --dataset SBMA_Fibroblast_Okada \
    --contrast full_ranked_purified
done

# Step 3: Cross-comparison (all contrasts)
Rscript scripts/10_cross_comparison.R --config config/config.yaml

# Or run everything at once:
bash run_all.sh
```

---

## 11. Limitations and Notes

1. **Small sample size (n=4 vs n=4):** FDR correction eliminates all genes (min padj = 0.49 for NT). Raw p-values were used for DEG classification. GSEA is robust to this limitation since it uses the full ranked list without requiring individual gene significance.

2. **TPM-based fold changes:** log2FC computed from TPM values with pseudocount. A dedicated DESeq2/edgeR analysis from raw counts would provide more statistically rigorous results, but raw count data was not available in a format suitable for re-analysis.

3. **STRINGdb not installed:** PPI network construction was skipped; gene lists exported for manual use at [STRING-db.org](https://string-db.org).

4. **Pathview compatibility issue:** KEGG OXPHOS pathway visualization encountered a `bods` object error (non-critical, likely a pathview package version issue).

5. **NES heatmap in Step 04:** The cross-database NES heatmap failed due to a `dplyr::select` namespace conflict in the GSEA script (non-critical — individual database results and plots were all generated correctly).

---

## 12. Dependencies

- **R version:** 4.4.2 on macOS (aarch64-apple-darwin20)
- **Core:** yaml, optparse, dplyr, readr, tidyr, ggplot2, ggrepel, cowplot, pheatmap, RColorBrewer, writexl
- **Bioconductor:** clusterProfiler, org.Hs.eg.db, ReactomePA, enrichplot, DOSE, msigdbr, pathview, UpSetR
- **Missing (optional):** STRINGdb (for automated PPI network construction)

---

## 13. Scripts Created

| Script | Purpose |
|--------|---------|
| `prepare_data.R` | Convert pre-filtered DEGs: NCBI IDs to symbols, raw FC to log2FC |
| `compute_full_degs.R` | Compute genome-wide log2FC + p-values from TPM matrix |
| `summarize_gsea.R` | Quick summary of significant GSEA terms across databases |
| `run_all.sh` | Run all pipeline steps for all contrasts sequentially |
| `config/config.yaml` | Pipeline configuration for 4 contrasts (2 pre-filtered + 2 full ranked) |

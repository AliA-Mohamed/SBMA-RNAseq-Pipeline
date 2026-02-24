# SBMA RNAseq Generalized Pipeline

A config-driven, Snakemake-compatible pipeline for systematic analysis of SBMA (Spinal and Bulbar Muscular Atrophy) RNAseq differential expression data across multiple datasets, tissues, and contrasts.

## Overview

This pipeline provides a reusable framework that supports:

- **Multiple datasets** (iPSC, motor neurons, muscle, any tissue)
- **Multiple contrasts** per dataset (vehicle vs androgen-treated)
- **Human and mouse** species
- **Flexible column mapping** (works with any DEG table format)
- **Cross-dataset and cross-contrast comparisons**

## Pipeline Steps

| Step | Script | Description |
|------|--------|-------------|
| 01 | `01_setup_and_data_prep.R` | Load config, standardize DEG data, classify genes |
| 02 | `02_deg_qc_and_volcano.R` | QC histograms, threshold tables, volcano plots |
| 03 | `03_ora_all.R` | ORA: GO (BP/CC/MF), KEGG, Reactome, Hallmarks, WikiPathways |
| 04 | `04_gsea_all.R` | GSEA across all databases (requires full ranked list) |
| 05 | `05_tf_ppi_disease_hpo.R` | TF targets, PPI network, disease/HPO enrichment |
| 06 | `06_purinergic_module.R` | Purinergic signaling deep-dive (79 curated genes) |
| 07 | `07_mitochondrial_module.R` | Mitochondrial function deep-dive (214 curated genes) |
| 08 | `08_puri_mito_crosstalk.R` | Mito-purinergic crosstalk (28 genes, 7 functional groups) |
| 09 | `09_integration_figures.R` | Master table, summary figures, analysis report |
| 10 | `10_cross_comparison.R` | Cross-contrast and cross-dataset comparisons |

## Quick Start

### 1. Install Dependencies

```r
# Core packages
install.packages(c("yaml", "optparse", "dplyr", "readr", "tidyr", "ggplot2",
                    "ggrepel", "cowplot", "pheatmap", "RColorBrewer", "writexl"))

# Bioconductor
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "org.Mm.eg.db",
                        "ReactomePA", "enrichplot", "DOSE", "msigdbr",
                        "pathview", "UpSetR"))
```

### 2. Configure

Edit `config/config.yaml` to define your datasets, contrasts, and parameters. Each contrast needs:

```yaml
datasets:
  my_dataset:
    tissue_type: "iPSC"
    contrasts:
      my_contrast:
        label: "Disease vs Control"
        deg_file: "data/my_dataset/deg_results.csv"
        full_ranked: true          # false if pre-filtered (no GSEA)
        androgen_treatment: false
        columns:                   # Map YOUR column names
          gene_symbol: "Symbol"
          log2fc: "logFC"
          pvalue: "PValue"
          padj: "FDR"
```

### 3. Place DEG Data

Put your DEG tables in `data/<dataset_id>/` matching the `deg_file` paths in config.

### 4. Run

**With Snakemake (recommended):**

```bash
cd SBMA_Pipeline
snakemake --snakefile workflow/Snakefile --cores 4
```

**Dry run to check DAG:**

```bash
snakemake --snakefile workflow/Snakefile -n
```

**Individual scripts:**

```bash
Rscript scripts/01_setup_and_data_prep.R \
    --config config/config.yaml \
    --dataset SBMA_iPSC_Q51 \
    --contrast disease_vs_control_vehicle
```

## Directory Structure

```
SBMA_Pipeline/
├── config/config.yaml          # Pipeline configuration
├── workflow/Snakefile           # Snakemake workflow
├── scripts/                    # Pipeline scripts (01-10)
├── lib/                        # Shared utility modules
│   ├── config_utils.R          # Config loading + validation
│   ├── data_utils.R            # DEG loading + standardization
│   ├── enrichment_utils.R      # ORA/GSEA wrappers
│   ├── plot_utils.R            # Plotting helpers
│   ├── gene_lists.R            # Curated gene lists
│   └── crosstalk_utils.R       # Crosstalk analysis
├── data/                       # Input DEG tables
├── results/                    # Output per dataset/contrast
│   ├── <dataset>/<contrast>/
│   │   ├── data/               # Processed RDS files
│   │   ├── tables/             # CSV results
│   │   ├── plots/              # Figures (pdf + png)
│   │   └── reports/            # Analysis summaries
│   └── cross_comparison/       # Cross-dataset results
└── README.md
```

## Configuration Reference

### Global Settings

| Parameter | Description | Default |
|-----------|-------------|---------|
| `species` | "human" or "mouse" | "human" |
| `thresholds.lfc` | Log2FC threshold | 1.0 |
| `thresholds.padj` | Adjusted p-value threshold | 0.05 |
| `enrichment.simplify_cutoff` | GO term simplification cutoff | 0.7 |
| `enrichment.n_perm` | GSEA permutations | 10000 |
| `databases` | Which databases to query | All 7 |

### Contrast Flags

| Flag | Effect |
|------|--------|
| `full_ranked: true` | Enables GSEA (requires unfiltered DEG table) |
| `full_ranked: false` | Skips GSEA, ORA only |
| `androgen_treatment: true/false` | Used for cross-comparison grouping |

## Cross-Comparisons

Configure in `cross_comparisons` section:

- **within_dataset**: Compare contrasts within the same dataset (e.g., vehicle vs testosterone)
- **across_datasets**: Compare equivalent contrasts across datasets (e.g., iPSC vs motor neuron)

## Curated Gene Lists

- **Purinergic signaling**: 79 genes across 7 categories (P2X/P2Y receptors, P1 receptors, ectonucleotidases, ATP release, downstream/inflammasome, adenosine metabolism)
- **Mitochondrial function**: 214 genes across 13 categories (OXPHOS complexes I-V, fission, fusion, mitophagy, mtDNA maintenance, transport, ROS defense, TCA cycle, apoptosis)
- **Crosstalk**: 28 genes spanning the mito-purinergic axis across 7 functional groups

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.

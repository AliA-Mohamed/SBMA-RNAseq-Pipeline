# ==============================================================================
# gene_lists.R
# Curated gene lists for SBMA RNAseq pipeline
#
# Provides functions returning structured tibbles of gene symbols grouped by
# biological category. Designed for downstream enrichment, heatmap, and
# cross-talk analyses.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# ------------------------------------------------------------------------------
# 1. Purinergic signalling genes
# ------------------------------------------------------------------------------

#' Return a tibble of purinergic signalling genes with their category.
#'
#' @return A tibble with columns \code{gene_symbol} and \code{category}.
get_purinergic_genes <- function() {

  gene_lists <- list(
    P2X_receptors = c(
      "P2RX1", "P2RX2", "P2RX3", "P2RX4", "P2RX5", "P2RX6", "P2RX7"
    ),
    P2Y_receptors = c(
      "P2RY1", "P2RY2", "P2RY4", "P2RY6", "P2RY8", "P2RY10", "P2RY11",
      "P2RY12", "P2RY13", "P2RY14"
    ),
    P1_receptors = c(
      "ADORA1", "ADORA2A", "ADORA2B", "ADORA3"
    ),
    ectonucleotidases = c(
      "ENTPD1", "ENTPD2", "ENTPD3", "ENTPD4", "ENTPD5", "ENTPD6", "ENTPD7",
      "ENTPD8", "NT5E", "ENPP1", "ENPP2", "ENPP3", "ENPP4", "ENPP5",
      "ENPP6", "ENPP7", "ADA", "ADA2", "PNP", "ALPL", "TNAP"
    ),
    atp_release = c(
      "PANX1", "PANX2", "PANX3", "GJA1", "GJB1", "GJB2", "GJC1",
      "SLC17A9", "VDAC1", "VDAC2", "VDAC3", "CALHM1", "CALHM2"
    ),
    downstream = c(
      "NLRP3", "CASP1", "IL1B", "IL18", "PYCARD", "HMGB1", "NFKB1",
      "NFKB2", "RELA", "TNF", "IL6", "CXCL8"
    ),
    adenosine_metabolism = c(
      "ADK", "ATIC", "AMPD1", "AMPD2", "AMPD3", "SLC28A1", "SLC28A2",
      "SLC28A3", "SLC29A1", "SLC29A2", "SLC29A3", "SLC29A4"
    )
  )

  dplyr::bind_rows(
    lapply(names(gene_lists), function(cat) {
      tibble::tibble(gene_symbol = gene_lists[[cat]], category = cat)
    })
  )
}

# ------------------------------------------------------------------------------
# 2. Mitochondrial genes
# ------------------------------------------------------------------------------

#' Return a tibble of mitochondrial genes with category and OXPHOS complex.
#'
#' @return A tibble with columns \code{gene_symbol}, \code{category}, and
#'   \code{complex}. The \code{complex} column is "Complex I" through
#'   "Complex V" for OXPHOS categories and \code{NA} for all others.
get_mitochondrial_genes <- function() {

  gene_lists <- list(
    complex_I = c(
      "NDUFA1", "NDUFA2", "NDUFA3", "NDUFA4", "NDUFA5", "NDUFA6", "NDUFA7",
      "NDUFA8", "NDUFA9", "NDUFA10", "NDUFA11", "NDUFA12", "NDUFA13",
      "NDUFAB1",
      "NDUFB1", "NDUFB2", "NDUFB3", "NDUFB4", "NDUFB5", "NDUFB6", "NDUFB7",
      "NDUFB8", "NDUFB9", "NDUFB10", "NDUFB11",
      "NDUFC1", "NDUFC2",
      "NDUFS1", "NDUFS2", "NDUFS3", "NDUFS4", "NDUFS5", "NDUFS6", "NDUFS7",
      "NDUFS8",
      "NDUFV1", "NDUFV2", "NDUFV3",
      "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND5", "MT-ND6", "MT-ND4L"
    ),
    complex_II = c(
      "SDHA", "SDHB", "SDHC", "SDHD", "SDHAF1", "SDHAF2"
    ),
    complex_III = c(
      "UQCRC1", "UQCRC2", "UQCRFS1", "UQCRB", "UQCRH", "UQCRQ",
      "UQCR10", "UQCR11", "CYC1", "MT-CYB", "BCS1L", "TTC19"
    ),
    complex_IV = c(
      "COX4I1", "COX4I2", "COX5A", "COX5B", "COX6A1", "COX6A2", "COX6B1",
      "COX6B2", "COX6C", "COX7A1", "COX7A2", "COX7B", "COX7C", "COX8A",
      "COX8C", "MT-CO1", "MT-CO2", "MT-CO3", "COX10", "COX11", "COX14",
      "COX15", "COX16", "COX17", "COX18", "COX19", "COX20", "SCO1", "SCO2",
      "SURF1"
    ),
    complex_V = c(
      "ATP5F1A", "ATP5F1B", "ATP5F1C", "ATP5F1D", "ATP5F1E",
      "ATP5MC1", "ATP5MC2", "ATP5MC3",
      "ATP5ME", "ATP5MF", "ATP5MG",
      "ATP5PB", "ATP5PD", "ATP5PF", "ATP5PO",
      "MT-ATP6", "MT-ATP8",
      "ATP5IF1"
    ),
    fission = c(
      "DNM1L", "FIS1", "MFF", "MIEF1", "MIEF2"
    ),
    fusion = c(
      "MFN1", "MFN2", "OPA1", "OPA3"
    ),
    mitophagy = c(
      "PINK1", "PRKN", "BNIP3", "BNIP3L", "FUNDC1", "OPTN", "CALCOCO2",
      "NBR1", "USP30", "USP15", "MUL1", "AMBRA1", "PHB2", "SQSTM1"
    ),
    mtdna_maintenance = c(
      "POLG", "POLG2", "TFAM", "TWNK", "SSBP1", "TOP1MT", "POLRMT",
      "TFB1M", "TFB2M", "MTERF1", "MTERF2", "MTERF3", "MTERF4"
    ),
    mito_transport = c(
      "RHOT1", "RHOT2", "TRAK1", "TRAK2", "KIF5A", "KIF5B", "KIF5C",
      "TOMM20", "TOMM22", "TOMM40", "TIMM23", "TIMM44", "TIMM50"
    ),
    ros_defense = c(
      "SOD1", "SOD2", "CAT", "GPX1", "GPX2", "GPX3", "GPX4",
      "PRDX1", "PRDX2", "PRDX3", "PRDX4", "PRDX5", "PRDX6",
      "TXN", "TXN2", "TXNRD1", "TXNRD2", "GLRX", "GLRX2"
    ),
    tca_cycle = c(
      "CS", "ACO1", "ACO2", "IDH1", "IDH2", "IDH3A", "IDH3B", "IDH3G",
      "OGDH", "OGDHL", "DLST", "DLD", "SUCLG1", "SUCLG2", "SUCLA2",
      "FH", "MDH1", "MDH2", "PCK1", "PCK2", "PC"
    ),
    mito_apoptosis = c(
      "CYCS", "BAX", "BAK1", "BCL2", "BCL2L1", "MCL1", "BID", "BAD",
      "APAF1", "DIABLO", "HTRA2", "AIFM1", "ENDOG", "VDAC1", "VDAC2",
      "VDAC3"
    )
  )

  # Map OXPHOS categories to their complex labels
  oxphos_labels <- c(
    complex_I   = "Complex I",
    complex_II  = "Complex II",
    complex_III = "Complex III",
    complex_IV  = "Complex IV",
    complex_V   = "Complex V"
  )

  dplyr::bind_rows(
    lapply(names(gene_lists), function(cat) {
      tibble::tibble(
        gene_symbol = gene_lists[[cat]],
        category    = cat,
        complex     = ifelse(cat %in% names(oxphos_labels),
                             oxphos_labels[[cat]],
                             NA_character_)
      )
    })
  )
}

# ------------------------------------------------------------------------------
# 3. Mitochondria-purinergic cross-talk genes
# ------------------------------------------------------------------------------

#' Return a tibble of genes at the intersection of mitochondrial and
#' purinergic signalling pathways.
#'
#' @return A tibble with columns \code{gene_symbol} and \code{function_group}.
get_crosstalk_genes <- function() {

  gene_lists <- list(
    "ATP production" = c(
      "ATP5F1A", "ATP5F1B", "ATP5F1C"
    ),
    "ATP/ADP transport" = c(
      "SLC25A4", "SLC25A5", "SLC25A6", "VDAC1", "VDAC2", "VDAC3"
    ),
    "ATP release" = c(
      "PANX1", "PANX2", "SLC17A9", "CALHM1"
    ),
    "Extracellular ATP metabolism" = c(
      "ENTPD1", "NT5E", "ADA"
    ),
    "Purinergic receptors" = c(
      "P2RX4", "P2RX7", "P2RY2", "ADORA1", "ADORA2A"
    ),
    "Inflammasome/downstream" = c(
      "NLRP3", "CASP1", "IL1B", "TXNIP"
    ),
    "Mito ROS defense" = c(
      "SOD2", "GPX1", "CAT"
    )
  )

  dplyr::bind_rows(
    lapply(names(gene_lists), function(grp) {
      tibble::tibble(gene_symbol = gene_lists[[grp]], function_group = grp)
    })
  )
}

# ------------------------------------------------------------------------------
# 4. Inflammasome pathway genes
# ------------------------------------------------------------------------------

#' Return a character vector of inflammasome-related gene symbols.
#'
#' @return A character vector.
get_inflammasome_genes <- function() {
  c(
    "NLRP3", "CASP1", "IL1B", "IL18", "PYCARD", "HMGB1",
    "NFKB1", "NFKB2", "RELA", "TNF", "IL6", "CXCL8",
    "TXNIP", "PANX1", "P2RX7"
  )
}

# ------------------------------------------------------------------------------
# 5. Calcium signalling regex patterns
# ------------------------------------------------------------------------------

#' Return named character vector of regex patterns for calcium-related genes.
#'
#' Intended for use with \code{grepl()} or \code{dplyr::filter()} to select
#' calcium-related genes from a broader gene list.
#'
#' @return A named character vector with entries \code{calcium_channels},
#'   \code{calcium_binding}, and \code{calcium_transport}.
get_calcium_patterns <- function() {
  c(
    calcium_channels  = "^(CACNA|CACNB|CACNG|TRPC|TRPV|TRPM|ORAI|STIM)",
    calcium_binding   = "^(CALM|CAMK|CAN|S100)",
    calcium_transport = "^(ATP2A|ATP2B|ATP2C|SLC8A|MCU|MICU)"
  )
}

# ------------------------------------------------------------------------------
# 6. Merge gene list with differential expression results
# ------------------------------------------------------------------------------

#' Left-join a curated gene list with DEG results and annotate significance.
#'
#' @param gene_list A tibble containing at least a \code{gene_symbol} column,
#'   plus any additional annotation columns (e.g. \code{category}).
#' @param deg_df A tibble of differential expression results with columns
#'   \code{gene_symbol}, \code{log2fc}, \code{pvalue}, and \code{padj}.
#' @param lfc_thresh Absolute log2 fold-change threshold for significance
#'   (default 1).
#' @param padj_thresh Adjusted p-value threshold for significance (default
#'   0.05).
#'
#' @return The input \code{gene_list} left-joined with \code{deg_df}, plus
#'   logical column \code{sig} and character column \code{direction}
#'   ("Up", "Down", or "NS").
merge_with_deg <- function(gene_list,
                           deg_df,
                           lfc_thresh  = 1,
                           padj_thresh = 0.05) {

  stopifnot(
    "gene_symbol" %in% colnames(gene_list),
    all(c("gene_symbol", "log2fc", "pvalue", "padj") %in% colnames(deg_df))
  )

  gene_list %>%
    dplyr::left_join(deg_df, by = "gene_symbol") %>%
    dplyr::mutate(
      sig = (abs(log2fc) > lfc_thresh & padj < padj_thresh),
      direction = dplyr::case_when(
        sig & log2fc > 0 ~ "Up",
        sig & log2fc < 0 ~ "Down",
        TRUE             ~ "NS"
      )
    )
}

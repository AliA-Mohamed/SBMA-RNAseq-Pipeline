# plot_utils.R
# Plotting helper functions for the SBMA RNAseq pipeline.
# ---------------------------------------------------------------------------

library(ggplot2)
library(ggrepel)
library(cowplot)
library(RColorBrewer)
library(dplyr)
library(tidyr)

# ---- Default color palette --------------------------------------------------

#' Default color configuration for all pipeline plots
#'
#' @return A named list of hex color strings.
#' @export
default_plot_colors <- function() {
  list(
    up           = "#bb0c00",
    down         = "#00AFBB",
    ns           = "grey60",
    heatmap_low  = "#2166AC",
    heatmap_mid  = "white",
    heatmap_high = "#B2182B"
  )
}


# ---- Publication theme ------------------------------------------------------

#' Clean ggplot2 theme for publication-quality figures
#'
#' Based on \code{theme_bw} with grid lines removed, legend positioned at the
#' bottom, and clean axis lines.
#'
#' @param base_size Numeric. Base font size (default 14).
#' @return A \code{ggplot2} theme object.
#' @export
get_theme_publication <- function(base_size = 14) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      # Remove grid lines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),

      # Clean panel border
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),

      # Axis
      axis.line       = element_line(colour = "black", linewidth = 0.4),
      axis.ticks      = element_line(colour = "black", linewidth = 0.4),
      axis.text        = element_text(colour = "black"),
      axis.title       = element_text(face = "bold"),

      # Legend
      legend.position  = "bottom",
      legend.background = element_blank(),
      legend.key        = element_blank(),

      # Strip (for faceted plots)
      strip.background = element_rect(fill = "grey95", colour = "black",
                                      linewidth = 0.4),
      strip.text       = element_text(face = "bold"),

      # Plot title
      plot.title   = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
}


# ---- Volcano plot -----------------------------------------------------------

#' Volcano plot for differential expression results
#'
#' @param df        Data frame with columns: \code{gene_symbol}, \code{log2fc},
#'                  \code{padj}, \code{diffexpressed}.
#' @param title     Character. Plot title.
#' @param lfc_thresh  Numeric. Log2 fold-change threshold for dashed lines
#'                    (default 1).
#' @param padj_thresh Numeric. Adjusted p-value threshold (default 0.05).
#' @param colors    Named list with elements \code{up}, \code{down}, \code{ns}.
#' @param top_n_label Integer. Number of top genes (by padj) to label
#'                    (default 15).
#' @return A \code{ggplot} object.
#' @export
make_volcano <- function(df,
                         title,
                         lfc_thresh   = 1,
                         padj_thresh  = 0.05,
                         colors       = default_plot_colors(),
                         top_n_label  = 15) {

  # -- Input validation -------------------------------------------------------
  required_cols <- c("gene_symbol", "log2fc", "padj", "diffexpressed")
  missing_cols  <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0L) {
    stop(sprintf("make_volcano: missing columns in df: %s",
                 paste(missing_cols, collapse = ", ")),
         call. = FALSE)
  }

  # -- Derived columns --------------------------------------------------------
  df <- df %>%
    mutate(
      neglog10padj = -log10(pmax(.data$padj, 1e-300)),
      diffexpressed = factor(.data$diffexpressed,
                             levels = c("Upregulated", "Downregulated", "NS"))
    )

  # -- Count significant genes ------------------------------------------------
  n_up   <- sum(df$diffexpressed == "Upregulated", na.rm = TRUE)
  n_down <- sum(df$diffexpressed == "Downregulated", na.rm = TRUE)

  # -- Select top genes to label ----------------------------------------------
  sig_df <- df %>%
    filter(.data$diffexpressed != "NS") %>%
    arrange(.data$padj) %>%
    head(top_n_label)

  # -- Axis limits for annotation placement -----------------------------------
  x_range <- range(df$log2fc, na.rm = TRUE)
  y_max   <- max(df$neglog10padj, na.rm = TRUE)

  # -- Color mapping ----------------------------------------------------------
  color_map <- c(
    "Upregulated"   = colors$up,
    "Downregulated" = colors$down,
    "NS"            = colors$ns
  )

  # -- Build plot -------------------------------------------------------------
  p <- ggplot(df, aes(x = .data$log2fc,
                      y = .data$neglog10padj,
                      colour = .data$diffexpressed)) +
    geom_point(alpha = 0.6, size = 1.2) +

    # Threshold lines
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh),
               linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_hline(yintercept = -log10(padj_thresh),
               linetype = "dashed", colour = "grey40", linewidth = 0.5) +

    # Gene labels
    geom_text_repel(
      data          = sig_df,
      aes(label = .data$gene_symbol),
      size          = 3,
      max.overlaps  = 20,
      show.legend   = FALSE,
      fontface      = "italic",
      segment.color = "grey50",
      segment.size  = 0.3
    ) +

    # Significance count annotations
    annotate("text",
             x     = x_range[2] * 0.85,
             y     = y_max * 0.95,
             label = paste(n_up, "UP"),
             colour = colors$up,
             fontface = "bold", size = 4, hjust = 1) +
    annotate("text",
             x     = x_range[1] * 0.85,
             y     = y_max * 0.95,
             label = paste(n_down, "DOWN"),
             colour = colors$down,
             fontface = "bold", size = 4, hjust = 0) +

    # Scales
    scale_colour_manual(values = color_map, name = NULL) +

    # Labels
    labs(
      title = title,
      x     = expression(log[2]~"Fold Change"),
      y     = expression(-log[10]~"adjusted P")
    ) +
    get_theme_publication() +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))

  p
}


# ---- Heatmap ----------------------------------------------------------------

#' Heatmap of log2 fold-changes across conditions
#'
#' Uses \pkg{pheatmap} when available; falls back to a \code{ggplot2}
#' tile-based heatmap otherwise.
#'
#' @param gene_df      Data frame with \code{gene_symbol}, a category column,
#'                     and one or more LFC value columns.
#' @param value_cols   Character vector. Column names holding LFC values.
#' @param category_col Character. Column name for row annotation categories.
#' @param title        Character. Plot title.
#' @param colors       Named list with \code{heatmap_low}, \code{heatmap_mid},
#'                     \code{heatmap_high}.
#' @return A \code{pheatmap} object or a \code{ggplot} object.
#' @export
make_heatmap <- function(gene_df,
                         value_cols,
                         category_col,
                         title,
                         colors = default_plot_colors()) {

  # -- Build matrix and annotation --------------------------------------------
  mat <- as.matrix(gene_df[, value_cols, drop = FALSE])
  rownames(mat) <- gene_df$gene_symbol

  annotation_row <- data.frame(
    Category = gene_df[[category_col]],
    row.names = gene_df$gene_symbol
  )
  colnames(annotation_row) <- category_col

  # -- Color palette ----------------------------------------------------------
  hm_palette <- colorRampPalette(
    c(colors$heatmap_low, colors$heatmap_mid, colors$heatmap_high)
  )(100)

  # -- Attempt pheatmap -------------------------------------------------------
  use_pheatmap <- requireNamespace("pheatmap", quietly = TRUE)

  if (use_pheatmap) {
    p <- pheatmap::pheatmap(
      mat,
      main            = title,
      color           = hm_palette,
      cluster_cols    = FALSE,
      cluster_rows    = TRUE,
      annotation_row  = annotation_row,
      fontsize_row    = 8,
      cellwidth       = 40,
      border_color    = NA,
      silent          = TRUE
    )
    return(p)
  }

  # -- ggplot2 fallback -------------------------------------------------------
  message("make_heatmap: pheatmap not available, using ggplot2 fallback.")

  plot_df <- gene_df %>%
    select(all_of(c("gene_symbol", category_col, value_cols))) %>%
    pivot_longer(
      cols      = all_of(value_cols),
      names_to  = "condition",
      values_to = "lfc"
    ) %>%
    mutate(
      gene_symbol = factor(.data$gene_symbol,
                           levels = rev(unique(gene_df$gene_symbol)))
    )

  # Determine symmetric scale limits
  lfc_max <- max(abs(plot_df$lfc), na.rm = TRUE)

  p <- ggplot(plot_df, aes(x = .data$condition,
                           y = .data$gene_symbol,
                           fill = .data$lfc)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    facet_grid(
      as.formula(paste(category_col, "~ .")),
      scales = "free_y",
      space  = "free_y"
    ) +
    scale_fill_gradient2(
      low      = colors$heatmap_low,
      mid      = colors$heatmap_mid,
      high     = colors$heatmap_high,
      midpoint = 0,
      limits   = c(-lfc_max, lfc_max),
      name     = expression(log[2]~FC)
    ) +
    labs(title = title, x = NULL, y = NULL) +
    get_theme_publication() +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1),
      axis.text.y  = element_text(size = 8),
      strip.text.y = element_text(angle = 0)
    )

  p
}


# ---- Forest plot ------------------------------------------------------------

#' Forest-style dot plot of log2 fold-changes
#'
#' @param gene_df      Data frame with \code{gene_symbol}, a category column,
#'                     LFC columns, and a significance column.
#' @param value_cols   Character vector. Column names holding LFC values.
#' @param category_col Character. Column name for faceting categories.
#' @param sig_col      Character. Column name indicating significance status.
#' @param colors       Named list (uses \code{up}, \code{down}, \code{ns}).
#' @return A \code{ggplot} object.
#' @export
make_forest_plot <- function(gene_df,
                             value_cols,
                             category_col,
                             sig_col,
                             colors = default_plot_colors()) {

  # -- Pivot to long format ---------------------------------------------------
  plot_df <- gene_df %>%
    select(all_of(c("gene_symbol", category_col, sig_col, value_cols))) %>%
    pivot_longer(
      cols      = all_of(value_cols),
      names_to  = "condition",
      values_to = "lfc"
    ) %>%
    mutate(
      gene_symbol = factor(.data$gene_symbol,
                           levels = rev(unique(gene_df$gene_symbol)))
    )

  # -- Color mapping ----------------------------------------------------------
  sig_levels <- sort(unique(plot_df[[sig_col]]))
  sig_colors <- setNames(
    colorRampPalette(c(colors$down, colors$ns, colors$up))(length(sig_levels)),
    sig_levels
  )

  # -- Build plot -------------------------------------------------------------
  p <- ggplot(plot_df, aes(x     = .data$lfc,
                           y     = .data$gene_symbol,
                           colour = .data[[sig_col]])) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40",
               linewidth = 0.5) +
    geom_point(size = 2.5, alpha = 0.8) +
    facet_grid(
      as.formula(paste(category_col, "~ condition")),
      scales = "free_y",
      space  = "free_y"
    ) +
    scale_colour_manual(values = sig_colors, name = "Significance") +
    labs(
      x = expression(log[2]~"Fold Change"),
      y = NULL
    ) +
    get_theme_publication() +
    theme(
      axis.text.y  = element_text(size = 8),
      strip.text.y = element_text(angle = 0)
    )

  p
}


# ---- LFC barplot ------------------------------------------------------------

#' Grouped barplot of log2 fold-changes with error bars
#'
#' @param df       Data frame with columns for x, y (mean LFC), fill grouping,
#'                 and optionally \code{sem} and \code{sig_label}.
#' @param x_col    Character. Column for x-axis (e.g., gene or category).
#' @param y_col    Character. Column for bar height (mean LFC).
#' @param fill_col Character. Column for fill grouping (e.g., condition).
#' @param title    Character. Plot title.
#' @param colors   Named list (not used directly; fill scale is automatic).
#' @return A \code{ggplot} object.
#' @export
make_barplot_lfc <- function(df,
                             x_col,
                             y_col,
                             fill_col,
                             title,
                             colors = default_plot_colors()) {

  # -- Expect sem and sig_label columns if available --------------------------
  has_sem  <- "sem" %in% colnames(df)
  has_sig  <- "sig_label" %in% colnames(df)

  p <- ggplot(df, aes(x    = .data[[x_col]],
                      y    = .data[[y_col]],
                      fill = .data[[fill_col]])) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.9) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.5)

  # Error bars (SEM)
  if (has_sem) {
    p <- p +
      geom_errorbar(
        aes(ymin = .data[[y_col]] - .data$sem,
            ymax = .data[[y_col]] + .data$sem),
        position = position_dodge(width = 0.8),
        width    = 0.25,
        linewidth = 0.4
      )
  }

  # Significance stars
  if (has_sig) {
    p <- p +
      geom_text(
        aes(
          label = .data$sig_label,
          y     = .data[[y_col]] + ifelse(has_sem, .data$sem, 0) +
                  max(abs(df[[y_col]]), na.rm = TRUE) * 0.03
        ),
        position = position_dodge(width = 0.8),
        vjust    = 0,
        size     = 4,
        show.legend = FALSE
      )
  }

  p <- p +
    labs(
      title = title,
      x     = NULL,
      y     = expression(log[2]~"Fold Change"),
      fill  = NULL
    ) +
    get_theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p
}


# ---- OXPHOS complex summary barplot -----------------------------------------

#' Barplot summarizing LFC across OXPHOS complexes
#'
#' @param complex_stats Data frame in long format with columns: \code{complex},
#'                      \code{condition}, \code{mean_lfc}, \code{sem},
#'                      \code{pval}.
#' @param colors        Named list (uses \code{up}, \code{down}).
#' @return A \code{ggplot} object.
#' @export
make_complex_summary <- function(complex_stats,
                                 colors = default_plot_colors()) {

  # -- Input validation -------------------------------------------------------
  required_cols <- c("complex", "condition", "mean_lfc", "sem", "pval")
  missing_cols  <- setdiff(required_cols, colnames(complex_stats))
  if (length(missing_cols) > 0L) {
    stop(sprintf("make_complex_summary: missing columns: %s",
                 paste(missing_cols, collapse = ", ")),
         call. = FALSE)
  }

  # -- Significance stars -----------------------------------------------------
  complex_stats <- complex_stats %>%
    mutate(
      sig_star = case_when(
        .data$pval < 0.001 ~ "***",
        .data$pval < 0.01  ~ "**",
        .data$pval < 0.05  ~ "*",
        TRUE               ~ ""
      ),
      bar_color = ifelse(.data$mean_lfc >= 0, "pos", "neg")
    )

  p <- ggplot(complex_stats, aes(x    = .data$complex,
                                 y    = .data$mean_lfc,
                                 fill = .data$bar_color)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.9) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.5) +

    # Error bars
    geom_errorbar(
      aes(ymin = .data$mean_lfc - .data$sem,
          ymax = .data$mean_lfc + .data$sem),
      position  = position_dodge(width = 0.8),
      width     = 0.25,
      linewidth = 0.4
    ) +

    # Significance stars
    geom_text(
      aes(
        label = .data$sig_star,
        y     = .data$mean_lfc + sign(.data$mean_lfc) *
                (.data$sem + max(abs(complex_stats$mean_lfc), na.rm = TRUE) * 0.05)
      ),
      position = position_dodge(width = 0.8),
      vjust    = ifelse(complex_stats$mean_lfc >= 0, 0, 1),
      size     = 5,
      show.legend = FALSE
    ) +

    # Facet by condition if multiple conditions present
    {
      if (length(unique(complex_stats$condition)) > 1L) {
        facet_wrap(~ condition, scales = "free_y")
      }
    } +

    scale_fill_manual(
      values = c("pos" = colors$up, "neg" = colors$down),
      guide  = "none"
    ) +
    labs(
      title = "OXPHOS Complex Mean LFC",
      x     = NULL,
      y     = expression("Mean" ~ log[2] ~ "Fold Change")
    ) +
    get_theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p
}


# ---- Save plot utility ------------------------------------------------------

#' Save a plot to disk in one or more formats
#'
#' Handles both \code{ggplot} objects (via \code{ggsave}) and \code{pheatmap}
#' objects (via base-R device functions).
#'
#' @param plot    A \code{ggplot} or \code{pheatmap} object.
#' @param path    Character. Output file path \emph{without} extension.
#' @param width   Numeric. Figure width in inches.
#' @param height  Numeric. Figure height in inches.
#' @param dpi     Numeric. Resolution for raster formats (default 300).
#' @param formats Character vector. File formats to produce
#'                (default \code{c("pdf", "png")}).
#' @return Invisible \code{NULL}. Called for side effects.
#' @export
save_plot <- function(plot,
                      path,
                      width,
                      height,
                      dpi     = 300,
                      formats = c("pdf", "png")) {

  # Ensure output directory exists
  out_dir <- dirname(path)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  is_pheatmap <- inherits(plot, "pheatmap")

  for (fmt in formats) {
    filepath <- paste0(path, ".", fmt)

    if (is_pheatmap) {
      # pheatmap objects must be drawn inside an open device
      if (fmt == "pdf") {
        pdf(filepath, width = width, height = height)
      } else if (fmt == "png") {
        png(filepath, width = width, height = height, units = "in", res = dpi)
      } else {
        warning(sprintf("save_plot: unsupported format '%s' for pheatmap; skipping.",
                        fmt))
        next
      }
      grid::grid.newpage()
      grid::grid.draw(plot$gtable)
      dev.off()
    } else {
      # ggplot or other grid-compatible objects
      ggsave(
        filename = filepath,
        plot     = plot,
        width    = width,
        height   = height,
        dpi      = dpi,
        bg       = "white"
      )
    }

    message(sprintf("Saved: %s", filepath))
  }

  invisible(NULL)
}

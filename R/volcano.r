#' Volcano plot for GWAS results
#'
#' @param x A GWASFormatter object or a data.frame/tibble containing GWAS results.
#' @param phenotype_label Optional label for the plot title.
#' @param ... Additional arguments passed to methods.
#' @return A ggplot2 object with the volcano plot.
#' @export
volcano <- function(x, phenotype_label = NULL, ...) {
  UseMethod("volcano")
}

#' @describeIn volcano Method for GWASFormatter objects
#' @export
volcano.GWASFormatter <- function(x, phenotype_label = NULL, ...) {
  df <- x$data %>%
    dplyr::select(
      BETA,
      PVALUE
    ) %>%
    dplyr::collect()

  volcano.data.frame(df, phenotype_label, ...)
}

#' @describeIn volcano Method for data.frame/tibble objects
#' @export
volcano.data.frame <- function(x, phenotype_label = NULL, ...) {
  # Determine x-axis variable - prefer posterior mean if available, otherwise use BETA
  x_var <- if ("pm" %in% names(x)) "pm" else "BETA"
  x_lab <- if (x_var == "pm") {
    "Posterior mean log odds-ratio"
  } else {
    "Effect size (BETA)"
  }
  # Determine y-axis variable - prefer lfsr if available, otherwise use PVALUE
  if ("lfsr" %in% names(x)) {
    y_var <- "lfsr"
    y_lab <- "-log10(local false-sign rate)"
    threshold_y <- -log10(0.05)
  } else if ("PVALUE" %in% names(x)) {
    y_var <- "PVALUE"
    y_lab <- "-log10(p-value)"
    threshold_y <- -log10(5e-8)
  } else {
    stop("Data must contain either 'lfsr' or 'PVALUE' column")
  }
  # Create the plot
  p <- ggplot2::ggplot(
    x,
    ggplot2::aes_string(x = x_var, y = paste0("-log10(", y_var, ")"))
  ) +
    ggplot2::geom_hline(
      yintercept = threshold_y,
      linetype = "dashed",
      alpha = 0.7
    ) +
    ggrastr::rasterize(
      ggplot2::geom_point(
        ggplot2::aes(color = if ("group" %in% names(x)) group else NULL),
        alpha = 0.75
      ),
      dpi = 300,
      dev = "ragg_png"
    ) +
    ggplot2::labs(
      x = x_lab,
      y = y_lab
    ) +
    cowplot::theme_cowplot()
  # Add color scale if group variable exists
  if ("group" %in% names(x)) {
    p <- p + ggplot2::scale_color_brewer(palette = "Set2")
  } else {
    p <- p + ggplot2::scale_color_manual(values = "black", guide = "none")
  }
  # Add text labels if label column exists
  if ("label" %in% names(x) && requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p +
      ggrepel::geom_text_repel(
        data = x[!is.na(x$label) & x$label != "", ],
        ggplot2::aes_string(label = "label"),
        size = 2.5,
        max.overlaps = 6,
        label.padding = 0.25,
        max.time = 4.5,
        nudge_x = -0.2,
        nudge_y = 1.1
      )
  }

  # Add title if provided
  if (!is.null(phenotype_label)) {
    p <- p + ggplot2::labs(title = phenotype_label)
  }

  return(p)
}

#' @describeIn volcano Alias for data.frame method
#' @export
volcano.tbl_df <- function(x, phenotype_label = NULL, ...) {
  volcano.data.frame(x, phenotype_label, ...)
}

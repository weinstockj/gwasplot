# Global variables to avoid NSE warnings in R CMD check
utils::globalVariables(c(
  "CHROM",
  "POS",
  "PVALUE",
  "CHROM_index",
  "chr_len",
  "tot",
  "POScum",
  "center",
  "aes",
  "gene_name",
  "locus_id",
  "is_lead",
  "label_text",
  "label_is_gene",
  "label_expr",
  "label_is_highlight",
  "label_y",
  "label_nudge_y",
  "LOGP",
  "LOGP_plot",
  ".label_cluster",
  ".cluster_gap",
  ".label_slot"
))

# Helper function to create chromosome lookup table
create_chrom_lookup = function() {
  tibble::tibble(
    CHROM = c(glue::glue("chr{1:22}"), "chrX"),
    CHROM_index = 1:23
  )
}

# Helper function to prepare GWAS data for plotting
prepare_manhattan_data = function(
  gwas_data,
  chrom_lookup,
  lower_logp_threshold = 3.0
) {
  scaling = 1e8

  extra_cols = intersect("gene_name", colnames(gwas_data))
  base_cols = c("CHROM", "POS", "PVALUE")

  base_data = gwas_data %>%
    dplyr::select(dplyr::all_of(c(base_cols, extra_cols))) %>%
    dplyr::inner_join(chrom_lookup, by = "CHROM")

  chrom_offsets = base_data %>%
    dplyr::group_by(CHROM_index, CHROM) %>%
    dplyr::summarise(
      chr_len = (max(POS, na.rm = TRUE) - min(POS, na.rm = TRUE)) / scaling,
      .groups = "drop"
    )

  if (inherits(chrom_offsets, "tbl_sql")) {
    chrom_offsets = chrom_offsets %>%
      dbplyr::window_order(CHROM_index) %>%
      dplyr::mutate(tot = cumsum(chr_len) - chr_len)
  } else {
    chrom_offsets = chrom_offsets %>%
      dplyr::arrange(CHROM_index) %>%
      dplyr::mutate(tot = cumsum(chr_len) - chr_len)
  }

  prepared = base_data %>%
    dplyr::inner_join(chrom_offsets, by = c("CHROM_index", "CHROM")) %>%
    dplyr::arrange(CHROM_index, POS) %>%
    dplyr::mutate(POScum = POS / scaling + tot) %>%
    dplyr::filter(-log10(PVALUE) > lower_logp_threshold)

  return(prepared)
}

# Helper function to calculate axis positions
calculate_axis_positions = function(prepared_data) {
  axis = prepared_data %>%
    dplyr::group_by(CHROM_index, CHROM) %>%
    dplyr::summarize(
      center = (max(POScum) + min(POScum)) / 2,
      .groups = "drop"
    ) %>% # integer overflow concern
    dplyr::arrange(CHROM_index, center)

  return(axis)
}

# Compress the upper y-axis while preserving original -log10(p) tick labels.
transform_manhattan_y = function(logp, y_axis_break = NULL, y_axis_break_scale = 0.2) {
  if (is.null(y_axis_break)) {
    return(logp)
  }

  if (
    !is.numeric(y_axis_break) ||
      length(y_axis_break) != 1 ||
      !is.finite(y_axis_break) ||
      y_axis_break <= 0
  ) {
    cli::cli_abort("{.arg y_axis_break} must be a single positive finite number or NULL.")
  }
  if (
    !is.numeric(y_axis_break_scale) ||
      length(y_axis_break_scale) != 1 ||
      !is.finite(y_axis_break_scale) ||
      y_axis_break_scale <= 0 ||
      y_axis_break_scale > 1
  ) {
    cli::cli_abort("{.arg y_axis_break_scale} must be a single number in (0, 1].")
  }

  y_axis_break + pmin(logp - y_axis_break, 0) +
    pmax(logp - y_axis_break, 0) * y_axis_break_scale
}

add_manhattan_y_coordinates = function(
  prepared_data,
  y_axis_break = NULL,
  y_axis_break_scale = 0.2
) {
  prepared_data %>%
    dplyr::mutate(
      LOGP = -log10(PVALUE),
      LOGP_plot = transform_manhattan_y(
        LOGP,
        y_axis_break = y_axis_break,
        y_axis_break_scale = y_axis_break_scale
      )
    )
}

calculate_y_axis = function(
  prepared_data,
  y_max_cap = 300,
  y_axis_break = NULL,
  y_axis_break_scale = 0.2
) {
  logp <- pmin(prepared_data$LOGP, y_max_cap)
  logp <- logp[is.finite(logp)]
  if (length(logp) == 0) {
    return(list(breaks = 0, labels = 0))
  }

  max_logp <- max(logp, na.rm = TRUE)
  min_logp <- min(logp, na.rm = TRUE)
  if (is.null(y_axis_break) || y_axis_break >= max_logp) {
    tick_start <- ceiling(min_logp / 4) * 4
    breaks <- if (tick_start <= max_logp) {
      seq(tick_start, max_logp, by = 4)
    } else {
      numeric(0)
    }
    original_breaks <- unique(c(min_logp, breaks))
    return(list(
      breaks = original_breaks,
      labels = original_breaks,
      limits = c(min_logp, max_logp)
    ))
  }

  lower_tick_start <- ceiling(min_logp / 4) * 4
  lower_tick_end <- floor(y_axis_break / 4) * 4
  lower_breaks <- if (lower_tick_start <= lower_tick_end) {
    seq(lower_tick_start, lower_tick_end, by = 4)
  } else {
    numeric(0)
  }
  upper_breaks <- pretty(c(y_axis_break, max_logp), n = 4)
  upper_breaks <- upper_breaks[upper_breaks > y_axis_break & upper_breaks <= max_logp]
  if (!any(abs(upper_breaks - max_logp) < .Machine$double.eps^0.5)) {
    upper_breaks <- c(upper_breaks, max_logp)
  }

  original_breaks <- unique(c(min_logp, lower_breaks, upper_breaks))
  list(
    breaks = transform_manhattan_y(
      original_breaks,
      y_axis_break = y_axis_break,
      y_axis_break_scale = y_axis_break_scale
    ),
    labels = original_breaks,
    limits = transform_manhattan_y(
      c(min_logp, max_logp),
      y_axis_break = y_axis_break,
      y_axis_break_scale = y_axis_break_scale
    )
  )
}

add_y_axis_break_slash = function(
  plot,
  prepared_data,
  y_axis_break,
  y_limits = NULL
) {
  if (is.null(y_axis_break) || y_axis_break >= max(prepared_data$LOGP, na.rm = TRUE)) {
    return(plot)
  }

  x_range <- range(prepared_data$POScum, na.rm = TRUE)
  x_span <- diff(x_range)
  if (!is.finite(x_span) || x_span <= 0) {
    x_span <- 1
  }

  y_range <- range(prepared_data$LOGP_plot, na.rm = TRUE)
  y_span <- diff(y_range)
  if (!is.finite(y_span) || y_span <= 0) {
    y_span <- 1
  }

  x_axis <- x_range[[1]]
  x_half_width <- x_span * 0.0075
  slash_y <- transform_manhattan_y(
    y_axis_break,
    y_axis_break = y_axis_break,
    y_axis_break_scale = 1
  )
  slash_gap <- y_span * 0.012
  slash_half_height <- y_span * 0.018

  axis_data <- tibble::tibble(
    x = x_axis,
    xend = x_axis,
    y = min(prepared_data$LOGP_plot, na.rm = TRUE),
    yend = max(prepared_data$LOGP_plot, na.rm = TRUE)
  )
  slash_data <- tibble::tibble(
    x = c(x_axis - x_half_width, x_axis - x_half_width),
    xend = c(x_axis + x_half_width, x_axis + x_half_width),
    y = c(slash_y - slash_half_height, slash_y + slash_gap - slash_half_height),
    yend = c(slash_y + slash_half_height, slash_y + slash_gap + slash_half_height)
  )

  plot +
    ggplot2::theme(axis.line.y = ggplot2::element_blank()) +
    ggplot2::coord_cartesian(ylim = y_limits, clip = "off") +
    ggplot2::geom_segment(
      data = axis_data,
      mapping = ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE,
      linewidth = 0.6,
      lineend = "square",
      color = "black"
    ) +
    ggplot2::geom_segment(
      data = slash_data,
      mapping = ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE,
      linewidth = 0.45,
      lineend = "round",
      color = "black"
    )
}

# Build cytoband-like labels (e.g., "chr4p16.3") for fallback labeling.
derive_cytoband_labels = function(df) {
  if (!exists("ideogram", inherits = TRUE)) {
    return(rep(NA_character_, nrow(df)))
  }

  bands <- ideogram %>%
    dplyr::select(chrom, start, end, name)

  joined <- df %>%
    dplyr::mutate(.row_id = dplyr::row_number()) %>%
    dplyr::left_join(
      bands,
      by = dplyr::join_by(CHROM == chrom, dplyr::between(POS, start, end))
    ) %>%
    dplyr::group_by(.row_id) %>%
    dplyr::summarise(
      band_name = dplyr::first(stats::na.omit(name)),
      .groups = "drop"
    ) %>%
    dplyr::right_join(
      df %>% dplyr::mutate(.row_id = dplyr::row_number()),
      by = ".row_id"
    ) %>%
    dplyr::arrange(.row_id)

  ifelse(
    is.na(joined$band_name),
    joined$CHROM,
    paste0(joined$CHROM, joined$band_name)
  )
}

# Escape text for plotmath and optionally return an italicized expression string.
as_plotmath_label = function(text, italic = FALSE) {
  escaped <- gsub("\\\\", "\\\\\\\\", text)
  escaped <- gsub("'", "\\\\'", escaped)
  if (italic) {
    paste0("italic('", escaped, "')")
  } else {
    paste0("'", escaped, "'")
  }
}

# Select points to label and build label text/expression metadata.
select_label_points = function(
  prepared_data,
  label_top_n,
  label_strategy = c("top_n", "lead_per_locus"),
  label_locus_window_kb = 500,
  highlight_genes = NULL,
  italic_gene_labels = TRUE,
  label_pvalue_threshold = 5e-8,
  force_highlight_labels = TRUE
) {
  if (is.null(label_top_n) || label_top_n <= 0) {
    return(NULL)
  }

  label_strategy <- match.arg(label_strategy)

  # Only label variants that clear the significance threshold. Highlighted genes
  # bypass the gate by default so explicitly requested genes are displayed.
  has_gene_col <- "gene_name" %in% names(prepared_data)
  is_highlight_row <- if (has_gene_col && !is.null(highlight_genes)) {
    !is.na(prepared_data$gene_name) &
      prepared_data$gene_name %in% highlight_genes
  } else {
    rep(FALSE, nrow(prepared_data))
  }

  candidates <- prepared_data %>%
    dplyr::filter(
      PVALUE < label_pvalue_threshold |
        (force_highlight_labels & is_highlight_row)
    ) %>%
    dplyr::arrange(PVALUE)

  if (nrow(candidates) == 0) {
    return(NULL)
  }

  if (label_strategy == "lead_per_locus" && nrow(candidates) > 0) {
    candidates <- identify_loci.data.frame(
      candidates,
      window_kb = label_locus_window_kb,
      pvalue_threshold = Inf
    ) %>%
      dplyr::filter(is_lead)
  }

  # Top N by p-value, while keeping explicitly highlighted genes even if they
  # would otherwise rank below the cutoff or be removed by lead-locus thinning.
  top_df <- candidates %>%
    dplyr::slice_min(PVALUE, n = label_top_n, with_ties = FALSE)
  highlight_df <- if (
    force_highlight_labels && has_gene_col && !is.null(highlight_genes)
  ) {
    prepared_data %>%
      dplyr::filter(!is.na(gene_name) & gene_name %in% highlight_genes)
  } else {
    candidates[FALSE, , drop = FALSE]
  }

  label_df <- dplyr::distinct(dplyr::bind_rows(top_df, highlight_df))

  if (nrow(label_df) == 0) {
    return(NULL)
  }

  has_gene <- "gene_name" %in% names(label_df)
  gene_vals <- if (has_gene) {
    label_df$gene_name
  } else {
    rep(NA_character_, nrow(label_df))
  }
  band_vals <- derive_cytoband_labels(label_df)

  label_df <- label_df %>%
    dplyr::mutate(
      label_is_gene = !is.na(gene_vals) & gene_vals != "",
      label_text = ifelse(label_is_gene, gene_vals, band_vals),
      label_text = ifelse(
        is.na(label_text) | label_text == "",
        CHROM,
        label_text
      ),
      label_is_highlight = label_is_gene &
        !is.null(highlight_genes) &
        gene_vals %in% highlight_genes,
      label_expr = ifelse(
        label_is_gene & italic_gene_labels,
        vapply(
          label_text,
          as_plotmath_label,
          FUN.VALUE = character(1),
          italic = TRUE
        ),
        vapply(
          label_text,
          as_plotmath_label,
          FUN.VALUE = character(1),
          italic = FALSE
        )
      )
    )

  label_df %>%
    dplyr::arrange(PVALUE) %>%
    dplyr::group_by(label_text) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
}

# Compute modest, data-aware starting nudges for Manhattan labels.
#
# `ggrepel` still does the final collision avoidance, but these nudges give it a
# better initial layout: isolated labels stay close to their SNPs, while labels
# in the same horizontal neighborhood get progressively stacked.
calculate_label_layout = function(
  label_df,
  prepared_data,
  y_max_cap = 300,
  label_nudge_y = NULL,
  label_max_y_nudge = NULL
) {
  if (is.null(label_df) || nrow(label_df) == 0) {
    return(label_df)
  }

  label_y <- if ("LOGP_plot" %in% names(label_df)) {
    label_df$LOGP_plot
  } else {
    -log10(label_df$PVALUE)
  }
  plotted_y <- if ("LOGP_plot" %in% names(prepared_data)) {
    prepared_data$LOGP_plot
  } else {
    pmin(-log10(prepared_data$PVALUE), y_max_cap)
  }
  plotted_y <- plotted_y[is.finite(plotted_y)]

  if (length(plotted_y) == 0) {
    y_span <- 1
  } else {
    y_span <- diff(range(plotted_y, na.rm = TRUE))
    if (!is.finite(y_span) || y_span <= 0) {
      y_span <- max(plotted_y, na.rm = TRUE)
    }
    if (!is.finite(y_span) || y_span <= 0) {
      y_span <- 1
    }
  }

  fixed_nudge <- if (is.null(label_nudge_y)) {
    NULL
  } else {
    max(0, label_nudge_y)
  }
  base_nudge <- max(0.25, min(0.8, y_span * 0.025))
  step_nudge <- max(0.35, min(0.9, y_span * 0.03))
  max_nudge <- if (is.null(label_max_y_nudge)) {
    max(1.25, min(4, y_span * 0.12))
  } else {
    max(0, label_max_y_nudge)
  }

  x_vals <- prepared_data$POScum[is.finite(prepared_data$POScum)]
  x_span <- if (length(x_vals) > 1) diff(range(x_vals, na.rm = TRUE)) else 1
  if (!is.finite(x_span) || x_span <= 0) {
    x_span <- 1
  }

  cluster_width <- max(0.05, x_span * 0.0125)

  label_df %>%
    dplyr::mutate(label_y = label_y) %>%
    dplyr::arrange(POScum, dplyr::desc(label_y)) %>%
    dplyr::mutate(
      .cluster_gap = c(Inf, diff(POScum)),
      .label_cluster = cumsum(.cluster_gap > cluster_width)
    ) %>%
    dplyr::group_by(.label_cluster) %>%
    dplyr::arrange(dplyr::desc(label_y), .by_group = TRUE) %>%
    dplyr::mutate(
      .label_slot = dplyr::row_number() - 1L,
      label_nudge_y = if (is.null(fixed_nudge)) {
        pmin(max_nudge, base_nudge + .label_slot * step_nudge)
      } else {
        fixed_nudge
      }
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(PVALUE)
}

# Helper function to create the Manhattan plot
create_manhattan_plot = function(
  prepared_data,
  axis_data,
  pvalue_threshold = 5e-8,
  y_max_cap = 300,
  y_axis_break = NULL,
  y_axis_break_scale = 0.2,
  label_top_n = NULL,
  label_strategy = c("top_n", "lead_per_locus"),
  label_locus_window_kb = 500,
  label_pvalue_threshold = NULL,
  label_nudge_y = NULL,
  force_highlight_labels = TRUE,
  label_repel_direction = c("y", "both", "x"),
  italic_gene_labels = TRUE,
  highlight_genes = NULL,
  highlight_color = "red3",
  label_color = "black",
  label_size = 3,
  highlight_label_size = NULL,
  label_segment_alpha = 0.65,
  label_max_y_nudge = NULL,
  base_size = 7
) {
  label_strategy <- match.arg(label_strategy)
  label_repel_direction <- match.arg(label_repel_direction)

  # Default the labeling cutoff to the genome-wide significance line.
  if (is.null(label_pvalue_threshold)) {
    label_pvalue_threshold <- pvalue_threshold
  }

  prepared_data <- add_manhattan_y_coordinates(
    prepared_data,
    y_axis_break = y_axis_break,
    y_axis_break_scale = y_axis_break_scale
  )
  y_axis <- calculate_y_axis(
    prepared_data,
    y_max_cap = y_max_cap,
    y_axis_break = y_axis_break,
    y_axis_break_scale = y_axis_break_scale
  )

  p = ggplot2::ggplot(prepared_data, aes(x = POScum, y = LOGP_plot)) +
    # Show all points
    ggrastr::rasterise(
      ggplot2::geom_point(
        aes(color = as.factor(CHROM_index)),
        alpha = 0.8,
        size = .3
      ),
      dev = "ragg_png",
      dpi = 300
    ) +
    ggplot2::scale_color_manual(
      values = c(rep(c("#7173C9", "#01035F"), 11), "#7173C9")
    ) +
    # custom X axis:
    ggplot2::scale_x_continuous(
      label = stringr::str_replace(axis_data$CHROM, "chr", ""),
      breaks = axis_data$center,
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = y_axis$breaks,
      labels = scales::number(y_axis$labels, accuracy = 1),
      expand = ggplot2::expansion(mult = c(0, 0), add = c(0, 1))
    ) + # remove space between plot area and x axis
    # Custom the theme:
    ggplot2::geom_hline(
      yintercept = transform_manhattan_y(
        -log10(pvalue_threshold),
        y_axis_break = y_axis_break,
        y_axis_break_scale = y_axis_break_scale
      ),
      linetype = "dashed",
      color = "grey"
    ) +
    ggplot2::labs(y = expression(-log[10](pvalue))) +
    ggplot2::theme_bw(base_size = base_size, base_family = "Helvetica") +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = base_size * 4 / 7),
      axis.text.y = element_text(size = base_size * 5 / 7),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.border = element_blank(),
      axis.line.y = element_line(size = .6)
    )

  if (!is.null(y_axis_break) && y_axis_break < max(prepared_data$LOGP, na.rm = TRUE)) {
    p <- add_y_axis_break_slash(
      p,
      prepared_data,
      y_axis_break,
      y_limits = y_axis$limits
    )
  } else {
    p <- p + ggplot2::coord_cartesian(ylim = y_axis$limits)
  }

  # Add labels for top points using gene_name when available; otherwise cytoband.
  if (!is.null(label_top_n)) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      label_df <- select_label_points(
        prepared_data = prepared_data,
        label_top_n = label_top_n,
        label_strategy = label_strategy,
        label_locus_window_kb = label_locus_window_kb,
        highlight_genes = highlight_genes,
        italic_gene_labels = italic_gene_labels,
        label_pvalue_threshold = label_pvalue_threshold,
        force_highlight_labels = force_highlight_labels
      )

      if (!is.null(label_df) && nrow(label_df) > 0) {
        label_df <- calculate_label_layout(
          label_df = label_df,
          prepared_data = prepared_data,
          y_max_cap = y_max_cap,
          label_nudge_y = label_nudge_y,
          label_max_y_nudge = label_max_y_nudge
        )

        common_repel_args <- list(
          mapping = ggplot2::aes(
            x = POScum,
            y = LOGP_plot,
            label = label_expr
          ),
          parse = TRUE,
          max.overlaps = Inf,
          seed = 42,
          direction = label_repel_direction,
          box.padding = 0.25,
          point.padding = 0.15,
          force = 0.25,
          force_pull = 1.5,
          max.time = 1,
          max.iter = 20000,
          segment.size = 0.25,
          segment.alpha = label_segment_alpha,
          segment.color = "grey35",
          min.segment.length = 0
        )

        normal_df <- label_df %>% dplyr::filter(!label_is_highlight)
        highlight_df <- label_df %>% dplyr::filter(label_is_highlight)

        if (nrow(normal_df) > 0) {
          p <- p + do.call(
            ggrepel::geom_text_repel,
            c(
              list(
                data = normal_df,
                colour = label_color,
                size = label_size,
                nudge_y = normal_df$label_nudge_y,
                inherit.aes = FALSE
              ),
              common_repel_args
            )
          )
        }
        if (nrow(highlight_df) > 0) {
          p <- p + do.call(
            ggrepel::geom_text_repel,
            c(
              list(
                data = highlight_df,
                colour = highlight_color,
                fontface = "bold",
                size = if (is.null(highlight_label_size)) label_size else highlight_label_size,
                nudge_y = highlight_df$label_nudge_y,
                inherit.aes = FALSE
              ),
              common_repel_args
            )
          )
        }
      }
    } else {
      cli::cli_warn(
        "Install the {.pkg ggrepel} package to enable Manhattan labels."
      )
    }
  }

  return(p)
}

# Helper function to save the Manhattan plot
save_manhattan_plot = function(plot, output, ...) {
  log_info("Now rendering")

  ggplot2::ggsave(
    output,
    plot,
    ...,
    bg = "white"
  )

  log_info("done plotting.")
}

#' Plot a Manhattan plot from a gwas object
#'
#' @param gwas A gwas object containing the data to plot.
#' @param output The output file name.
#' @param lower_logp_threshold The lower threshold for the -log10(p-value) to plot. Default is 3.0.
#' @param y_axis_break Optional `-log10(p)` value where the upper y-axis should
#'   be compressed. Values below the break are unchanged; values above are
#'   compressed by `y_axis_break_scale`, with a slash marker drawn at the break.
#'   Default NULL (no compression).
#' @param y_axis_break_scale Compression factor for values above
#'   `y_axis_break`. Smaller values create more visual room below the break.
#'   Must be in (0, 1]. Default 0.2.
#' @param label_top_n Integer. If set and a `gene_name` column is present, annotate the
#'   top N variants (by p-value) with labels using ggrepel. If `gene_name` is missing,
#'   labels fall back to cytoband-like text (e.g., `chr4p16.3`). Default is NULL (no labels).
#' @param label_strategy Label selection mode: `"top_n"` (global top N by p-value)
#'   or `"lead_per_locus"` (one lead per locus, then top N). Default `"lead_per_locus"`.
#' @param label_locus_window_kb Locus window size in kilobases used when
#'   `label_strategy = "lead_per_locus"`. Default 500.
#' @param label_pvalue_threshold P-value cutoff for labeling: only variants with
#'   `PVALUE < label_pvalue_threshold` are eligible for labels. Genes in
#'   `highlight_genes` are labeled regardless when
#'   `force_highlight_labels = TRUE`. Default NULL, which uses the genome-wide
#'   significance line (5e-8).
#' @param label_nudge_y Vertical distance (in `-log10(p)` units) to lift gene
#'   labels above their points. Default NULL uses adaptive per-label nudges.
#'   Set to a number to force the same upward nudge for every label, or 0 to
#'   disable the upward nudge.
#' @param force_highlight_labels Logical. If TRUE, genes in `highlight_genes`
#'   are labeled even when they do not pass `label_pvalue_threshold` or
#'   `label_top_n`. Default TRUE.
#' @param label_repel_direction Direction passed to `ggrepel`: `"y"` keeps
#'   labels aligned above variants, `"both"` lets labels move horizontally and
#'   vertically, and `"x"` repels horizontally. Default `"y"`.
#' @param italic_gene_labels Logical. If TRUE, gene labels are italicized; cytoband
#'   fallback labels remain plain text. Default TRUE.
#' @param highlight_genes Optional character vector of gene symbols to highlight
#'   (e.g., novel genes). Matching labels are colored with `highlight_color`.
#' @param highlight_color Color for highlighted labels. Default `"red3"`.
#' @param label_color Color for non-highlighted labels. Default `"black"`.
#' @param label_size Text size for non-highlighted labels. Default 3.
#' @param highlight_label_size Text size for highlighted labels. Default NULL,
#'   which uses `label_size`.
#' @param label_segment_alpha Alpha transparency for label leader lines.
#'   Default 0.65.
#' @param label_max_y_nudge Maximum adaptive starting nudge for labels in
#'   `-log10(p)` units. `NULL` chooses a conservative value from the plotted
#'   y-range. Lower values keep labels closer to variants. Default NULL.
#' @param ... Additional arguments passed to `ggsave`.
#' @return NULL
#' @export
manhattan = function(
  gwas,
  output,
  lower_logp_threshold = 3.0,
  y_axis_break = NULL,
  y_axis_break_scale = 0.2,
  label_top_n = NULL,
  label_strategy = "lead_per_locus",
  label_locus_window_kb = 500,
  label_pvalue_threshold = NULL,
  label_nudge_y = NULL,
  force_highlight_labels = TRUE,
  label_repel_direction = "y",
  italic_gene_labels = TRUE,
  highlight_genes = NULL,
  highlight_color = "red3",
  label_color = "black",
  label_size = 3,
  highlight_label_size = NULL,
  label_segment_alpha = 0.65,
  label_max_y_nudge = NULL,
  ...
) {
  UseMethod("manhattan")
}

#' @export
manhattan.tbl_df = function(
  gwas,
  output,
  lower_logp_threshold = 3.0,
  y_axis_break = NULL,
  y_axis_break_scale = 0.2,
  label_top_n = NULL,
  label_strategy = "lead_per_locus",
  label_locus_window_kb = 500,
  label_pvalue_threshold = NULL,
  label_nudge_y = NULL,
  force_highlight_labels = TRUE,
  label_repel_direction = "y",
  italic_gene_labels = TRUE,
  highlight_genes = NULL,
  highlight_color = "red3",
  label_color = "black",
  label_size = 3,
  highlight_label_size = NULL,
  label_segment_alpha = 0.65,
  label_max_y_nudge = NULL,
  ...
) {
  manhattan.data.frame(
    gwas,
    output,
    lower_logp_threshold,
    y_axis_break = y_axis_break,
    y_axis_break_scale = y_axis_break_scale,
    label_top_n = label_top_n,
    label_strategy = label_strategy,
    label_locus_window_kb = label_locus_window_kb,
    label_pvalue_threshold = label_pvalue_threshold,
    label_nudge_y = label_nudge_y,
    force_highlight_labels = force_highlight_labels,
    label_repel_direction = label_repel_direction,
    italic_gene_labels = italic_gene_labels,
    highlight_genes = highlight_genes,
    highlight_color = highlight_color,
    label_color = label_color,
    label_size = label_size,
    highlight_label_size = highlight_label_size,
    label_segment_alpha = label_segment_alpha,
    label_max_y_nudge = label_max_y_nudge,
    ...
  )
}

#' @export
manhattan.data.frame = function(
  gwas,
  output,
  lower_logp_threshold = 3.0,
  y_axis_break = NULL,
  y_axis_break_scale = 0.2,
  label_top_n = NULL,
  label_strategy = "lead_per_locus",
  label_locus_window_kb = 500,
  label_pvalue_threshold = NULL,
  label_nudge_y = NULL,
  force_highlight_labels = TRUE,
  label_repel_direction = "y",
  italic_gene_labels = TRUE,
  highlight_genes = NULL,
  highlight_color = "red3",
  label_color = "black",
  label_size = 3,
  highlight_label_size = NULL,
  label_segment_alpha = 0.65,
  label_max_y_nudge = NULL,
  base_size = 7,
  ...
) {
  log_info("Now preparing to plot")

  # Create chromosome lookup
  chrom_lookup = create_chrom_lookup()

  # Prepare data
  prepared = prepare_manhattan_data(gwas, chrom_lookup, lower_logp_threshold)

  # Calculate axis positions
  axis = calculate_axis_positions(prepared)

  log_info(glue("Done preparing to plot {nrow(prepared)} SNPs."))

  # Create plot
  p = create_manhattan_plot(
    prepared,
    axis,
    y_axis_break = y_axis_break,
    y_axis_break_scale = y_axis_break_scale,
    label_top_n = label_top_n,
    label_strategy = label_strategy,
    label_locus_window_kb = label_locus_window_kb,
    label_pvalue_threshold = label_pvalue_threshold,
    label_nudge_y = label_nudge_y,
    force_highlight_labels = force_highlight_labels,
    label_repel_direction = label_repel_direction,
    italic_gene_labels = italic_gene_labels,
    highlight_genes = highlight_genes,
    highlight_color = highlight_color,
    label_color = label_color,
    label_size = label_size,
    highlight_label_size = highlight_label_size,
    label_segment_alpha = label_segment_alpha,
    label_max_y_nudge = label_max_y_nudge,
    base_size = base_size
  )

  # Save plot
  save_manhattan_plot(p, output, ...)
}

#' @export
manhattan.GWASFormatter = function(
  gwas,
  output,
  lower_logp_threshold = 3.0,
  y_axis_break = NULL,
  y_axis_break_scale = 0.2,
  label_top_n = NULL,
  label_strategy = "lead_per_locus",
  label_locus_window_kb = 500,
  label_pvalue_threshold = NULL,
  label_nudge_y = NULL,
  force_highlight_labels = TRUE,
  label_repel_direction = "y",
  italic_gene_labels = TRUE,
  highlight_genes = NULL,
  highlight_color = "red3",
  label_color = "black",
  label_size = 3,
  highlight_label_size = NULL,
  label_segment_alpha = 0.65,
  label_max_y_nudge = NULL,
  base_size = 7,
  ...
) {
  log_info("Now preparing to plot")

  # Create chromosome lookup and copy to database connection
  chrom_lookup = copy_to_if_missing(
    gwas$con,
    create_chrom_lookup(),
    "chrom_lookup"
  )

  # Prepare data
  prepared = prepare_manhattan_data(
    gwas$data,
    chrom_lookup,
    lower_logp_threshold
  )

  # Collect prepared data from database
  prepared = prepared %>%
    dplyr::collect(.)

  # Calculate axis positions after collecting the prepared points once
  axis = calculate_axis_positions(prepared)

  log_info(glue("Done preparing to plot {nrow(prepared)} SNPs."))

  # Create plot (no y_max_cap for GWASFormatter to avoid pmin issue with database)
  p = create_manhattan_plot(
    prepared,
    axis,
    y_max_cap = max(-log10(prepared$PVALUE)),
    y_axis_break = y_axis_break,
    y_axis_break_scale = y_axis_break_scale,
    label_top_n = label_top_n,
    label_strategy = label_strategy,
    label_locus_window_kb = label_locus_window_kb,
    label_pvalue_threshold = label_pvalue_threshold,
    label_nudge_y = label_nudge_y,
    force_highlight_labels = force_highlight_labels,
    label_repel_direction = label_repel_direction,
    italic_gene_labels = italic_gene_labels,
    highlight_genes = highlight_genes,
    highlight_color = highlight_color,
    label_color = label_color,
    label_size = label_size,
    highlight_label_size = highlight_label_size,
    label_segment_alpha = label_segment_alpha,
    label_max_y_nudge = label_max_y_nudge,
    base_size = base_size
  )

  # Save plot
  save_manhattan_plot(p, output, ...)
}

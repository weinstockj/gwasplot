
# Global variables to avoid NSE warnings in R CMD check
utils::globalVariables(c(
  "CHROM", "POS", "PVALUE", "CHROM_index", "chr_len", "tot", "POScum",
  "center", "aes", "gene_name", "locus_id", "is_lead", "label_text",
  "label_expr", "label_is_gene", "label_is_highlight"
))

# Helper function to create chromosome lookup table
create_chrom_lookup = function() {
  tibble::tibble(
    CHROM = c(glue::glue("chr{1:22}"), "chrX"),
    CHROM_index = 1:23
  )
}

# Helper function to prepare GWAS data for plotting
prepare_manhattan_data = function(gwas_data, chrom_lookup, lower_logp_threshold = 3.0) {
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
    dplyr::summarize(center = (max(POScum) + min(POScum)) / 2, .groups = "drop") %>% # integer overflow concern
    dplyr::arrange(CHROM_index, center)

  return(axis)
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
select_label_points = function(prepared_data, label_top_n, label_strategy = c("top_n", "lead_per_locus"),
                               label_locus_window_kb = 500, highlight_genes = NULL,
                               italic_gene_labels = TRUE, label_pvalue_threshold = 5e-8) {
  if (is.null(label_top_n) || label_top_n <= 0) {
    return(NULL)
  }

  label_strategy <- match.arg(label_strategy)

  # Only label variants that clear the significance threshold. Highlighted genes
  # bypass the gate so an explicitly requested gene is always labeled.
  has_gene_col <- "gene_name" %in% names(prepared_data)
  is_highlight_row <- if (has_gene_col && !is.null(highlight_genes)) {
    !is.na(prepared_data$gene_name) & prepared_data$gene_name %in% highlight_genes
  } else {
    rep(FALSE, nrow(prepared_data))
  }

  candidates <- prepared_data %>%
    dplyr::filter(PVALUE < label_pvalue_threshold | is_highlight_row) %>%
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

  # Top N by p-value, but always keep explicitly highlighted genes even if they
  # would otherwise rank below the cutoff.
  candidate_highlight <- if (has_gene_col && !is.null(highlight_genes)) {
    !is.na(candidates$gene_name) & candidates$gene_name %in% highlight_genes
  } else {
    rep(FALSE, nrow(candidates))
  }

  top_df <- candidates %>%
    dplyr::slice_min(PVALUE, n = label_top_n, with_ties = FALSE)
  highlight_df <- candidates[candidate_highlight, , drop = FALSE]

  label_df <- dplyr::distinct(dplyr::bind_rows(top_df, highlight_df))

  if (nrow(label_df) == 0) {
    return(NULL)
  }

  has_gene <- "gene_name" %in% names(label_df)
  gene_vals <- if (has_gene) label_df$gene_name else rep(NA_character_, nrow(label_df))
  band_vals <- derive_cytoband_labels(label_df)

  label_df <- label_df %>%
    dplyr::mutate(
      label_is_gene = !is.na(gene_vals) & gene_vals != "",
      label_text = ifelse(label_is_gene, gene_vals, band_vals),
      label_text = ifelse(is.na(label_text) | label_text == "", CHROM, label_text),
      label_is_highlight = label_is_gene & !is.null(highlight_genes) & gene_vals %in% highlight_genes,
      label_expr = ifelse(
        label_is_gene & italic_gene_labels,
        vapply(label_text, as_plotmath_label, FUN.VALUE = character(1), italic = TRUE),
        vapply(label_text, as_plotmath_label, FUN.VALUE = character(1), italic = FALSE)
      )
    )

  label_df
}

# Helper function to create the Manhattan plot
create_manhattan_plot = function(prepared_data, axis_data, pvalue_threshold = 5e-8,
                                 y_max_cap = 300, label_top_n = NULL,
                                 label_strategy = c("top_n", "lead_per_locus"),
                                 label_locus_window_kb = 500,
                                 label_pvalue_threshold = NULL,
                                 italic_gene_labels = TRUE,
                                 highlight_genes = NULL,
                                 highlight_color = "red3",
                                 label_color = "black",
                                 base_size = 7) {
  label_strategy <- match.arg(label_strategy)

  # Default the labeling cutoff to the genome-wide significance line.
  if (is.null(label_pvalue_threshold)) {
    label_pvalue_threshold <- pvalue_threshold
  }

  p = ggplot2::ggplot(prepared_data, aes(x=POScum, y=-log10(PVALUE))) +
    # Show all points
    ggrastr::rasterise(
      ggplot2::geom_point(aes(color=as.factor(CHROM_index)), alpha=0.8, size=.3),
      dev = "ragg_png",
      dpi = 300
     ) +
    ggplot2::scale_color_manual(values = c(rep(c("#7173C9", "#01035F"), 11), "#7173C9")) +
    # custom X axis:
    ggplot2::scale_x_continuous(
      label = stringr::str_replace(axis_data$CHROM, "chr", ""),
      breaks= axis_data$center, expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, max(pmin(-log10(prepared_data$PVALUE), y_max_cap)), by = 4),
      expand = ggplot2::expansion(mult = c(0, 0), add = c(0, 1))
    ) +     # remove space between plot area and x axis
    # Custom the theme:
    ggplot2::geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", color = "grey") +
    ggplot2::labs(y = expression(-log[10](pvalue))) +
    ggplot2::theme_bw(base_size = base_size, base_family = "Helvetica") +
    ggplot2::theme(
      legend.position="none",
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = base_size * 3 / 7),
      axis.text.y = element_text(size = base_size * 5 / 7),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.border = element_blank(),
      axis.line.y = element_line(size = .6)
    )

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
        label_pvalue_threshold = label_pvalue_threshold
      )

      if (!is.null(label_df) && nrow(label_df) > 0) {
        common_repel_args <- list(
          mapping = ggplot2::aes(x = POScum, y = -log10(PVALUE), label = label_expr),
          parse = TRUE,
          size = 3,
          max.overlaps = Inf,
          seed = 42,
          segment.size = 0.25,
          segment.alpha = 0.65,
          segment.color = "grey35",
          min.segment.length = 0
        )

        normal_df <- label_df %>% dplyr::filter(!label_is_highlight)
        highlight_df <- label_df %>% dplyr::filter(label_is_highlight)

        if (nrow(normal_df) > 0) {
          p <- p + do.call(
            ggrepel::geom_text_repel,
            c(list(data = normal_df, colour = label_color, inherit.aes = FALSE), common_repel_args)
          )
        }
        if (nrow(highlight_df) > 0) {
          p <- p + do.call(
            ggrepel::geom_text_repel,
            c(list(data = highlight_df, colour = highlight_color, fontface = "bold", inherit.aes = FALSE), common_repel_args)
          )
        }
      }
    } else {
      cli::cli_warn("Install the {.pkg ggrepel} package to enable Manhattan labels.")
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
#' @param label_top_n Integer. If set and a `gene_name` column is present, annotate the
#'   top N variants (by p-value) with labels using ggrepel. If `gene_name` is missing,
#'   labels fall back to cytoband-like text (e.g., `chr4p16.3`). Default is NULL (no labels).
#' @param label_strategy Label selection mode: `"top_n"` (global top N by p-value)
#'   or `"lead_per_locus"` (one lead per locus, then top N). Default `"lead_per_locus"`.
#' @param label_locus_window_kb Locus window size in kilobases used when
#'   `label_strategy = "lead_per_locus"`. Default 500.
#' @param label_pvalue_threshold P-value cutoff for labeling: only variants with
#'   `PVALUE < label_pvalue_threshold` are eligible for labels (genes in
#'   `highlight_genes` are always labeled regardless). Default NULL, which uses
#'   the genome-wide significance line (5e-8).
#' @param italic_gene_labels Logical. If TRUE, gene labels are italicized; cytoband
#'   fallback labels remain plain text. Default TRUE.
#' @param highlight_genes Optional character vector of gene symbols to highlight
#'   (e.g., novel genes). Matching labels are colored with `highlight_color`.
#' @param highlight_color Color for highlighted labels. Default `"red3"`.
#' @param label_color Color for non-highlighted labels. Default `"black"`.
#' @param ... Additional arguments passed to `ggsave`.
#' @return NULL
#' @export
manhattan = function(gwas, output, lower_logp_threshold = 3.0, label_top_n = NULL,
                     label_strategy = "lead_per_locus", label_locus_window_kb = 500,
                     label_pvalue_threshold = NULL,
                     italic_gene_labels = TRUE, highlight_genes = NULL,
                     highlight_color = "red3", label_color = "black", ...) {
  UseMethod("manhattan")
}

#' @export
manhattan.tbl_df = function(gwas, output, lower_logp_threshold = 3.0, label_top_n = NULL,
                            label_strategy = "lead_per_locus", label_locus_window_kb = 500,
                            label_pvalue_threshold = NULL,
                            italic_gene_labels = TRUE, highlight_genes = NULL,
                            highlight_color = "red3", label_color = "black", ...) {
  manhattan.data.frame(
    gwas, output, lower_logp_threshold,
    label_top_n = label_top_n,
    label_strategy = label_strategy,
    label_locus_window_kb = label_locus_window_kb,
    label_pvalue_threshold = label_pvalue_threshold,
    italic_gene_labels = italic_gene_labels,
    highlight_genes = highlight_genes,
    highlight_color = highlight_color,
    label_color = label_color,
    ...
  )
}

#' @export
manhattan.data.frame = function(gwas, output, lower_logp_threshold = 3.0, label_top_n = NULL,
                                label_strategy = "lead_per_locus", label_locus_window_kb = 500,
                                label_pvalue_threshold = NULL,
                                italic_gene_labels = TRUE, highlight_genes = NULL,
                                highlight_color = "red3", label_color = "black",
                                base_size = 7, ...) {
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
    prepared, axis,
    label_top_n = label_top_n,
    label_strategy = label_strategy,
    label_locus_window_kb = label_locus_window_kb,
    label_pvalue_threshold = label_pvalue_threshold,
    italic_gene_labels = italic_gene_labels,
    highlight_genes = highlight_genes,
    highlight_color = highlight_color,
    label_color = label_color,
    base_size = base_size
  )

  # Save plot
  save_manhattan_plot(p, output, ...)
}

#' @export
manhattan.GWASFormatter = function(gwas, output, lower_logp_threshold = 3.0, label_top_n = NULL,
                                   label_strategy = "lead_per_locus", label_locus_window_kb = 500,
                                   label_pvalue_threshold = NULL,
                                   italic_gene_labels = TRUE, highlight_genes = NULL,
                                   highlight_color = "red3", label_color = "black",
                                   base_size = 7, ...) {
  log_info("Now preparing to plot")

  # Create chromosome lookup and copy to database connection
  chrom_lookup = copy_to_if_missing(gwas$con, create_chrom_lookup(), "chrom_lookup")
  
  # Prepare data
  prepared = prepare_manhattan_data(gwas$data, chrom_lookup, lower_logp_threshold)
  
  # Collect prepared data from database
  prepared = prepared %>%
    dplyr::collect(.)

  # Calculate axis positions after collecting the prepared points once
  axis = calculate_axis_positions(prepared)

  log_info(glue("Done preparing to plot {nrow(prepared)} SNPs."))

  # Create plot (no y_max_cap for GWASFormatter to avoid pmin issue with database)
  p = create_manhattan_plot(prepared, axis, y_max_cap = max(-log10(prepared$PVALUE)),
                            label_top_n = label_top_n,
                            label_strategy = label_strategy,
                            label_locus_window_kb = label_locus_window_kb,
                            label_pvalue_threshold = label_pvalue_threshold,
                            italic_gene_labels = italic_gene_labels,
                            highlight_genes = highlight_genes,
                            highlight_color = highlight_color,
                            label_color = label_color,
                            base_size = base_size)

  # Save plot
  save_manhattan_plot(p, output, ...)
}

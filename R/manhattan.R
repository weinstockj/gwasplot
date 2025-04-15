manhattan = function(gwas, output_prefix, lower_logp_threshold = 3.0) {

  chrom_lookup = tibble::tibble(
    CHROM = c(glue::glue("chr{1:22}"), "chrX"),
    CHROM_index = 1:23
  ) %>%
    dplyr::copy_to(gwas$con, ., "chrom_lookup", overwrite = TRUE)

  scaling = 1e8
  pvalue_threshold = 5e-8

  log_info("Now preparing to plot")

  prepared = gwas$data %>%
    dplyr::select(CHROM, POS, PVALUE) %>%
    dplyr::inner_join(chrom_lookup, by = "CHROM") %>%
    # Compute chromosome size
    dplyr::group_by(CHROM_index, CHROM) %>%
    dplyr::summarise(chr_len = (max(POS) - min(POS)) / scaling) %>%
    dplyr::ungroup(.) %>%
    # Calculate cumulative start position of each chromosome
    dplyr::arrange(CHROM_index) %>%
    dplyr::mutate(tot = cumsum(chr_len) - chr_len) %>%
    # Add this info to the initial dataset
    dplyr::right_join(gwas$data %>% dplyr::select(CHROM, POS, PVALUE), by = "CHROM") %>%
    #     # Add a cumulative position of each SNP
    dplyr::arrange(CHROM_index, POS) %>%
    dplyr::mutate(POScum = POS / scaling + tot) %>%
    # Filter SNP to make the plot lighter
    dplyr::filter(-log10(PVALUE) > lower_logp_threshold)


  axis = prepared %>%
    dplyr::group_by(CHROM_index, CHROM) %>%
    dplyr::summarize(center = (max(POScum) + min(POScum)) / 2) %>% # integer overflow concern
    dplyr::ungroup(.) %>%
    dplyr::arrange(CHROM_index, center) %>%
    dplyr::collect(.)

  prepared = prepared %>%
    dplyr::collect(.)

  log_info(glue("Done preparing to plot {nrow(prepared)} SNPs."))

  p = ggplot2::ggplot(prepared, aes(x=POScum, y=-log10(PVALUE))) +
    # Show all points
    # ggrastr::geom_point_rast(aes(color=as.factor(CHROM_index)), alpha=0.8, size=.3, raster.dpi = 300) +
    ggplot2::geom_point(aes(color=as.factor(CHROM_index)), alpha=0.8, size=.3) +
    ggplot2::scale_color_manual(values = c(rep(c("#7173C9", "#01035F"), 11), "#7173C9")) +
    # custom X axis:
    ggplot2::scale_x_continuous(
      label = stringr::str_replace(axis$CHROM, "chr", ""),
      breaks= axis$center, expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, max(-log10(prepared$PVALUE)), by = 4),
      expand = ggplot2::expansion(mult = c(0, 0), add = c(0, 1))
    ) +     # remove space between plot area and x axis
    # Custom the theme:
    # ggplot2::geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "grey") +
    ggplot2::geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", color = "grey") +
    ggplot2::labs(y = expression(-log[10](pvalue))) +
    ggplot2::theme_bw(base_size = 7, base_family = "Helvetica") +
    ggplot2::theme(
      legend.position="none",
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 3),
      axis.text.y = element_text(size = 5),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.border = element_blank(),
      axis.line.y = element_line(size = .6)
    )


  log_info("Now rendering")
  # fname_pdf = file.path(output_prefix, "manhattan.pdf")
  fname_png = glue("{output_prefix}_manhattan.png")
  # ggplot2::ggsave(fname_pdf, p, units = "mm", width = 89, height = 70, dpi = 300)
  ggplot2::ggsave(
    fname_png,
    p,
    device = ragg::agg_png,
    units = "mm",
    width = 89,
    height = 70,
    dpi = 300,
    bg = "white"
  )

  log_info("done plotting.")
}

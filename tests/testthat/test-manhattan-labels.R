test_that("calculate_label_layout keeps adaptive nudges bounded", {
  label_df <- tibble::tibble(
    CHROM = "chr1",
    POS = c(1L, 2L, 3L),
    PVALUE = c(1e-12, 1e-11, 1e-10),
    CHROM_index = 1L,
    chr_len = 1,
    tot = 0,
    POScum = c(1.00, 1.01, 1.02),
    label_expr = c("'A'", "'B'", "'C'"),
    label_is_highlight = FALSE
  )

  prepared_data <- dplyr::bind_rows(
    label_df,
    tibble::tibble(
      CHROM = "chr1",
      POS = 4L,
      PVALUE = 1e-4,
      CHROM_index = 1L,
      chr_len = 1,
      tot = 0,
      POScum = 3,
      label_expr = "'D'",
      label_is_highlight = FALSE
    )
  )

  laid_out <- calculate_label_layout(
    label_df,
    prepared_data,
    label_max_y_nudge = 0.8
  )

  expect_true(all(laid_out$label_nudge_y <= 0.8))
  expect_true(all(laid_out$label_nudge_y >= 0))
  expect_equal(nrow(laid_out), nrow(label_df))
})

test_that("calculate_label_layout stacks nearby labels more than isolated labels", {
  label_df <- tibble::tibble(
    CHROM = "chr1",
    POS = 1:3,
    PVALUE = c(1e-12, 1e-11, 1e-10),
    CHROM_index = 1L,
    chr_len = 1,
    tot = 0,
    POScum = c(1.00, 1.01, 2.50),
    label_expr = c("'A'", "'B'", "'C'"),
    label_is_highlight = FALSE
  )

  laid_out <- calculate_label_layout(
    label_df,
    label_df,
    label_max_y_nudge = 3
  )

  close_labels <- laid_out[laid_out$POScum < 2, ]
  isolated_label <- laid_out[laid_out$POScum > 2, ]

  expect_gt(max(close_labels$label_nudge_y), min(close_labels$label_nudge_y))
  expect_equal(isolated_label$label_nudge_y, min(laid_out$label_nudge_y))
})

test_that("highlight genes are labeled by default even when weak", {
  prepared_data <- tibble::tibble(
    CHROM = "chr1",
    POS = 1:3,
    PVALUE = c(1e-12, 1e-11, 1e-4),
    CHROM_index = 1L,
    chr_len = 1,
    tot = 0,
    POScum = c(1, 2, 3),
    gene_name = c("STRONG1", "STRONG2", "WEAK")
  )

  label_df <- select_label_points(
    prepared_data,
    label_top_n = 1,
    label_strategy = "top_n",
    highlight_genes = "WEAK",
    label_pvalue_threshold = 5e-8
  )

  expect_true("WEAK" %in% label_df$gene_name)
  expect_true(label_df$label_is_highlight[label_df$gene_name == "WEAK"])
})

test_that("highlight label forcing can be disabled", {
  prepared_data <- tibble::tibble(
    CHROM = "chr1",
    POS = 1:3,
    PVALUE = c(1e-12, 1e-11, 1e-4),
    CHROM_index = 1L,
    chr_len = 1,
    tot = 0,
    POScum = c(1, 2, 3),
    gene_name = c("STRONG1", "STRONG2", "WEAK")
  )

  label_df <- select_label_points(
    prepared_data,
    label_top_n = 1,
    label_strategy = "top_n",
    highlight_genes = "WEAK",
    label_pvalue_threshold = 5e-8,
    force_highlight_labels = FALSE
  )

  expect_false("WEAK" %in% label_df$gene_name)
})

test_that("forced highlighted genes are labeled once at their strongest variant", {
  prepared_data <- tibble::tibble(
    CHROM = "chr1",
    POS = 1:4,
    PVALUE = c(1e-12, 1e-4, 1e-9, 1e-6),
    CHROM_index = 1L,
    chr_len = 1,
    tot = 0,
    POScum = 1:4,
    gene_name = c("TOP", "HILITE", "HILITE", "HILITE")
  )

  label_df <- select_label_points(
    prepared_data,
    label_top_n = 1,
    label_strategy = "top_n",
    highlight_genes = "HILITE",
    label_pvalue_threshold = 5e-8
  )

  hilite <- label_df[label_df$gene_name == "HILITE", ]
  expect_equal(nrow(hilite), 1)
  expect_equal(hilite$PVALUE, 1e-9)
})

test_that("forced highlighted genes survive lead-locus thinning", {
  prepared_data <- tibble::tibble(
    CHROM = "chr22",
    POS = c(1000L, 1200L, 1400L),
    PVALUE = c(1e-20, 1e-6, 1e-5),
    CHROM_index = 22L,
    chr_len = 1,
    tot = 0,
    POScum = c(1, 1.000002, 1.000004),
    gene_name = c("LEAD", "VSIG4", "VSIG4")
  )

  label_df <- select_label_points(
    prepared_data,
    label_top_n = 1,
    label_strategy = "lead_per_locus",
    label_locus_window_kb = 500,
    highlight_genes = "VSIG4",
    label_pvalue_threshold = 5e-8
  )

  highlight <- label_df[label_df$gene_name == "VSIG4", ]
  expect_equal(nrow(highlight), 1)
  expect_equal(highlight$PVALUE, 1e-6)
})

test_that("transform_manhattan_y compresses only values above the break", {
  y <- transform_manhattan_y(
    c(10, 50, 70, 110),
    y_axis_break = 50,
    y_axis_break_scale = 0.2
  )

  expect_equal(y, c(10, 50, 54, 62))
})

test_that("transform_manhattan_y preserves every finite point order", {
  logp <- c(4, 12, 16, 24, 40, 80, 112)
  transformed <- transform_manhattan_y(
    logp,
    y_axis_break = 16,
    y_axis_break_scale = 0.2
  )

  expect_equal(length(transformed), length(logp))
  expect_false(anyNA(transformed))
  expect_true(all(diff(transformed) > 0))
})

test_that("calculate_y_axis keeps original labels for compressed ticks", {
  prepared_data <- tibble::tibble(
    LOGP = c(4, 12, 48, 70, 110),
    LOGP_plot = transform_manhattan_y(
      LOGP,
      y_axis_break = 50,
      y_axis_break_scale = 0.2
    )
  )

  axis <- calculate_y_axis(
    prepared_data,
    y_axis_break = 50,
    y_axis_break_scale = 0.2
  )

  expect_true(max(axis$breaks) < max(axis$labels))
  expect_true(max(axis$labels) >= 110)
  expect_equal(min(axis$labels), 4)
  expect_equal(axis$limits[[1]], 4)
})

test_that("calculate_y_axis starts uncompressed plots at observed minimum", {
  prepared_data <- tibble::tibble(LOGP = c(3.2, 5, 9.5), LOGP_plot = LOGP)

  axis <- calculate_y_axis(prepared_data)

  expect_equal(axis$limits, c(3.2, 9.5))
  expect_equal(axis$labels[[1]], 3.2)
  expect_equal(axis$breaks[[1]], 3.2)
})

test_that("compressed Manhattan plot includes a break slash layer", {
  prepared_data <- tibble::tibble(
    CHROM = "chr1",
    POS = 1:4,
    PVALUE = c(1e-110, 1e-80, 1e-30, 1e-8),
    CHROM_index = 1L,
    chr_len = 1,
    tot = 0,
    POScum = 1:4
  )
  axis_data <- tibble::tibble(CHROM = "chr1", CHROM_index = 1L, center = 2.5)

  p <- create_manhattan_plot(
    prepared_data,
    axis_data,
    y_axis_break = 50,
    y_axis_break_scale = 0.2
  )
  built <- ggplot2::ggplot_build(p)

  slash_layer <- built$data[[length(built$data)]]
  axis_layer <- built$data[[length(built$data) - 1]]
  axis_x <- unique(axis_layer$x)

  expect_equal(nrow(axis_layer), 1)
  expect_equal(axis_layer$x, axis_layer$xend)
  expect_equal(axis_layer$y, min(transform_manhattan_y(
    -log10(prepared_data$PVALUE),
    y_axis_break = 50,
    y_axis_break_scale = 0.2
  )))
  expect_equal(nrow(slash_layer), 2)
  expect_true(all(c("x", "xend", "y", "yend") %in% names(slash_layer)))
  expect_true(all(slash_layer$x < axis_x))
  expect_true(all(slash_layer$xend > axis_x))
})

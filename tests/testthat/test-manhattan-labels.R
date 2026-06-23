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

test_that("transform_manhattan_y compresses only values above the break", {
  y <- transform_manhattan_y(
    c(10, 50, 70, 110),
    y_axis_break = 50,
    y_axis_break_scale = 0.2
  )

  expect_equal(y, c(10, 50, 54, 62))
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
  expect_equal(nrow(slash_layer), 2)
  expect_true(all(c("x", "xend", "y", "yend") %in% names(slash_layer)))
})

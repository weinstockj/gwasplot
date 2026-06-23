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

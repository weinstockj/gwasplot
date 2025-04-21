#' @title QQ plot for p-values
#' @param pvalues A numeric vector of p-values or a list of numeric vectors of p-values.
#' @param should.thin Logical indicating whether to thin the points for plotting. Default is TRUE.
#' @param thin.obs.places Number of decimal places to round the observed p-values for thinning. Default is 2.
#' @param thin.exp.places Number of decimal places to round the expected p-values for thinning. Default is 2.
#' @export
qqunif_plot = function(pvalues,
                       should.thin = TRUE,
                       thin.obs.places = 2,
                       thin.exp.places = 2,
                       xlab = expression(paste("Expected (", -log[10], " p-value)")),
                       ylab = expression(paste("Observed (", -log[10], " p-value)")),
                       draw.conf = TRUE,
                       conf.points = 1000,
                       conf.col = "lightgray",
                       conf.alpha = .05,
                       already.transformed = FALSE,
                       point.size = 1.5,
                       point.alpha = 0.7,
                       point.color = "#3366FF",
                       ...) {
  # Error checking
  if (length(pvalues) == 0)
    stop("pvalue vector is empty, can't draw plot")
  if (!(class(pvalues) == "numeric" ||
        (class(pvalues) == "list" &&
         all(sapply(pvalues, class) == "numeric"))))
    stop("pvalue vector is not numeric, can't draw plot")
  if (any(is.na(unlist(pvalues))))
    stop("pvalue vector contains NA values, can't draw plot")
  if (already.transformed == FALSE) {
    if (any(unlist(pvalues) == 0))
      stop("pvalue vector contains zeros, can't draw plot")
  } else {
    if (any(unlist(pvalues) < 0))
      stop("-log10 pvalue vector contains negative values, can't draw plot")
  }

  # Process data
  plot_data <- NULL
  grp <- NULL
  n <- 1
  
  if (is.list(pvalues)) {
    # Handle list of p-values (multiple datasets)
    nn <- sapply(pvalues, length)
    rs <- cumsum(nn)
    re <- rs - nn + 1
    n <- min(nn)
    
    if (!is.null(names(pvalues))) {
      grp = factor(rep(names(pvalues), nn), levels = names(pvalues))
    } else {
      grp = factor(rep(1:length(pvalues), nn))
    }
    
    # Process each dataset
    plot_data <- data.frame(
      observed = numeric(sum(nn)),
      expected = numeric(sum(nn)),
      group = grp
    )
    
    for (i in 1:length(pvalues)) {
      if (!already.transformed) {
        plot_data$observed[rs[i]:re[i]] <- -log10(pvalues[[i]])
        plot_data$expected[rs[i]:re[i]] <- -log10((rank(pvalues[[i]], ties.method = "first") - 0.5) / nn[i])
      } else {
        plot_data$observed[rs[i]:re[i]] <- pvalues[[i]]
        plot_data$expected[rs[i]:re[i]] <- -log10((nn[i] + 1 - rank(pvalues[[i]], ties.method = "first") - 0.5) / (nn[i] + 1))
      }
    }
  } else {
    # Handle single vector of p-values
    n <- length(pvalues) + 1
    
    plot_data <- data.frame(
      expected = numeric(length(pvalues)),
      observed = numeric(length(pvalues))
    )
    
    if (!already.transformed) {
      plot_data$expected <- -log10((rank(pvalues, ties.method = "first") - 0.5) / n)
      plot_data$observed <- -log10(pvalues)
    } else {
      plot_data$expected <- -log10((n - rank(pvalues, ties.method = "first") - 0.5) / n)
      plot_data$observed <- pvalues
    }
  }
  
  # Thin points if requested
  if (should.thin) {
    if (!is.null(grp)) {
      plot_data <- dplyr::distinct(
        dplyr::mutate(plot_data,
          observed = round(observed, thin.obs.places),
          expected = round(expected, thin.exp.places)
        )
      )
    } else {
      plot_data <- dplyr::distinct(
        dplyr::mutate(plot_data,
          observed = round(observed, thin.obs.places),
          expected = round(expected, thin.exp.places)
        )
      )
    }
  }
  
  # Calculate confidence intervals if needed
  if (draw.conf) {
    conf.points <- min(conf.points, n - 1)
    conf_data <- qqconf(n, conf.points, conf.alpha)
    
    # Convert matrix to data frame for ggplot
    conf_df <- data.frame(
      expected = conf_data[,1],
      observed = conf_data[,2]
    )
  }
  
  # Create the plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = expected, y = observed))
  
  # Add confidence interval
  if (draw.conf) {
    p <- p + ggplot2::geom_polygon(data = conf_df, 
                                  ggplot2::aes(x = expected, y = observed), 
                                  fill = conf.col, 
                                  alpha = 0.5)
  }
  
  # Add points
  if (!is.null(grp)) {
    p <- p + ggplot2::geom_point(ggplot2::aes(color = group), 
                                alpha = point.alpha,
                                size = point.size)
  } else {
    p <- p + ggplot2::geom_point(color = point.color, 
                                alpha = point.alpha,
                                size = point.size)
  }
  
  # Add diagonal line
  p <- p + ggplot2::geom_abline(slope = 1, intercept = 0, 
                              color = "darkgray", 
                              linetype = "dashed")
  
  # Adjust theme and appearance
  p <- p + ggplot2::theme_bw(base_size = 12, base_family = "Helvetica") +
    ggplot2::coord_equal() +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "lightgray", linetype = "dotted"),
      panel.border = ggplot2::element_rect(color = "black", fill = NA),
      axis.line = ggplot2::element_line(color = "black"),
      axis.text = ggplot2::element_text(color = "black"),
      axis.ticks = ggplot2::element_line(color = "black"),
      legend.position = if (is.null(grp)) "none" else "right"
    )
  
  # Return the plot
  return(p)
}

#' Save a QQ plot of GWAS p-values
#' 
#' @param gwas A gwas object containing the p-values to plot.
#' @param output_prefix The prefix for the output file name.
#' @param width Width of the output figure in inches.
#' @param height Height of the output figure in inches.
#' @param dpi DPI of the output figure.
#' @return NULL
#' @export
qqunif_plot_save = function(gwas, output_prefix, width = 5, height = 5, dpi = 300) {
  log_info("Creating QQ plot")
  fname <- glue::glue("{output_prefix}_qqplot.png")
  
  # Sample data to avoid memory issues with large datasets
  con <- db_connect() 
  tbl = dplyr::tbl(con, "summary_stats")
  
  p_values <- dplyr::pull(tbl, PVALUE) %>%
    dplyr::collect(.)
  
  # Create and save the plot
  p <- qqunif_plot(p_values)
  
  ggplot2::ggsave(
    filename = fname,
    plot = p,
    device = ragg::agg_png,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
  
  log_info(glue::glue("QQ plot saved to {fname}"))
  
  return(invisible(NULL))
}

#' @title QQ plot for p-values
#' @param pvalues A numeric vector of p-values or a list of numeric vectors of p-values.
#' @param should.thin Logical indicating whether to thin the points for plotting. Default is TRUE.
#' @param thin.obs.places Number of decimal places to round the observed p-values for thinning. Default is 2.
#' @param thin.exp.places Number of decimal places to round the expected p-values for thinning. Default is 2.
#' @export
qqunif_plot = function(pvalues,
                       should.thin = T,
                       thin.obs.places = 2,
                       thin.exp.places = 2,
                       xlab = expression(paste("Expected (", -log[10], " p-value)")),
                       ylab = expression(paste("Observed (", -log[10], " p-value)")),
                       draw.conf = TRUE,
                       conf.points = 1000,
                       conf.col = "lightgray",
                       conf.alpha = .05,
                       already.transformed = FALSE,
                       pch = 20,
                       aspect = "iso",
                       prepanel = prepanel.qqunif,
                       par.settings = list(superpose.symbol = list(pch =
                                                                     pch)),
                       ...) {
  #error checking
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


  grp <- NULL
  n <- 1
  exp.x <- c()
  if (is.list(pvalues)) {
    nn <- sapply(pvalues, length)
    rs <- cumsum(nn)
    re <- rs - nn + 1
    n <- min(nn)
    if (!is.null(names(pvalues))) {
      grp = factor(rep(names(pvalues), nn), levels = names(pvalues))
      names(pvalues) <- NULL
    } else {
      grp = factor(rep(1:length(pvalues), nn))
    }
    pvo <- pvalues
    pvalues <- numeric(sum(nn))
    exp.x <- numeric(sum(nn))
    for (i in 1:length(pvo)) {
      if (!already.transformed) {
        pvalues[rs[i]:re[i]] <- -log10(pvo[[i]])
        exp.x[rs[i]:re[i]] <- -log10((rank(pvo[[i]], ties.method = "first") -
                                        .5) / nn[i])
      } else {
        pvalues[rs[i]:re[i]] <- pvo[[i]]
        exp.x[rs[i]:re[i]] <- -log10((nn[i] + 1 - rank(pvo[[i]], ties.method =
                                                         "first") - .5) / (nn[i] + 1))
      }
    }
  } else {
    n <- length(pvalues) + 1
    if (!already.transformed) {
      exp.x <- -log10((rank(pvalues, ties.method = "first") - .5) / n)
      pvalues <- -log10(pvalues)
    } else {
      exp.x <- -log10((n - rank(pvalues, ties.method = "first") - .5) / n)
    }
  }


  #this is a helper function to draw the confidence interval
  panel.qqconf <- function(n,
                           conf.points = 1000,
                           conf.col = "gray",
                           conf.alpha = .05,
                           ...) {
    conf.points = min(conf.points, n - 1)

    # mpts <- matrix(nrow = conf.points * 2, ncol = 2)
    # for (i in seq(from = 1, to = conf.points)) {
    #   mpts[i, 1] <- -log10((i - .5) / n)
    #   mpts[i, 2] <- -log10(qbeta(1 - conf.alpha / 2, i, n - i))
    #   mpts[conf.points * 2 + 1 - i, 1] <- -log10((i - .5) / n)
    #   mpts[conf.points * 2 + 1 - i, 2] <- -log10(qbeta(conf.alpha / 2, i, n -
    #                                                      i))
    # }
    mpts = qqconf(n, conf.points, conf.alpha)

    grid::grid.polygon(
      x = mpts[, 1],
      y = mpts[, 2],
      gp = grid::gpar(fill = conf.col, lty = 0),
      default.units = "native"
    )
  }

  #reduce number of points to plot
  if (should.thin == T) {
    if (!is.null(grp)) {
      thin <- dplyr::distinct(data.frame(
        pvalues = round(pvalues, thin.obs.places),
        exp.x = round(exp.x, thin.exp.places),
        grp = grp
      ))
      grp = thin$grp
    } else {
      thin <- dplyr::distinct(data.frame(
        pvalues = round(pvalues, thin.obs.places),
        exp.x = round(exp.x, thin.exp.places)
      ))
    }
    pvalues <- thin$pvalues
    exp.x <- thin$exp.x
  }
  gc()

  prepanel.qqunif = function(x, y, ...) {
    A = list()
    A$xlim = range(x, y) * 1.02
    A$xlim[1] = 0
    A$ylim = A$xlim
    return(A)
  }

  #draw the plot
  lattice::xyplot(
    pvalues ~ exp.x,
    groups = grp,
    xlab = xlab,
    ylab = ylab,
    aspect = aspect,
    prepanel = prepanel,
    scales = list(axs = "i"),
    pch = pch,
    panel = function(x, y, ...) {
      if (draw.conf) {
        panel.qqconf(
          n,
          conf.points = conf.points,
          conf.col = conf.col,
          conf.alpha = conf.alpha
        )
      }

      lattice::panel.xyplot(x, y, ...)

      lattice::panel.abline(0, 1)

    },
    par.settings = par.settings,
    ...
  )

}

qqunif_plot_save = function(gwas, output_prefix) {
  fname = glue("{output_prefix}_qqplot.png")
  ragg::agg_png(file = fname,
      width = 5,
      height = 5, units = "in")
  # png(file = fname, type = "cairo")
  # lattice::trellis.device(file = fname)
  print(qqunif_plot(dplyr::pull(compute_sample(gwas, 1e7), PVALUE)))
  dev.off()
}

##' Bin summary by covariate
##'
##' Drawing a bin summary plot of covariates given a `multinethist` object with a user-specified order.
##'
##' @param object a `multinethist` object from [multinethist()].
##' @param covariate a vector for univariate covariate. If it is a factor, a stacked bar chart is drawn. If it is numeric, a violin plot is drawn.
##' @param idx_order A numeric vector for index label order, which must be a permutation of `object$cluster`. If `NA`, it uses `1:max(object$cluster)`.
##' @param main title of summary plot. If NA, the plot has no title.
##' @param ylab label of y-axis. If NA, y-axis label is "covariate"
##' @param legend_title title of legend. If NA, the legend title is "covariate"
##' @param stat variables pass to [ggplot2::geom_bar()]. Only used for a factor covariate.
##' @param position variables pass to [ggplot2::geom_bar()]. Only used for a factor covariate.
##' @param ... currently unused.
##' @details
##' When `covariate` is a factor, a stacked bar chart is drawn with bins ordered by `idx_order`.
##' When `covariate` is numeric, a violin plot is drawn.
##' @returns
##' A `ggplot` object. Printed as a side effect. Returns the plot invisibly.
##' @importFrom ggplot2 ggplot aes geom_bar geom_violin labs ggtitle ylab theme element_rect element_blank
##' @importFrom rlang .data
##' @examples
##' {
##' set.seed(42)
##' data(polblog)
##' nethist_polblog <- multinethist(polblog)
##' x <- factor(c(rep("Liberal", 586), rep("Conservative", 638)))
##' covariate_plot(nethist_polblog, x)
##' }
##' @export
covariate_plot <- function(object, covariate,
                           idx_order = 1:max(object$cluster),
                           main = NA,
                           ylab = NA,
                           legend_title = NA,
                           stat = "count", position = "stack") {
  UseMethod("covariate_plot")
}

##' @exportS3Method
covariate_plot.nethist <- function(object, covariate,
                                   idx_order = 1:max(object$cluster),
                                   main = NA,
                                   ylab = NA,
                                   legend_title = NA,
                                   stat = "count", position = "stack") {
  return(covariate_plot.multinethist(object, covariate, idx_order, main, ylab,
                                     legend_title, stat, position))
}

##' @exportS3Method
covariate_plot.multinethist <- function(object, covariate,
                                        idx_order = 1:max(object$cluster),
                                        main = NA,
                                        ylab = NA,
                                        legend_title = NA,
                                        stat = "count", position = "stack") {
  k <- max(object$cluster)
  if (!.is_valid_order(idx_order, 1:k)) {
    warning(paste0("idx_order is invalid. Set idx_order = 1:", k))
    idx_order <- 1:k
  }
  if (length(covariate) != length(object$cluster))
    stop("Length of covariate must equal the number of vertices (length of object$cluster).")
  
  df <- data.frame(cluster   = factor(object$cluster, levels = idx_order),
                   covariate = covariate)
  if (is.factor(covariate)) {
    p <- summary_covariate_factor(df, legend_title, stat, position)
  } else if (is.vector(covariate) & is.numeric(covariate)) {
    p <- summary_covariate_numeric(df)
  } else {
    stop("covariate must be either a numeric vector or a factor.")
  }
  
  if (!is.na(main)) p <- p + ggplot2::ggtitle(main)
  if (!is.na(ylab)) p <- p + ggplot2::ylab(ylab)
  p <- p + ggplot2::theme(
    panel.background      = ggplot2::element_rect(fill = "transparent"),
    plot.background       = ggplot2::element_rect(fill = "transparent", color = NA),
    panel.grid.major      = ggplot2::element_blank(),
    panel.grid.minor      = ggplot2::element_blank(),
    legend.background     = ggplot2::element_rect(fill = "transparent"),
    legend.box.background = ggplot2::element_rect(fill = "transparent")
  )
  print(p)
  invisible(p)
}

summary_covariate_factor <- function(df, legend_title, stat, position) {
  p <- ggplot2::ggplot(
    data    = df,
    mapping = ggplot2::aes(x = .data$cluster, fill = .data$covariate)
  ) + ggplot2::geom_bar(position = position, stat = stat)
  if (!is.na(legend_title)) p <- p + ggplot2::labs(fill = legend_title)
  p
}

summary_covariate_numeric <- function(df) {
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$cluster, y = .data$covariate)) +
    ggplot2::geom_violin()
  p
}

##' @rdname covariate_plot
##' @export
summary_plot <- function(object, covariate, ...) {
  .Deprecated("covariate_plot")
  covariate_plot(object, covariate, ...)
}

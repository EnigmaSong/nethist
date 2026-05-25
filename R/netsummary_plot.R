##' Network summary plots
##'
##' Draw a network summary plot proposed by Maugis et al. (2017). To count k-cycles, Alon et al. (1997) is used.
##'
##' @param A an adjacency matrix, igraph object, or network object to draw a network summary plot. It must be an undirected and simple graph.
##' @param subsample_sizes a numeric vector of vertex subsample sizes. If `NA`, the subsample size is selected automatically.
##' @param max_cycle_order an integer value of the maximum cycle size. Must be `>=3` and `<=7`.
##' @param n_rep an integer value of subsampling replication. If `NA`, `n_rep` is automatically selected by `alpha`.
##' @param n_subsample_sizes number of different subsample sizes for automatic selection. It is only used when `subsample_sizes = NA`.
##' @param alpha a pre-specified level used in determining `n_rep` and `subsample_sizes` when they are not specified. It must be in (0,1). Default is 0.05. Smaller `alpha` gives larger `n_rep` and `subsample_sizes`.
##' @param y_max Upper limit of y-axis of the plot. Must be 0 < `y_max` <= 1. If `NA`, the upper limit is automatically selected.
##' @param save_plot A logical indicating whether to save the generated figure. If `TRUE`, the plot is saved via [ggplot2::ggsave()] using the specified file name. Otherwise, the plot is displayed.
##' @param filename file name to save the generated figure.
##' @param width a numeric value of the width of the generated figure in inch. It is only used when `save_plot = TRUE`.
##' @param height a numeric value of the height of the generated figure in inch. It is only used when `save_plot = TRUE`.
##' @param max_subsample_size integer. Upper bound on the automatically selected subsample size. Larger values improve statistical accuracy but increase computation time. Default is 250. 
##' @param ... `r lifecycle::badge("deprecated")` Pass \code{R}, \code{Ns}, \code{y.max}, or \code{save.plot} via the renamed arguments \code{n_rep}, \code{n_subsample_sizes}, \code{y_max}, and \code{save_plot} instead.
##' @returns
##' A `ggplot` object. Printed as a side effect when `save_plot = FALSE`.
##' Returns the plot invisibly when `save_plot = TRUE`.
##' @details
##' Vertex sampling is done by simple random sampling without replacement.
##'
##' The automatically selected subsample size is capped at `max_subsample_size`
##' to limit computation time.
##'
##' Each violin shows the distribution of the subsampled statistic, and a dot
##' marks the mean.
##'
##' The following input classes are supported: [base::matrix],
##' [Matrix::dgCMatrix-class], [igraph::igraph], [network::network].
##' @references Maugis et al. (2017). Topology reveals universal features for network comparison. arXiv: 1705.05677
##' @references Alon et al. (1997). Finding and counting given length cycles. Algorithmica 17, 209–223 (1997). https://doi.org/10.1007/BF02523189
##' @examples
##' {
##' set.seed(2022)
##' #Generating Erdos-Renyi graph
##' n <- 400
##' #igraph object
##' A <- igraph::sample_gnp(n, 0.05)
##' netsummary_plot(A)
##'
##' #sparse adjacency matrix
##' A2 <- igraph::as_adjacency_matrix(A)
##' netsummary_plot(A2)
##'
##' #dense adjacency matrix
##' A2 <- igraph::as_adjacency_matrix(A, sparse = FALSE)
##' netsummary_plot(A2)
##'
##' #user-specified n_rep and subsample_sizes
##' netsummary_plot(A, n_rep = 500, subsample_sizes = 150)
##'
##' #user-specified alpha
##' netsummary_plot(A, alpha = 0.1)
##'
##' #network object
##' A3 <- network::as.network(igraph::as_adjacency_matrix(A, sparse = FALSE))
##' netsummary_plot(A3)
##'
##' #user-specified max_subsample_size
##' netsummary_plot(A, max_subsample_size = 100)
##'
##' #saving the plot with user-specified file name
##' \dontrun{
##' netsummary_plot(A, save_plot = TRUE, filename = "myfig.pdf")
##' }
##' }
##' @importFrom ggtext element_markdown
##' @importFrom ggplot2 ggplot aes geom_violin stat_summary ylim ylab scale_x_discrete theme ggsave rel
##' @importFrom rlang .data
##' @export
##'
netsummary_plot <- function(A,
                            subsample_sizes = NA,
                            max_cycle_order = 4,
                            n_rep = NA,
                            n_subsample_sizes = 11,
                            alpha = 0.05,
                            y_max = NA,
                            save_plot = FALSE,
                            filename = "myplot.pdf",
                            width = 7, height = 5,
                            max_subsample_size = 250,
                            ...) {
  UseMethod("netsummary_plot")
}

##' @exportS3Method
netsummary_plot.igraph <- function(A, ...) {
  netsummary_plot.default(
    igraph::as_adjacency_matrix(A, sparse = FALSE), ...)
}

##' @exportS3Method
netsummary_plot.matrix <- function(A, ...) {
  netsummary_plot.default(A, ...)
}

##' @exportS3Method
netsummary_plot.dgCMatrix <- function(A, ...) {
  netsummary_plot.default(as.matrix(A), ...)
}

##' @exportS3Method
netsummary_plot.network <- function(A, ...) {
  netsummary_plot.default(
    as.matrix(A, matrix.type = "adjacency"), ...)
}

##' @exportS3Method
netsummary_plot.default <- function(A,
                                    subsample_sizes = NA,
                                    max_cycle_order = 4,
                                    n_rep = NA,
                                    n_subsample_sizes = 11,
                                    alpha = 0.05,
                                    y_max = NA,
                                    save_plot = FALSE,
                                    filename = "myplot.pdf",
                                    width = 7, height = 5,
                                    max_subsample_size = 250,
                                    ...) {
  # Handle deprecated argument names
  dots <- list(...)
  old_map <- c(R          = "n_rep",
               Ns         = "n_subsample_sizes",
               y.max      = "y_max",
               save.plot  = "save_plot")
  hits <- intersect(names(dots), names(old_map))
  if (length(hits)) {
    warning(paste0(
      "Argument(s) ", paste(hits, collapse = ", "),
      " are deprecated. Use ",
      paste(old_map[hits], collapse = ", "), " instead."),
      call. = FALSE)
    for (nm in hits) assign(old_map[nm], dots[[nm]])
  }

  if (!is.matrix(A)) stop("A must be a matrix.")
  if (nrow(A) != ncol(A)) stop("A must be a square matrix.")
  if (any(diag(A) != 0))
    stop("A has self-loops. Network A must be a simple graph (no self-loops).")
  if (!isSymmetric(unname(A)))
    stop("A is not symmetric. Network A must be an undirected graph.")
  if (any(A[A != 0] != 1))
    stop("A is not a simple graph. All non-zero entries must be 1 (binary adjacency matrix).")

  if ((max_cycle_order < 3) | (max_cycle_order %% 1 != 0)) {
    stop("order_cycle must be >= 3 integers.")
  } else if (max_cycle_order > 7) {
    message("order_cycle >7 is not implemented. Use an integer between 3 and 7.")
    max_cycle_order <- 7
  }

  n <- nrow(A)

  if (is.na(n_rep)) {
    n_rep <- ceiling((1 / (2 * alpha) *
                        stats::qnorm(1 - alpha / (2 * (max_cycle_order - 1))))^2)
    message(paste("Use n_rep =", n_rep))
  }

  if (is.na(subsample_sizes)) {
    subsample_sizes <- auto_select_subsample_sizes(
      A, n_subsample_sizes, k_max = max_cycle_order,
      R = n_rep, alpha = alpha, delta = 0.05,
      max_subsample_size = max_subsample_size)
  }

  result <- .net_summary_subsample_adj(A, subsample_sizes, max_cycle_order, n_rep)
  colnames(result) <- c("v-shape", "triangle", "square", "pentagon",
                         "hexagon", "septagon")[1:(max_cycle_order - 1)]
  if (!is.na(y_max) & ((y_max > 1) | (y_max < 0))) {
    warning("y_max: Use a number between 0 and 1")
    y_max <- NA
  }
  if (is.na(y_max)) {
    y_max <- max(result)
  }

  result <- data.frame(result)
  melt_result <- suppressMessages(reshape2::melt(result,
                           value.name = "summaries", variable.name = "subgraphs"))
  xlabels <- get_netsummary_xlables(max_cycle_order)
  p <- ggplot2::ggplot(melt_result, ggplot2::aes(.data$subgraphs, .data$summaries))
  p <- p + ggplot2::geom_violin() + ggplot2::ylim(0, y_max)
  p <- p + ggplot2::stat_summary(fun = mean, geom = "point", size = 2)
  p <- p + ggplot2::ylab("Prevalence and local variability")
  p <- p + ggplot2::scale_x_discrete(name = NULL, labels = xlabels)
  p <- p + ggplot2::theme(
    axis.text.x = ggtext::element_markdown(color = "black", size = rel(1))
  )

  if (save_plot) {
    ggplot2::ggsave(filename, plot = p, width = width, height = height, unit = "in")
    return(invisible(p))
  }
  return(p)
}

##' @rdname netsummary_plot
##' @export
violin_netsummary <- function(A, ...) {
  .Deprecated("netsummary_plot")
  netsummary_plot(A, ...)
}

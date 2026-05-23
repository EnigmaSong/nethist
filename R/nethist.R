##' Network histogram estimation for single-layer networks
##'
##' Estimating a network histogram for a single-layer network and returning
##' the indices of partitions.
##'
##' @param A An adjacency matrix or graph object. Accepted formats: a
##'   \code{matrix}, sparse \code{dgCMatrix}, \code{igraph} object, or
##'   \code{network} object. Must be an undirected simple graph.
##' @param h A bandwidth parameter. If \code{NA}, the bandwidth is selected
##'   by Olhede and Wolfe (2014). If specified, the user-supplied value is used.
##' @param method Type of loss function. One of \code{"PLL"} (default, profile
##'   log-likelihood) or \code{"LSE"} (least squares). See Details.
##' @param control A control object from \code{\link{nethist_control}}.
##'   Governs \code{max_itr}, \code{greedy_swap_rule},
##'   \code{greedy_stop_threshold}, and \code{verbose}.
##' @param ... `r lifecycle::badge("deprecated")` Pass \code{max_itr},
##'   \code{greedy_swap_rule}, \code{greedy_stop_threshold}, or \code{verbose}
##'   via \code{control = nethist_control(...)} instead.
##' @returns An object of class \code{"nethist"} with the following fields:
##' \itemize{
##'   \item \code{cluster} an integer vector of length \eqn{n} with block
##'     assignments.
##'   \item \code{thetahat} a \eqn{k \times k} probability matrix ordered by
##'     block labels.
##'   \item \code{rho_hat} estimated sparsity parameter.
##'   \item \code{normalized_LL} normalized log-likelihood.
##'   \item \code{MSE} mean squared error.
##'   \item \code{method} loss function used (\code{"PLL"} or \code{"LSE"}).
##' }
##' @details
##' \code{method = "PLL"} is for Olhede and Wolfe (2014).
##' \code{method = "LSE"} is for Gao et al. (2015).
##'
##' Note that \code{cluster} labels are not ordered: nodes in cluster 1 are not
##' necessarily more similar to cluster 2 than to cluster 10. Users may specify
##' a custom display order in \code{\link{plot.nethist}}.
##' @examples
##' \donttest{
##' set.seed(42)
##' data(polblog)
##' fit <- nethist(polblog)
##' fit
##' plot(fit)
##'
##' fit_h <- nethist(polblog, h = 72)
##' }
##' @seealso \code{\link{multinethist}}, \code{\link{plot.nethist}},
##'   \code{\link{nethist_control}}
##' @references Olhede, S. C. & Wolfe, P. J. (2014). Network Histograms and
##'   Universality of Blockmodel Approximation. PNAS, 111(41), 14722-14727.
##'   doi:10.1073/pnas.1400374111
##' @references Gao, C., Lu, Y., & Zhou, H. H. (2015). Rate-Optimal Graphon
##'   Estimation. The Annals of Statistics, 43(6), 2624-2652.
##'   doi:10.1214/15-AOS1354
##' @export
nethist <- function(A, h = NA,
                    method = "PLL",
                    control = nethist_control(), ...) {
  UseMethod("nethist")
}

##' @exportS3Method
nethist.igraph <- function(A, h = NA,
                            method = "PLL",
                            control = nethist_control(), ...) {
  A <- igraph::as_adjacency_matrix(A, sparse = FALSE)
  nethist.default(A, h = h, method = method, control = control, ...)
}

##' @exportS3Method
nethist.matrix <- function(A, h = NA,
                            method = "PLL",
                            control = nethist_control(), ...) {
  nethist.default(A, h = h, method = method, control = control, ...)
}

##' @exportS3Method
nethist.dgCMatrix <- function(A, h = NA,
                               method = "PLL",
                               control = nethist_control(), ...) {
  A <- as.matrix(A)
  nethist.default(A, h = h, method = method, control = control, ...)
}

##' @exportS3Method
nethist.combined_networks <- function(A, h = NA,
                                      method = "PLL",
                                      control = nethist_control(), ...) {
  stop("A is a combined_networks (multilayer) object. Use multinethist() instead.")
}

##' @exportS3Method
nethist.network <- function(A, h = NA,
                             method = "PLL",
                             control = nethist_control(), ...) {
  A <- as.matrix(A, matrix.type = "adjacency")
  nethist.default(A, h = h, method = method, control = control, ...)
}

##' @exportS3Method
nethist.default <- function(A, h = NA,
                             method = "PLL",
                             control = nethist_control(), ...) {
  if (!is.matrix(A)) stop("A must be a matrix.")
  if (nrow(A) != ncol(A)) stop("A must be a square matrix.")
  if (any(diag(A) != 0))
    stop("A has self-loops. Network A must be a simple graph (no self-loops).")
  if (!isSymmetric(unname(A)))
    stop("A is not symmetric. Network A must be an undirected graph.")
  if (any(A[A != 0] != 1))
    stop("A is not a simple graph. All non-zero entries must be 1 (binary adjacency matrix).")
  control <- .resolve_control(control, list(...))
  result <- multinethist.array(array(A, dim = c(nrow(A), ncol(A), 1)),
                               h, common_f = FALSE, method, control = control)
  result$homogeneous <- NULL
  return(result)
}

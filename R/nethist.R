##' @rdname multinethist
##' @param ... `r lifecycle::badge("deprecated")` Pass \code{max_itr},
##'   \code{greedy_swap_rule}, \code{greedy_stop_threshold}, or \code{verbose}
##'   via \code{control = nethist_control(...)} instead.
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
  multinethist.array(array(A, dim = c(nrow(A), ncol(A), 1)),
                     h, common_f = FALSE, method, control = control)
}

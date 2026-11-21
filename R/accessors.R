# -- Accessor generics ---------------------------------------------------------

##' Extract components of a network histogram fit
##'
##' Accessor functions for fitted network histograms. They return the
##' quantities stored in a fit, so that code does not rely on the internal
##' structure of the object, which may change in future versions.
##'
##' \code{block_id()} returns the block assignment of each vertex.
##' \code{bandwidth()} returns the bandwidth, that is, the group size used to
##' construct the network histogram.
##' \code{shape_id()} returns, for a hybrid network histogram, the shape
##' assigned to each block as a symmetric \eqn{k \times k} matrix, where
##' \eqn{k} is the number of blocks. Shapes are numbered in increasing order
##' of their estimated connection probability, so shape 1 has the lowest
##' probability.
##'
##' The estimated block probabilities are returned by
##' \code{\link[=coef.nethist]{coef()}}.
##'
##' @param object a fit returned by \code{\link{nethist}},
##'   \code{\link{multinethist}}, or \code{\link{hnethist}}. For
##'   \code{shape_id()}, a fit returned by \code{\link{hnethist}}.
##' @param ... currently unused.
##' @returns
##' \code{block_id()} returns an integer vector of length \eqn{n} with the block
##' label of each vertex.
##' \code{bandwidth()} returns a single number.
##' \code{shape_id()} returns an integer \eqn{k \times k} matrix whose entry
##' \eqn{(a, b)} is the shape of the block formed by groups \eqn{a} and
##' \eqn{b}.
##' @examples
##' \donttest{
##' set.seed(42)
##' data(polblog)
##' fit <- hnethist(polblog)
##' bandwidth(fit)
##' table(block_id(fit))
##' shape_id(fit)
##' }
##' @seealso \code{\link{coef.nethist}}, \link{network_histogram-class}
##' @name network_histogram-accessors
NULL

##' @rdname network_histogram-accessors
##' @export
block_id <- function(object, ...) {
  UseMethod("block_id")
}

##' @rdname network_histogram-accessors
##' @exportS3Method
block_id.network_histogram <- function(object, ...) {
  return(as.integer(object$cluster))
}

##' @rdname network_histogram-accessors
##' @export
bandwidth <- function(object, ...) {
  UseMethod("bandwidth")
}

##' @rdname network_histogram-accessors
##' @exportS3Method
bandwidth.network_histogram <- function(object, ...) {
  return(object$h)
}

##' @rdname network_histogram-accessors
##' @export
shape_id <- function(object, ...) {
  UseMethod("shape_id")
}

##' @rdname network_histogram-accessors
##' @exportS3Method
shape_id.default <- function(object, ...) {
  stop("shape_id() is only defined for hnethist fits.", call. = FALSE)
}

##' @rdname network_histogram-accessors
##' @exportS3Method
shape_id.hnethist <- function(object, ...) {
  k       <- length(unique(object$cluster))
  centers <- object$blockcluster$centers[, 1]
  # kmeans labels are arbitrary, so rank shapes by their probability.
  rank_of <- match(seq_along(centers), order(centers))
  ids     <- rank_of[object$blockcluster$cluster]

  mat <- matrix(0L, k, k)
  mat[upper.tri(mat, diag = TRUE)] <- ids
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
  return(mat)
}

# -- coef methods --------------------------------------------------------------

##' Estimated block probabilities of a network histogram fit
##'
##' Returns the estimated block probability matrix of a fitted network
##' histogram. Unlike \code{\link[=fitted.nethist]{fitted()}}, which expands
##' the estimate to pairs of vertices, \code{coef()} returns one value per
##' pair of blocks.
##'
##' For an \code{hnethist} fit, the probabilities are those of the selected
##' hybrid model, in which blocks of the same shape share one value.
##'
##' @param object a fit returned by \code{\link{nethist}},
##'   \code{\link{hnethist}}, or \code{\link{multinethist}}.
##' @param type one of \code{"nethist"} (default) or \code{"prob"}.
##'   \code{"nethist"} returns the network histogram, the block probabilities
##'   divided by the edge density \code{rho_hat}. \code{"prob"} returns the
##'   block probabilities.
##' @param ... currently unused.
##' @returns
##' For \code{nethist} and \code{hnethist} fits, a \eqn{k \times k} matrix,
##' where \eqn{k} is the number of blocks. For \code{multinethist} fits, a
##' \eqn{k \times k \times L} array, where \eqn{L} is the number of layers.
##' @examples
##' \donttest{
##' set.seed(42)
##' data(polblog)
##' fit <- nethist(polblog)
##' coef(fit, type = "prob")
##' }
##' @seealso \code{\link{fitted.nethist}}, \link{network_histogram-accessors}
##' @exportS3Method
coef.nethist <- function(object, type = "nethist", ...) {
  type <- match.arg(type, c("nethist", "prob"))
  # A fit that collapses to a single block stores thetahat as a scalar.
  th <- as.matrix(object$thetahat)
  if (type == "nethist") {
    th <- th / object$rho_hat
  }
  return(th)
}

##' @rdname coef.nethist
##' @exportS3Method
coef.multinethist <- function(object, type = "nethist", ...) {
  type <- match.arg(type, c("nethist", "prob"))
  th <- object$thetahat
  # A fit that collapses to a single block stores one value per layer.
  if (is.null(dim(th))) {
    th <- array(th, dim = c(1L, 1L, length(th)))
  }
  if (type == "nethist") {
    for (l in seq_len(dim(th)[3L])) {
      th[, , l] <- th[, , l] / object$rho_hat[l]
    }
  }
  return(th)
}

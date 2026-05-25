# -- fitted methods for nethist objects ---------------------------------------

##' Fitted values for a nethist model
##'
##' Returns a submatrix of the fitted block-probability or graphon matrix
##' for a `nethist` or `hnethist` object.
##'
##' @param object a `nethist` or `hnethist` object from [nethist()] or
##'   [hnethist()].
##' @param set1 integer vector of vertex indices for rows. `NULL` uses all vertices.
##' @param set2 integer vector of vertex indices for columns. `NULL` uses all
##'   vertices.
##' @param type one of `"nethist"` (default) or `"prob"`. `"nethist"` returns
##'   the graphon estimate `thetahat / rho_hat`; `"prob"` returns the raw
##'   block probability matrix `thetahat`.
##' @param ... currently unused.
##' @returns A numeric matrix of dimension `|set1| x |set2|`.
##' @seealso [nethist()], [hnethist()], [multinethist()]
##' @examples
##' \donttest{
##' set.seed(1)
##' A <- igraph::as_adjacency_matrix(
##'   igraph::sample_gnp(60, 0.3), sparse = FALSE)
##' fit <- nethist(A, h = 10L)
##' fitted(fit, set1 = 1:10, set2 = 1:10)
##' fitted(fit, type = "prob")
##' }
##' @exportS3Method
fitted.nethist <- function(object, set1 = NULL, set2 = NULL,
                           type = "nethist", ...) {
  type <- match.arg(type, c("nethist", "prob"))
  n    <- length(object$cluster)
  s1   <- if (is.null(set1)) seq_len(n) else set1
  s2   <- if (is.null(set2)) seq_len(n) else set2

  if (is.null(set1) && is.null(set2) && n > 3000L)
    warning("Returning full ", n, " x ", n, " matrix. ",
            "Consider using set1/set2 to limit scope.",
            call. = FALSE)

  th <- switch(type,
    nethist = object$thetahat / object$rho_hat,
    prob    = object$thetahat
  )
  return(th[object$cluster[s1], object$cluster[s2], drop = FALSE])
}

##' Fitted values for a multinethist model
##'
##' Returns a subarray of the fitted block-probability or graphon matrix
##' for a `multinethist` object.
##'
##' @param object a `multinethist` object from [multinethist()].
##' @param set1 integer vector of vertex indices for rows. `NULL` uses all vertices.
##' @param set2 integer vector of vertex indices for columns. `NULL` uses all
##'   vertices.
##' @param layer integer vector of layer indices. `NULL` uses all layers.
##' @param type one of `"nethist"` (default) or `"prob"`. `"nethist"` returns
##'   `thetahat[,,l] / rho_hat[l]` per layer; `"prob"` returns raw `thetahat`.
##' @param drop logical. If `TRUE` (default) and a single layer is selected,
##'   the layer dimension is dropped and a matrix is returned. If `FALSE`, a
##'   3-dimensional array is always returned.
##' @param ... currently unused.
##' @returns A numeric matrix of dimension `|set1| x |set2|` when a single
##'   layer is selected and `drop = TRUE`, otherwise a 3-dimensional array of
##'   dimension `|set1| x |set2| x |layer|`.
##' @seealso [multinethist()], [nethist()]
##' @examples
##' \donttest{
##' data(IndianVil)
##' fit <- multinethist(IndianVil, h = 20L)
##' fitted(fit, set1 = 1:10, set2 = 1:10, layer = 1)
##' fitted(fit, layer = c(1, 2), drop = FALSE)
##' }
##' @exportS3Method
fitted.multinethist <- function(object, set1 = NULL, set2 = NULL,
                                layer = NULL, type = "nethist",
                                drop = TRUE, ...) {
  type <- match.arg(type, c("nethist", "prob"))
  n    <- length(object$cluster)
  L    <- dim(object$thetahat)[3L]
  s1   <- if (is.null(set1)) seq_len(n) else set1
  s2   <- if (is.null(set2)) seq_len(n) else set2
  ll   <- if (is.null(layer)) seq_len(L) else layer

  if (is.null(set1) && is.null(set2) && n > 3000L)
    warning("Returning full ", n, " x ", n, " x ", length(ll),
            " array. Consider using set1/set2 to limit scope.",
            call. = FALSE)

  th <- object$thetahat[, , ll, drop = FALSE]
  if (type == "nethist") {
    rho <- object$rho_hat[ll]
    for (i in seq_along(ll))
      th[, , i] <- th[, , i] / rho[i]
  }

  result <- th[object$cluster[s1], object$cluster[s2], ,
               drop = FALSE]
  if (drop && length(ll) == 1L)
    return(result[, , 1L])
  return(result)
}

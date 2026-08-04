##' 3D histogram plot for nethist objects
##'
##' Drawing [plot3D::hist3D()] using a `nethist`, `multinethist`, or `hnethist`
##' object with a user-specified order.
##'
##' @param x a `nethist`, `multinethist`, or `hnethist` object.
##' @param idx_order A numeric vector for index label order, which must be a
##'   permutation of `x$cluster`. If `NA`, it uses `1:max(x$cluster)`.
##' @param type One of `"nethist"` (default) or `"prob"`. `"MNhist"` is a
##'   deprecated alias for `"nethist"`.
##' @param ... Other arguments passed to [plot3D::hist3D()].
##' @returns Called for its side effects (plotting). Returns `NULL` invisibly.
##' @examples
##' \donttest{
##' set.seed(42)
##' data(IndianVil)
##' mnhist_Ind_vil <- multinethist(IndianVil)
##' plot3d(mnhist_Ind_vil)
##'
##' data(polblog)
##' fit <- nethist(polblog)
##' plot3d(fit)
##' }
##' @importFrom plot3D hist3D
##' @export
plot3d <- function(x, idx_order = 1:max(x$cluster),
                   type = "nethist",
                   ...) {
  UseMethod("plot3d")
}

##' @rdname plot3d
##' @exportS3Method
plot3d.nethist <- function(x,
                           idx_order = 1:max(x$cluster),
                           type = "nethist",
                           ...) {
  k <- max(x$cluster)
  if (!.is_valid_order(idx_order, 1:k)) {
    warning(paste0("idx_order is invalid. Set idx_order = 1:", k))
    idx_order <- 1:k
  }
  mat <- switch(type,
    nethist = x$thetahat[idx_order, idx_order] / x$rho_hat,
    prob    = x$thetahat[idx_order, idx_order],
    stop("type must be one of nethist or prob."))
  hist3D(z = mat, ...)
  return(invisible(NULL))
}

##' @rdname plot3d
##' @exportS3Method
plot3d.hnethist <- function(x,
                            idx_order = 1:max(x$cluster),
                            type = "nethist",
                            ...) {
  plot3d.nethist(x, idx_order = idx_order, type = type, ...)
}

##' @rdname plot3d
##' @exportS3Method
plot3d.multinethist <- function(x,
                                idx_order = 1:max(x$cluster),
                                type = "nethist",
                                ...) {
  if (type == "MNhist") {
    warning("type = \"MNhist\" is deprecated. Use type = \"nethist\" instead.",
            call. = FALSE)
    type <- "nethist"
  }
  k <- max(x$cluster)
  if (!.is_valid_order(idx_order, 1:k)) {
    warning(paste0("idx_order is invalid. Set idx_order = 1:", k))
    idx_order <- 1:k
  }
  n_loop <- ifelse(x$homogeneous & (type == "nethist"), 1,
                   dim(x$thetahat)[3])
  for (l in seq_len(n_loop)) {
    mat <- switch(type,
      nethist = x$thetahat[idx_order, idx_order, l] / x$rho_hat[l],
      prob    = x$thetahat[idx_order, idx_order, l],
      stop("type must be one of nethist or prob."))
    hist3D(z = mat, main = paste("Layer", l), ...)
  }
  return(invisible(NULL))
}

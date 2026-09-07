# -- Class hierarchy ----------------------------------------------------------

##' Network histogram class hierarchy
##'
##' All fitted objects returned by \code{\link{nethist}},
##' \code{\link{multinethist}}, and \code{\link{hnethist}} inherit from the
##' common base class \code{network_histogram}:
##'
##' \itemize{
##' \item \code{nethist} fits have class
##'   \code{c("nethist", "network_histogram")}.
##' \item \code{multinethist} fits have class
##'   \code{c("multinethist", "network_histogram")}.
##' \item \code{hnethist} fits have class
##'   \code{c("hnethist", "nethist", "network_histogram")}.
##' }
##'
##' Every \code{network_histogram} object carries \code{cluster},
##' \code{thetahat}, \code{rho_hat}, \code{normalized_LL}, \code{method}, and
##' \code{h}. Methods that only depend on those fields are defined once on the
##' base class, so \code{\link{logLik.network_histogram}} and
##' \code{\link{summary.network_histogram}} are shared by all three. Methods
##' whose behaviour depends on the layer dimension, such as \code{plot()},
##' \code{plot3d()}, and \code{fitted()}, remain specific to the subclasses.
##'
##' @name network_histogram-class
##' @seealso \code{\link{logLik.network_histogram}},
##'   \code{\link{summary.network_histogram}}, \code{\link{print.nethist}}
NULL

# multinethist() serves both the single-layer and the multilayer case.
##' @noRd
.network_histogram_class <- function(n_layers) {
  leaf <- if (n_layers > 1) "multinethist" else "nethist"
  return(c(leaf, "network_histogram"))
}

# -- Class-specific model dimensions ------------------------------------------

# Free block probabilities of the fitted model.
##' @noRd
.model_df <- function(x) UseMethod(".model_df")

##' @noRd
.model_df.nethist <- function(x) {
  k <- length(unique(x$cluster))
  return(k * (k + 1) / 2)
}

##' @noRd
.model_df.multinethist <- function(x) {
  k <- length(unique(x$cluster))
  # A common graphon across layers shares one set of block probabilities.
  if (isTRUE(x$homogeneous)) return(k * (k + 1) / 2)
  return(length(x$rho_hat) * k * (k + 1) / 2)
}

##' @noRd
.model_df.hnethist <- function(x) {
  return(x$s)
}

# Number of Bernoulli edge variables.
##' @noRd
.model_nobs <- function(x) {
  n <- length(x$cluster)
  return(length(x$rho_hat) * n * (n - 1) / 2)
}

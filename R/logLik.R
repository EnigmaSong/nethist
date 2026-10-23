# The fitted object stores only the per-edge normalized log-likelihood, so the
# edge count is recovered from rho_hat rather than kept as an extra field, which
# keeps objects fitted by earlier versions usable.
##' @noRd
.nethist_nedges <- function(x) {
  n <- length(x$cluster)
  return(round(sum(x$rho_hat) * n * (n - 1) / 2))
}

##' Log-likelihood of network histogram fits
##'
##' Extracts the Bernoulli log-likelihood of a fitted network histogram. The
##' method is defined on the common base class, so it applies to
##' \code{nethist}, \code{multinethist}, and \code{hnethist} objects alike; see
##' \link{network_histogram-class}.
##'
##' @param object a \code{network_histogram} object, that is, a fit returned by
##'   \code{\link{nethist}}, \code{\link{multinethist}}, or
##'   \code{\link{hnethist}}.
##' @param normalized logical. If \code{FALSE} (default), the unnormalized
##'   log-likelihood is returned as a \code{logLik} object. If \code{TRUE}, the
##'   per-edge normalized value stored in \code{object$normalized_LL} is
##'   returned as a plain numeric.
##' @param ... currently unused.
##' @returns
##' With \code{normalized = FALSE}, an object of class \code{logLik} carrying
##' \code{df} and \code{nobs} attributes, so that [stats::AIC()] and
##' [stats::BIC()] work directly on the fit. With \code{normalized = TRUE}, an
##' unclassed numeric equal to \code{object$normalized_LL}.
##' @details
##' The normalized value is not given class \code{logLik} on purpose. Passing a
##' per-edge normalized value to \code{AIC()} or \code{BIC()} would return a
##' meaningless number without any error.
##'
##' The degrees of freedom count the free block probabilities: \code{k * (k + 1)
##' / 2} for \code{nethist}, where \code{k} is the number of blocks; the same
##' for \code{multinethist} fitted with \code{common_f = TRUE}, and \code{L}
##' times that for \code{L} layers otherwise; and \code{s}, the number of
##' selected shapes, for \code{hnethist}.
##'
##' The number of observations is the number of Bernoulli edge variables, i.e.
##' \code{n * (n - 1) / 2} for a single layer and \code{L} times that for
##' \code{L} layers.
##'
##' The log-likelihood is evaluated at the fitted probability matrix regardless
##' of \code{method}, so it is also available for \code{method = "LSE"} fits.
##' @examples
##' \donttest{
##' set.seed(42)
##' data(polblog)
##' fit <- suppressMessages(nethist(polblog))
##' logLik(fit)
##' logLik(fit, normalized = TRUE)
##' AIC(fit)
##' BIC(fit)
##' }
##' @seealso \link{network_histogram-class},
##'   \code{\link{summary.network_histogram}}
##' @rdname logLik.network_histogram
##' @exportS3Method
logLik.network_histogram <- function(object, normalized = FALSE, ...) {
  if (normalized) return(object$normalized_LL)
  val <- object$normalized_LL * .nethist_nedges(object)
  return(structure(val, class = "logLik",
                   df = .model_df(object), nobs = .model_nobs(object)))
}

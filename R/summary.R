# -- Shared display helpers ---------------------------------------------------

##' @noRd
.nethist_method_label <- function(x) {
  return(ifelse(x$method == "PLL", "Profile Likelihood",
                "Least Squared Error"))
}

# thetahat is a matrix for nethist, a 3d array for multinethist, and a bare
# scalar or vector when the fit collapses to a single block.
##' @noRd
.nethist_theta_desc <- function(x) {
  d <- dim(x$thetahat)
  if (is.null(d)) {
    return(paste(length(x$thetahat),
                 ngettext(length(x$thetahat), "probability value",
                          "probability values")))
  }
  if (length(d) == 2L) {
    return(paste0(d[1L], " x ", d[2L], " probability matrix"))
  }
  return(paste0(paste(d, collapse = " x "), " probability array"))
}

# multinethist() does not compute an LSE, so MSE can be absent.
##' @noRd
.nethist_loss_line <- function(x, digits) {
  if (x$method == "PLL" || is.null(x$MSE)) {
    return(paste0("Normalized log-likelihood: ",
                  format(x$normalized_LL, digits = digits)))
  }
  return(paste0("Mean squared error: ", format(x$MSE, digits = digits)))
}

# Nodes / Blocks / Bandwidth, plus Layers for multilayer fits.
##' @noRd
.nethist_dims_line <- function(x) {
  fields <- c(paste("Nodes:", length(x$cluster)),
              paste("Blocks:", length(unique(x$cluster))),
              paste("Bandwidth:", x$h))
  if (inherits(x, "multinethist")) {
    fields <- c(fields, paste("Layers:", length(x$rho_hat)))
  }
  return(paste(fields, collapse = "   "))
}

##' @noRd
.nethist_title <- function(x) {
  if (inherits(x, "hnethist")) {
    return("Hybrid network histogram fit (hnethist)")
  }
  if (inherits(x, "multinethist")) {
    return("Multilayer network histogram fit (multinethist)")
  }
  return("Network histogram fit (nethist)")
}

# -- summary methods ----------------------------------------------------------

##' Summarize network histogram fits
##'
##' Detailed summary of a fitted network histogram. The \code{print} methods of
##' the fitted objects give a short identification only; the model fit
##' statistics are reported here.
##'
##' The \code{summary} method is defined on the common base class, so it
##' applies to \code{nethist}, \code{multinethist}, and \code{hnethist} objects
##' alike. The subclass methods add the fields specific to them and are
##' otherwise deferred to the base method via \code{NextMethod()}; see
##' \link{network_histogram-class}.
##'
##' @param object a \code{network_histogram} object, that is, a fit returned by
##'   \code{\link{nethist}}, \code{\link{multinethist}}, or
##'   \code{\link{hnethist}}.
##' @param x a \code{summary.network_histogram} object.
##' @param digits number of significant digits used when printing numeric
##'   quantities.
##' @param ... additional arguments. Currently unused.
##' @returns
##' \code{summary()} returns an object inheriting from
##' \code{summary.network_histogram}: a list holding the fitted object together
##' with the number of nodes, the number of blocks, the log-likelihood, the
##' degrees of freedom, \code{AIC}, and \code{BIC}. Multilayer fits add
##' \code{n_layers}, and \code{hnethist} fits add the selected number of shapes
##' \code{s} and the number of candidate models \code{n_candidate}. The
##' \code{print} method returns its input invisibly.
##' @details
##' The estimated probabilities themselves are not printed, since the matrix
##' grows with the number of blocks. Use \code{object$thetahat} or
##' \code{\link[=fitted.nethist]{fitted()}} to obtain them.
##' @examples
##' \donttest{
##' set.seed(42)
##' data(polblog)
##' fit <- suppressMessages(nethist(polblog))
##' fit
##' summary(fit)
##' }
##' @seealso \link{network_histogram-class}, \code{\link{print.nethist}},
##'   \code{\link{logLik.network_histogram}}
##' @rdname summary.network_histogram
##' @exportS3Method
summary.network_histogram <- function(object, ...) {
  ll <- stats::logLik(object)
  return(structure(
    list(object    = object,
         n_nodes   = length(object$cluster),
         n_blocks  = length(unique(object$cluster)),
         logLik    = as.numeric(ll),
         df        = attr(ll, "df"),
         AIC       = stats::AIC(object),
         BIC       = stats::BIC(object)),
    class = c("summary.nethist", "summary.network_histogram")
  ))
}

##' @rdname summary.network_histogram
##' @exportS3Method
summary.multinethist <- function(object, ...) {
  out <- NextMethod()
  out$n_layers <- length(object$rho_hat)
  return(structure(out, class = c("summary.multinethist",
                                  "summary.network_histogram")))
}

##' @rdname summary.network_histogram
##' @exportS3Method
summary.hnethist <- function(object, ...) {
  out <- NextMethod()
  out$s           <- object$s
  out$n_candidate <- length(object$details)
  return(structure(out, class = c("summary.hnethist",
                                  "summary.network_histogram")))
}

# -- print method for summary objects -----------------------------------------

##' @rdname summary.network_histogram
##' @exportS3Method
print.summary.network_histogram <- function(x, digits = getOption("digits"),
                                            ...) {
  object <- x$object
  cat("\n", .nethist_title(object), "\n\n", sep = "")
  cat(.nethist_dims_line(object), "\n", sep = "")
  if (inherits(object, "multinethist")) {
    cat("Common graphon across layers: ", isTRUE(object$homogeneous),
        "\n", sep = "")
  }
  cat("Method: ", .nethist_method_label(object), "\n", sep = "")

  if (length(object$rho_hat) > 1L) {
    cat("Edge density (rho_hat) by layer:\n")
    print(format(object$rho_hat, digits = digits), quote = FALSE)
    cat("\n")
  } else {
    cat("Edge density (rho_hat): ", format(object$rho_hat, digits = digits),
        "\n\n", sep = "")
  }

  if (!is.null(x$s)) {
    cat("Shape selection:\n")
    cat("  Selected shapes: ", x$s, " of ", x$n_candidate,
        " candidates\n", sep = "")
    cat("  BIC of selected model: ", format(object$BIC, digits = digits),
        "\n\n", sep = "")
  }

  cat("Model fit:\n")
  cat("  Normalized log-likelihood: ",
      format(object$normalized_LL, digits = digits), "\n", sep = "")
  cat("  Log-likelihood: ", format(x$logLik, digits = digits),
      "  (df = ", x$df, ")\n", sep = "")
  cat("  AIC: ", format(x$AIC, digits = digits),
      "   BIC: ", format(x$BIC, digits = digits), "\n", sep = "")
  if (!is.null(object$MSE)) {
    cat("  Mean squared error: ", format(object$MSE, digits = digits),
        "\n", sep = "")
  }

  cat("\nthetahat: ", .nethist_theta_desc(object),
      " (use fitted() or $thetahat)\n", sep = "")
  return(invisible(x))
}

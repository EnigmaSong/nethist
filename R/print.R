##' Print nethist objects
##'
##' Compact one-screen description of \code{nethist}, \code{multinethist}, and
##' \code{hnethist} objects. The estimated probabilities are reported by their
##' dimensions only; use \code{\link{summary.network_histogram}} for the model
##' fit statistics and \code{object$thetahat} for the values themselves.
##'
##' @param x a \code{nethist}, \code{multinethist}, or \code{hnethist} object.
##' @param digits number of significant digits used when printing the loss.
##' @param ... additional arguments. Currently unused.
##' @returns The input object, invisibly.
##' @examples
##' \donttest{
##' set.seed(42)
##' data(polblog)
##' fit <- suppressMessages(nethist(polblog))
##' print(fit)
##' summary(fit)
##' }
##' @seealso \code{\link{summary.network_histogram}}
##' @exportS3Method
print.nethist <- function(x, digits = getOption("digits"), ...){
  cat("\nNetwork histogram fit (nethist)\n\n")
  cat(.nethist_dims_line(x), "\n", sep = "")
  cat("Method: ", .nethist_method_label(x), "\n", sep = "")
  cat("thetahat: ", .nethist_theta_desc(x), "\n", sep = "")
  cat(.nethist_loss_line(x, digits), "\n", sep = "")
  cat("\nUse summary() for details.\n")
  return(invisible(x))
}

##' @rdname print.nethist
##' @exportS3Method
print.multinethist <- function(x, digits = getOption("digits"), ...){
  cat("\nMultilayer network histogram fit (multinethist)\n\n")
  cat(.nethist_dims_line(x), "\n", sep = "")
  cat("Common graphon across layers: ", isTRUE(x$homogeneous), "\n", sep = "")
  cat("Method: ", .nethist_method_label(x), "\n", sep = "")
  cat("thetahat: ", .nethist_theta_desc(x), "\n", sep = "")
  cat(.nethist_loss_line(x, digits), "\n", sep = "")
  cat("\nUse summary() for details.\n")
  return(invisible(x))
}

##' @rdname print.nethist
##' @exportS3Method
print.hnethist <- function(x, digits = getOption("digits"), ...){
  cat("\nHybrid network histogram fit (hnethist)\n\n")
  cat(.nethist_dims_line(x), "\n", sep = "")
  cat("Method: ", .nethist_method_label(x), "\n", sep = "")
  cat("Selected shapes: ", x$s, " (BIC: ", format(x$BIC, digits = digits),
      ")\n", sep = "")
  cat("thetahat: ", .nethist_theta_desc(x), "\n", sep = "")
  cat(.nethist_loss_line(x, digits), "\n", sep = "")
  cat("\nUse summary() for details.\n")
  return(invisible(x))
}

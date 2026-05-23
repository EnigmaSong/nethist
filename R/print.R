##' @noRd
.print_nethist_footer <- function(x) {
  cat(paste0("\nMethod: ", ifelse(x$method == "PLL",
                                  "Profile Likelihood",
                                  "Least Squared Error"),
             "\n"))
  cat(paste0("\n", ifelse(x$method == "PLL",
                          paste("normalized likelihood:", x$normalized_LL, sep = "\n"),
                          paste("Mean squared error:", x$MSE, sep = "\n")), "\n"))
  cat("\nAvailable components:\n", sep = "\n")
  print(names(x))
}

##' Print nethist objects
##'
##' Prints the estimated probability matrix and model summary for
##' \code{nethist}, \code{multinethist}, and \code{hnethist} objects.
##'
##' @param x a \code{nethist}, \code{multinethist}, or \code{hnethist} object.
##' @param ... additional arguments passed to [print()].
##' @returns The input object, invisibly.
##' @examples
##' \donttest{
##' set.seed(42)
##' data(polblog)
##' fit <- suppressMessages(nethist(polblog))
##' print(fit)
##' }
##' @exportS3Method
print.nethist <- function(x, ...){
  cat("\nthetahat:\n")
  print(x$thetahat, ...)
  .print_nethist_footer(x)
  return(invisible(x))
}

##' @rdname print.nethist
##' @exportS3Method
print.multinethist <- function(x, ...){
  cat("\nTheta_hat:\n")
  print(x$thetahat[,,], ...)
  .print_nethist_footer(x)
  return(invisible(x))
}

##' @rdname print.nethist
##' @exportS3Method
print.hnethist <- function(x, ...){
  cat("\nTheta_hat:\n")
  print(x$thetahat, ...)
  cat(paste0("\nSelected shapes: ", x$s,
             " (BIC: ", round(x$BIC, 4), ")\n"))
  .print_nethist_footer(x)
  return(invisible(x))
}

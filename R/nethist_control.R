##' Control parameters for network histogram algorithms
##'
##' Constructs a control object for \code{\link{nethist}},
##' \code{\link{multinethist}}, and \code{\link{hnethist}}.
##'
##' @param algorithm character. Optimization algorithm. Currently only
##'   \code{"greedy"} is implemented.
##' @param max_itr integer. Maximum number of iterations. Default is
##'   \eqn{5 \times 10^6}.
##' @param greedy_swap_rule character. Node-pair selection rule for the greedy
##'   search. At each iteration, two nodes are drawn and their block labels
##'   are swapped if the move improves the objective. Currently only
##'   \code{"single_random"} (one pair drawn uniformly at random) is
##'   implemented.
##' @param greedy_stop_threshold integer. Early stopping criterion for the
##'   greedy search; the algorithm terminates if the objective has not improved
##'   for this many consecutive iterations. Default is 20,000.
##' @param verbose logical. Print progress messages during fitting.
##' @param ... Accepts deprecated argument names \code{swap_rule} and
##'   \code{consecutive_iter_threshold} with a warning.
##' @returns An object of class \code{"nethist_control"}.
##' @examples
##' # default control object
##' ctrl <- nethist_control()
##' print(ctrl)
##'
##' # reduce iteration limit for quick testing
##' ctrl <- nethist_control(max_itr = 1e4, greedy_stop_threshold = 100)
##'
##' \donttest{
##' data(polblog)
##' fit <- nethist(polblog, control = nethist_control(max_itr = 1e4))
##' }
##' @seealso \code{\link{nethist}}, \code{\link{multinethist}},
##'   \code{\link{hnethist}}
##' @export
nethist_control <- function(
  algorithm             = "greedy",
  max_itr               = 5e6,
  greedy_swap_rule      = "single_random",
  greedy_stop_threshold = 2e4,
  verbose               = FALSE,
  ...
) {
  # backward compat: old argument names
  dots <- list(...)
  old_to_new <- c(swap_rule                  = "greedy_swap_rule",
                  consecutive_iter_threshold = "greedy_stop_threshold")
  hits <- intersect(names(dots), names(old_to_new))
  if (length(hits)) {
    new_hits <- old_to_new[hits]
    warning(paste0(
      "Argument(s) ", paste(hits, collapse = ", "),
      " have been renamed to ", paste(new_hits, collapse = ", "),
      ". Please update your code."
    ), call. = FALSE)
    for (nm in hits) assign(old_to_new[nm], dots[[nm]])
  }

  if (!is.character(algorithm) || length(algorithm) != 1L)
    stop("algorithm must be a character string.")
  if (is.na(pmatch(algorithm, c("greedy"))))
    stop("algorithm must be one of: greedy.")
  if (!is.numeric(max_itr) || length(max_itr) != 1L || max_itr < 1L)
    stop("max_itr must be a positive number.")
  if (!is.character(greedy_swap_rule) || length(greedy_swap_rule) != 1L)
    stop("greedy_swap_rule must be a character string.")
  if (!is.numeric(greedy_stop_threshold) ||
      length(greedy_stop_threshold) != 1L ||
      greedy_stop_threshold < 1L)
    stop("greedy_stop_threshold must be a positive number.")
  if (!is.logical(verbose) || length(verbose) != 1L)
    stop("verbose must be TRUE or FALSE.")
  structure(
    list(
      algorithm             = algorithm,
      max_itr               = max_itr,
      greedy_swap_rule      = greedy_swap_rule,
      greedy_stop_threshold = greedy_stop_threshold,
      verbose               = verbose
    ),
    class = "nethist_control"
  )
}

##' @rdname nethist_control
##' @param x a \code{nethist_control} object.
##' @exportS3Method
print.nethist_control <- function(x, ...) {
  cat("nethist control parameters:\n")
  cat("  algorithm            :", x$algorithm, "\n")
  cat("  max_itr              :", x$max_itr, "\n")
  cat("  greedy_swap_rule     :", x$greedy_swap_rule, "\n")
  cat("  greedy_stop_threshold:", x$greedy_stop_threshold, "\n")
  cat("  verbose              :", x$verbose, "\n")
  invisible(x)
}

# Merge deprecated direct arguments into control.
# Warns once per call if any deprecated args are found in dots.
.resolve_control <- function(control, dots) {
  # maps old direct-arg names to current control field names
  deprecated_map <- c(max_itr                    = "max_itr",
                      swap_rule                  = "greedy_swap_rule",
                      consecutive_iter_threshold = "greedy_stop_threshold",
                      verbose                    = "verbose")
  hits <- intersect(names(dots), names(deprecated_map))
  if (length(hits)) {
    warning(paste0(
      "Passing ", paste(hits, collapse = ", "),
      " directly is deprecated. Use nethist_control() instead."
    ), call. = FALSE)
    for (nm in hits) control[[deprecated_map[nm]]] <- dots[[nm]]
  }
  return(control)
}

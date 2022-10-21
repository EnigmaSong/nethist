#' @noRd

# Falling factorial
.ffct <- function(x, k) {
  sapply(x, function(y,k) {
    prod(y:(y - k + 1))
  }, k= k)
}

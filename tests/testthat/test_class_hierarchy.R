# Tests for the network_histogram class hierarchy and the methods defined on it.

set.seed(42)
G      <- igraph::sample_gnp(60, p = 0.2)
A      <- igraph::as_adjacency_matrix(G, sparse = FALSE)
fit_nh <- suppressMessages(nethist(A, h = 12))
fit_hn <- suppressMessages(hnethist(A, h = 12))

set.seed(42)
A_multi <- array(0, dim = c(60L, 60L, 2L))
for (l in 1:2) {
  A_l <- igraph::as_adjacency_matrix(
    igraph::sample_gnp(60L, 0.2, directed = FALSE),
    sparse = FALSE
  )
  diag(A_l) <- 0
  A_multi[, , l] <- A_l
}
fit_mn <- suppressMessages(multinethist(A_multi, h = 12))

## Class vectors ----
test_that("every fit inherits from network_histogram", {
  expect_true(inherits(fit_nh, "network_histogram"))
  expect_true(inherits(fit_mn, "network_histogram"))
  expect_true(inherits(fit_hn, "network_histogram"))
})

test_that("the base class is the last element of every class vector", {
  for (fit in list(fit_nh, fit_mn, fit_hn)) {
    expect_identical(utils::tail(class(fit), 1L), "network_histogram")
  }
})

test_that("class vectors have the documented shape", {
  expect_equal(class(fit_nh), c("nethist", "network_histogram"))
  expect_equal(class(fit_mn), c("multinethist", "network_histogram"))
  expect_equal(class(fit_hn), c("hnethist", "nethist", "network_histogram"))
})

test_that("a single-layer array is returned as a nethist", {
  # multinethist() serves both cases, so a one-layer fit is not a multinethist.
  arr <- array(A, dim = c(nrow(A), ncol(A), 1L))
  fit <- suppressMessages(multinethist(arr, h = 12))
  expect_equal(class(fit), c("nethist", "network_histogram"))
})

## Shared methods are registered on the base class ----
test_that("logLik and summary are defined once on the base class", {
  registered <- function(generic) {
    tbl <- attr(utils::methods(class = "network_histogram"), "info")
    return(generic %in% tbl$generic)
  }
  expect_true(registered("logLik"))
  expect_true(registered("summary"))
})

test_that("logLik dispatches through the base class for every subclass", {
  for (fit in list(fit_nh, fit_mn, fit_hn)) {
    expect_s3_class(logLik(fit), "logLik")
    expect_false(is.na(attr(logLik(fit), "df")))
    expect_false(is.na(attr(logLik(fit), "nobs")))
  }
})

test_that("summary dispatches through the base class for every subclass", {
  for (fit in list(fit_nh, fit_mn, fit_hn)) {
    expect_s3_class(summary(fit), "summary.network_histogram")
  }
})

test_that("subclass summary methods extend the base method", {
  # summary.multinethist() and summary.hnethist() call NextMethod(), so the
  # fields built by the base method must survive.
  base_fields <- names(summary(fit_nh))
  expect_true(all(base_fields %in% names(summary(fit_mn))))
  expect_true(all(base_fields %in% names(summary(fit_hn))))
  expect_true("n_layers" %in% names(summary(fit_mn)))
  expect_true(all(c("s", "n_candidate") %in% names(summary(fit_hn))))
})

test_that("one print method serves every summary object", {
  for (fit in list(fit_nh, fit_mn, fit_hn)) {
    output <- capture.output(print(summary(fit)))
    expect_true(any(grepl("Model fit:", output)))
  }
})

## Class-specific behaviour is not shared ----
test_that("degrees of freedom follow the subclass, not the base class", {
  k <- length(unique(fit_nh$cluster))
  expect_equal(attr(logLik(fit_nh), "df"), k * (k + 1) / 2)
  expect_equal(attr(logLik(fit_hn), "df"), fit_hn$s)
  k_mn <- length(unique(fit_mn$cluster))
  expect_equal(attr(logLik(fit_mn), "df"), 2L * k_mn * (k_mn + 1) / 2)
})

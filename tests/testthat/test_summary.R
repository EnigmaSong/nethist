# Tests for summary() methods on nethist, multinethist, and hnethist objects.

set.seed(42)
G       <- igraph::sample_gnp(100, p = 0.1)
hist_G  <- nethist(G, h = 10)

set.seed(42)
n_layers <- 3L
A_multi  <- array(0, dim = c(60L, 60L, n_layers))
for (l in seq_len(n_layers)) {
  A_l <- igraph::as_adjacency_matrix(
    igraph::sample_gnp(60L, 0.3, directed = FALSE),
    sparse = FALSE
  )
  diag(A_l) <- 0
  A_multi[, , l] <- A_l
}
hist_mn <- suppressMessages(multinethist(A_multi, h = 12))
hist_hn <- suppressMessages(hnethist(igraph::as_adjacency_matrix(G,
                                                                sparse = FALSE),
                                     h = 10))

## summary.nethist ----
test_that("summary.nethist returns a summary.nethist object", {
  expect_s3_class(summary(hist_G), "summary.nethist")
})

test_that("summary.nethist reports the fit dimensions", {
  s <- summary(hist_G)
  expect_equal(s$n_nodes, length(hist_G$cluster))
  expect_equal(s$n_blocks, length(unique(hist_G$cluster)))
})

test_that("summary.nethist carries logLik, df, AIC, and BIC", {
  s  <- summary(hist_G)
  ll <- logLik(hist_G)
  expect_equal(s$logLik, as.numeric(ll))
  expect_equal(s$df, attr(ll, "df"))
  expect_equal(s$AIC, AIC(hist_G))
  expect_equal(s$BIC, BIC(hist_G))
})

test_that("summary.nethist keeps the fitted object", {
  expect_identical(summary(hist_G)$object, hist_G)
})

test_that("print.summary.nethist does not error", {
  expect_no_error(capture.output(print(summary(hist_G))))
})

test_that("print.summary.nethist returns its input invisibly", {
  s <- summary(hist_G)
  capture.output(result <- withVisible(print(s)))
  expect_false(result$visible)
  expect_identical(result$value, s)
})

test_that("print.summary.nethist reports the fit statistics", {
  output <- capture.output(print(summary(hist_G)))
  expect_true(any(grepl("Log-likelihood", output)))
  expect_true(any(grepl("AIC", output)))
  expect_true(any(grepl("BIC", output)))
  expect_true(any(grepl("rho_hat", output)))
})

test_that("print.summary.nethist respects digits", {
  short <- capture.output(print(summary(hist_G), digits = 3))
  long  <- capture.output(print(summary(hist_G), digits = 10))
  expect_false(identical(short, long))
})

test_that("print.summary.nethist does not print the whole thetahat", {
  # The reviewer asked for compact output; only the dimensions are reported.
  output <- capture.output(print(summary(hist_G)))
  expect_true(any(grepl("probability matrix", output)))
  expect_lt(length(output), 20L)
})

## summary.multinethist ----
test_that("summary.multinethist returns a summary.multinethist object", {
  s <- summary(hist_mn)
  expect_s3_class(s, "summary.multinethist")
  expect_equal(s$n_layers, n_layers)
})

test_that("print.summary.multinethist reports layers and per-layer density", {
  output <- capture.output(print(summary(hist_mn)))
  expect_true(any(grepl("Layers: 3", output)))
  expect_true(any(grepl("rho_hat", output)))
})

test_that("print.summary.multinethist omits MSE when it is unavailable", {
  # multinethist() does not compute an LSE, so no MSE line should be printed.
  skip_if_not(is.null(hist_mn$MSE))
  output <- capture.output(print(summary(hist_mn)))
  expect_false(any(grepl("Mean squared error", output)))
})

## summary.hnethist ----
test_that("summary.hnethist reports the selected and candidate shapes", {
  s <- summary(hist_hn)
  expect_s3_class(s, "summary.hnethist")
  expect_equal(s$s, hist_hn$s)
  expect_equal(s$n_candidate, length(hist_hn$details))
})

test_that("summary.hnethist uses the hnethist degrees of freedom", {
  expect_equal(summary(hist_hn)$df, hist_hn$s)
})

test_that("print.summary.hnethist reports the shape selection", {
  output <- capture.output(print(summary(hist_hn)))
  expect_true(any(grepl("Shape selection", output)))
  expect_true(any(grepl("Selected shapes", output)))
})

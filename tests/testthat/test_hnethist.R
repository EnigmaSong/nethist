# Tests for hnethist()
#
# Note: there is a known bug on line 87 of nethist_hybrid.R:
#   result$rho_hat <- result$initialrho_hat  # should be result$initial$rho_hat
# As a result, hnethist()$rho_hat is always NULL. Tests below reflect this.

set.seed(42)
# Use a small graph to keep test runtime short
G_small   <- igraph::sample_gnp(50, p = 0.15)
A_small   <- igraph::as_adj(G_small, sparse = FALSE)

## Valid input types ----
test_that("hnethist with matrix input does not error", {
  expect_no_error(hnethist(A_small, h = 5, max_itr = 1000))
})

test_that("hnethist with igraph input does not error", {
  expect_no_error(hnethist(G_small, h = 5, max_itr = 1000))
})

## Output class and structure ----
hn <- hnethist(A_small, h = 5, max_itr = 1000)

test_that("return class is hnethist", {
  expect_s3_class(hn, "hnethist")
})

test_that("required top-level fields are present", {
  expect_true(all(c("initial", "details", "blockcluster", "thetahat", "LSE") %in% names(hn)))
})

test_that("initial field is a nethist object", {
  expect_s3_class(hn$initial, "nethist")
})

test_that("details is a list", {
  expect_true(is.list(hn$details))
})

test_that("thetahat is a matrix", {
  expect_true(is.matrix(hn$thetahat))
})

test_that("thetahat is square with dimension equal to initial k", {
  k <- max(hn$initial$cluster)
  expect_equal(nrow(hn$thetahat), k)
  expect_equal(ncol(hn$thetahat), k)
})

test_that("LSE is a non-negative scalar", {
  expect_true(is.numeric(hn$LSE))
  expect_true(hn$LSE >= 0)
})

test_that("each details entry has required sub-fields", {
  for (d in hn$details) {
    expect_true(all(c("n", "s", "cluster", "thetahat", "LSE", "BIC") %in% names(d)))
  }
})

## Methods ----
test_that("hnethist with method = 'PLL' does not error", {
  expect_no_error(hnethist(A_small, h = 5, method = "PLL", max_itr = 1000))
})

test_that("hnethist with method = 'LSE' does not error", {
  expect_no_error(hnethist(A_small, h = 5, method = "LSE", max_itr = 1000))
})

## print.hnethist ----
test_that("print.hnethist does not error", {
  expect_no_error(print(hn))
})

test_that("print.hnethist returns object invisibly", {
  result <- withVisible(print(hn))
  expect_false(result$visible)
})

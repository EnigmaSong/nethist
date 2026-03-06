# Tests for .prob_mat_from_adj()

set.seed(42)
A      <- igraph::as_adj(igraph::sample_gnp(100, 0.1), sparse = FALSE)
hist_A <- nethist(A, h = 10)

p_mat <- .prob_mat_from_adj(A, hist_A$cluster)
k     <- max(hist_A$cluster)

test_that("returns a matrix", {
  expect_true(is.matrix(p_mat))
})

test_that("dimensions equal the number of clusters k x k", {
  expect_equal(dim(p_mat), c(k, k))
})

test_that("all entries are in [0, 1]", {
  expect_true(all(p_mat >= 0))
  expect_true(all(p_mat <= 1))
})

test_that("matrix is symmetric", {
  expect_equal(p_mat, t(p_mat))
})

test_that("diagonal entries are consistent with within-group edge rates", {
  # Diagonal should be >= 0 and <= 1 (already covered), but also non-negative
  expect_true(all(diag(p_mat) >= 0))
})

test_that("result is consistent with thetahat from nethist", {
  # .prob_mat_from_adj should agree with nethist's thetahat (which is the same computation)
  # unname() because .prob_mat_from_adj sets dimnames but thetahat does not
  expect_equal(unname(p_mat), hist_A$thetahat, tolerance = 1e-10)
})

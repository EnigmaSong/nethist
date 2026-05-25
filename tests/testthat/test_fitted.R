# -- Fixtures ----------------------------------------------------------------
withr::local_seed(42)
n  <- 60L
h  <- 10L
A  <- igraph::as_adjacency_matrix(
  igraph::sample_gnp(n, 0.3, directed = FALSE), sparse = FALSE)
diag(A) <- 0L

fit_n  <- suppressMessages(nethist(A, h = h))
fit_hn <- suppressMessages(hnethist(A, h = h))

L <- 4L
A3 <- array(0L, c(n, n, L))
for (l in seq_len(L)) {
  tmp <- igraph::as_adjacency_matrix(
    igraph::sample_gnp(n, 0.2 + 0.1 * l, directed = FALSE),
    sparse = FALSE)
  diag(tmp) <- 0L
  A3[, , l] <- tmp
}
fit_m <- suppressMessages(multinethist(A3, h = h))

k  <- length(unique(fit_n$cluster))
kL <- dim(fit_m$thetahat)

# -- fitted.nethist: return dimensions ---------------------------------------
test_that("fitted.nethist default returns n x n matrix", {
  mat <- fitted(fit_n)
  expect_equal(dim(mat), c(n, n))
})

test_that("fitted.nethist set1/set2 returns correct submatrix dimensions", {
  mat <- fitted(fit_n, set1 = 1:10, set2 = 1:20)
  expect_equal(dim(mat), c(10L, 20L))
})

test_that("fitted.nethist single-node set returns 1 x m matrix", {
  mat <- fitted(fit_n, set1 = 1L, set2 = 1:5)
  expect_equal(dim(mat), c(1L, 5L))
})

# -- fitted.nethist: type ----------------------------------------------------
test_that("fitted.nethist type='nethist' returns thetahat / rho_hat", {
  mat <- fitted(fit_n, type = "nethist")
  expected <- fit_n$thetahat / fit_n$rho_hat
  expect_equal(mat, expected[fit_n$cluster, fit_n$cluster])
})

test_that("fitted.nethist type='prob' returns thetahat", {
  mat <- fitted(fit_n, type = "prob")
  expected <- fit_n$thetahat
  expect_equal(mat, expected[fit_n$cluster, fit_n$cluster])
})

test_that("fitted.nethist invalid type raises error", {
  expect_error(fitted(fit_n, type = "graphon"))
})

# -- fitted.nethist: values in range -----------------------------------------
test_that("fitted.nethist type='prob' values are in [0, 1]", {
  mat <- fitted(fit_n, type = "prob")
  expect_true(all(mat >= 0 & mat <= 1))
})

test_that("fitted.nethist type='nethist' values are non-negative", {
  mat <- fitted(fit_n, type = "nethist")
  expect_true(all(mat >= 0))
})

# -- fitted.nethist: hnethist dispatch ---------------------------------------
test_that("fitted dispatches to fitted.nethist for hnethist objects", {
  mat <- fitted(fit_hn)
  expect_equal(dim(mat), c(n, n))
  expect_true(all(mat >= 0))
})

# -- fitted.multinethist: return dimensions ----------------------------------
test_that("fitted.multinethist single layer with drop=TRUE returns matrix", {
  mat <- fitted(fit_m, layer = 1L)
  expect_equal(dim(mat), c(n, n))
  expect_true(is.matrix(mat))
})

test_that("fitted.multinethist single layer with drop=FALSE returns array", {
  arr <- fitted(fit_m, layer = 1L, drop = FALSE)
  expect_equal(dim(arr), c(n, n, 1L))
})

test_that("fitted.multinethist multiple layers returns 3D array", {
  arr <- fitted(fit_m, layer = c(1L, 2L))
  expect_equal(dim(arr), c(n, n, 2L))
})

test_that("fitted.multinethist set1/set2/layer all specified", {
  arr <- fitted(fit_m, set1 = 1:10, set2 = 1:15, layer = c(1L, 3L))
  expect_equal(dim(arr), c(10L, 15L, 2L))
})

# -- fitted.multinethist: type -----------------------------------------------
test_that("fitted.multinethist type='prob' returns raw thetahat slice", {
  mat <- fitted(fit_m, layer = 1L, type = "prob")
  expected <- fit_m$thetahat[, , 1L]
  expect_equal(mat, expected[fit_m$cluster, fit_m$cluster])
})

test_that("fitted.multinethist type='nethist' divides by rho_hat", {
  mat <- fitted(fit_m, layer = 1L, type = "nethist")
  expected <- fit_m$thetahat[, , 1L] / fit_m$rho_hat[1L]
  expect_equal(mat, expected[fit_m$cluster, fit_m$cluster])
})

test_that("fitted.multinethist type='prob' values are in [0, 1]", {
  arr <- fitted(fit_m, type = "prob")
  expect_true(all(arr >= 0 & arr <= 1))
})

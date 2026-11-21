# Tests for reproducibility across runs and platforms.
#
# Fitted partitions may differ across platforms because of floating-point
# rounding, so these tests do not compare partitions with fixed values. They
# check that a fit is reproducible on the same platform with the same seed, and
# that quantities which do not depend on the fitted partition follow from the
# input alone.

# -- Fixtures ----------------------------------------------------------------
set.seed(42)
n_v <- 60L
h_v <- 12L
A <- igraph::as_adjacency_matrix(
  igraph::sample_gnp(n_v, 0.3, directed = FALSE),
  sparse = FALSE
)
diag(A) <- 0

set.seed(42)
n_layers <- 3L
A_multi <- array(0, dim = c(n_v, n_v, n_layers))
for (l in seq_len(n_layers)) {
  A_l <- igraph::as_adjacency_matrix(
    igraph::sample_gnp(n_v, 0.3, directed = FALSE),
    sparse = FALSE
  )
  diag(A_l) <- 0
  A_multi[, , l] <- A_l
}

fit_twice <- function(seed, fit_fun) {
  first  <- withr::with_seed(seed, suppressMessages(fit_fun()))
  second <- withr::with_seed(seed, suppressMessages(fit_fun()))
  return(list(first = first, second = second))
}

# Group sizes implied by n vertices and bandwidth h: equal groups of size h,
# with the remainder added to the last group.
expected_sizes <- function(n, h) {
  n_groups <- n %/% h
  return(sort(c(rep(h, n_groups - 1), h + n - n_groups * h)))
}

# -- Reproducibility on the same platform -------------------------------------

test_that("nethist() gives identical fits with the same seed", {
  for (method in c("PLL", "LSE")) {
    fits <- fit_twice(42, function() nethist(A, h = h_v, method = method))
    expect_identical(block_id(fits$first), block_id(fits$second))
    expect_identical(coef(fits$first, type = "prob"),
                     coef(fits$second, type = "prob"))
    expect_identical(logLik(fits$first, normalized = TRUE),
                     logLik(fits$second, normalized = TRUE))
  }
})

test_that("hnethist() gives identical fits with the same seed", {
  fits <- fit_twice(42, function() hnethist(A, h = h_v))
  expect_identical(block_id(fits$first), block_id(fits$second))
  expect_identical(shape_id(fits$first), shape_id(fits$second))
  expect_identical(BIC(fits$first), BIC(fits$second))
})

test_that("multinethist() gives identical fits with the same seed", {
  for (common_f in c(FALSE, TRUE)) {
    fits <- fit_twice(42, function() {
      multinethist(A_multi, h = h_v, common_f = common_f)
    })
    expect_identical(block_id(fits$first), block_id(fits$second))
    expect_identical(coef(fits$first, type = "prob"),
                     coef(fits$second, type = "prob"))
    expect_identical(logLik(fits$first, normalized = TRUE),
                     logLik(fits$second, normalized = TRUE))
  }
})

test_that("data-driven bandwidth selection is reproducible", {
  fits <- fit_twice(42, function() nethist(A))
  expect_identical(bandwidth(fits$first), bandwidth(fits$second))
  expect_identical(block_id(fits$first), block_id(fits$second))
})

# -- Quantities that do not depend on the fitted partition --------------------

test_that("group sizes follow from the number of vertices and the bandwidth", {
  fit_nh  <- withr::with_seed(42, suppressMessages(nethist(A, h = h_v)))
  fit_mnh <- withr::with_seed(42, suppressMessages(
    multinethist(A_multi, h = h_v)))
  fit_sel <- withr::with_seed(42, suppressMessages(nethist(A)))
  for (fit in list(fit_nh, fit_mnh, fit_sel)) {
    sizes <- sort(as.vector(table(block_id(fit))))
    expect_equal(sizes, expected_sizes(n_v, bandwidth(fit)))
  }
})

test_that("shape labels of hnethist are 1 to s with no gaps", {
  fit <- withr::with_seed(42, suppressMessages(hnethist(A, h = h_v)))
  ids <- shape_id(fit)
  expect_setequal(unique(as.vector(ids)), seq_len(nrow(summary(fit)$shapes)))
})

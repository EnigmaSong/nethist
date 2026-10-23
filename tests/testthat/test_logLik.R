# Tests for logLik() methods on nethist, multinethist, and hnethist objects.

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

fit_nh  <- suppressMessages(nethist(A, h = h_v))
fit_mnh <- suppressMessages(multinethist(A_multi, h = h_v))
fit_cf  <- suppressMessages(multinethist(A_multi, h = h_v, common_f = TRUE))
fit_hnh <- suppressMessages(hnethist(A, h = h_v))

# Bernoulli log-likelihood of a single layer at the fitted probabilities.
layer_LL <- function(A_l, thetahat, cluster) {
  theta_mat <- thetahat[cluster, cluster]
  eps       <- .Machine$double.eps
  theta_mat <- pmax(pmin(theta_mat, 1 - eps), eps)
  upper     <- upper.tri(A_l)
  return(sum(A_l[upper] * log(theta_mat[upper]) +
             (1 - A_l[upper]) * log(1 - theta_mat[upper])))
}

n_pairs <- n_v * (n_v - 1) / 2

# -- nethist -----------------------------------------------------------------
test_that("logLik.nethist returns a logLik object", {
  ll <- logLik(fit_nh)
  expect_s3_class(ll, "logLik")
  expect_length(as.numeric(ll), 1L)
})

test_that("logLik.nethist equals the Bernoulli log-likelihood", {
  expect_equal(as.numeric(logLik(fit_nh)),
               layer_LL(A, fit_nh$thetahat, fit_nh$cluster))
})

test_that("logLik.nethist has df = k(k+1)/2 and nobs = n(n-1)/2", {
  ll <- logLik(fit_nh)
  k  <- length(unique(fit_nh$cluster))
  expect_equal(attr(ll, "df"), k * (k + 1) / 2)
  expect_equal(attr(ll, "nobs"), n_pairs)
})

test_that("normalized = TRUE returns the stored per-edge value", {
  expect_identical(logLik(fit_nh, normalized = TRUE), fit_nh$normalized_LL)
})

test_that("normalized = TRUE is not classed as logLik", {
  # Guards against AIC()/BIC() silently accepting a normalized value.
  expect_false(inherits(logLik(fit_nh, normalized = TRUE), "logLik"))
})

test_that("unnormalized value is the normalized value times the edge count", {
  m <- sum(A) / 2
  expect_equal(as.numeric(logLik(fit_nh)),
               logLik(fit_nh, normalized = TRUE) * m)
})

test_that("AIC and BIC agree with the definition", {
  ll <- logLik(fit_nh)
  expect_equal(AIC(fit_nh), -2 * as.numeric(ll) + 2 * attr(ll, "df"))
  expect_equal(BIC(fit_nh),
               -2 * as.numeric(ll) + attr(ll, "df") * log(attr(ll, "nobs")))
})

test_that("logLik works for method = 'LSE' fits", {
  fit_lse <- suppressMessages(nethist(A, h = h_v, method = "LSE"))
  expect_s3_class(logLik(fit_lse), "logLik")
  expect_equal(as.numeric(logLik(fit_lse)),
               layer_LL(A, fit_lse$thetahat, fit_lse$cluster))
})

# -- multinethist ------------------------------------------------------------
test_that("logLik.multinethist equals the summed layer log-likelihoods", {
  expected <- sum(vapply(seq_len(n_layers), function(l) {
    layer_LL(A_multi[, , l], fit_mnh$thetahat[, , l], fit_mnh$cluster)
  }, numeric(1)))
  expect_equal(as.numeric(logLik(fit_mnh)), expected)
})

test_that("logLik.multinethist counts one parameter set per layer", {
  ll <- logLik(fit_mnh)
  k  <- length(unique(fit_mnh$cluster))
  expect_equal(attr(ll, "df"), n_layers * k * (k + 1) / 2)
  expect_equal(attr(ll, "nobs"), n_layers * n_pairs)
})

test_that("common_f = TRUE shares one parameter set across layers", {
  ll <- logLik(fit_cf)
  k  <- length(unique(fit_cf$cluster))
  expect_equal(attr(ll, "df"), k * (k + 1) / 2)
  expect_equal(attr(ll, "nobs"), n_layers * n_pairs)
})

test_that("logLik.multinethist normalized value matches the stored field", {
  expect_identical(logLik(fit_mnh, normalized = TRUE), fit_mnh$normalized_LL)
})

# -- hnethist ----------------------------------------------------------------
test_that("logLik.hnethist uses the number of selected shapes as df", {
  ll <- logLik(fit_hnh)
  expect_equal(attr(ll, "df"), fit_hnh$s)
  expect_equal(attr(ll, "nobs"), n_pairs)
})

test_that("logLik.hnethist equals the Bernoulli log-likelihood", {
  expect_equal(as.numeric(logLik(fit_hnh)),
               layer_LL(A, fit_hnh$thetahat, fit_hnh$cluster))
})

test_that("BIC() reproduces the BIC used for shape selection", {
  expect_equal(BIC(fit_hnh), fit_hnh$BIC)
})

# -- Degenerate single-block fits --------------------------------------------
test_that("logLik works when multinethist returns a single block", {
  # h >= n/2 makes multinethist() return thetahat as a vector, not an array.
  set.seed(42)
  n_s <- 20L
  A_s <- array(0, dim = c(n_s, n_s, n_layers))
  for (l in seq_len(n_layers)) {
    A_l <- igraph::as_adjacency_matrix(
      igraph::sample_gnp(n_s, 0.3, directed = FALSE),
      sparse = FALSE
    )
    diag(A_l) <- 0
    A_s[, , l] <- A_l
  }
  fit_s <- suppressWarnings(suppressMessages(multinethist(A_s, h = n_s)))
  ll    <- logLik(fit_s)
  expect_equal(attr(ll, "df"), n_layers)
  expect_equal(attr(ll, "nobs"), n_layers * n_s * (n_s - 1) / 2)
  expect_false(is.na(as.numeric(ll)))
})

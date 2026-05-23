# -- Fixture -----------------------------------------------------------------
set.seed(42)
n <- 50L
h <- 10L
A <- igraph::as_adjacency_matrix(
  igraph::sample_gnp(n, 0.4, directed = FALSE),
  sparse = FALSE
)
diag(A) <- 0

# -- Input validation --------------------------------------------------------
test_that("self-loop raises error", {
  A_sl <- A
  diag(A_sl) <- 1
  expect_error(nethist:::hnethist(A_sl, h = h))
})

test_that("asymmetric matrix raises error", {
  A_asym <- A
  A_asym[1, 2] <- 1 - A_asym[1, 2]
  expect_error(nethist:::hnethist(A_asym, h = h))
})

test_that("unsupported input type raises error", {
  expect_error(nethist:::hnethist(as.data.frame(A), h = h))
})

# -- Valid input types -------------------------------------------------------
test_that("matrix input runs without error", {
  expect_no_error(suppressMessages(nethist:::hnethist(A, h = h)))
})

test_that("igraph input runs without error", {
  g <- igraph::graph_from_adjacency_matrix(A, mode = "undirected")
  expect_no_error(suppressMessages(nethist:::hnethist(g, h = h)))
})

# -- Output structure --------------------------------------------------------
result <- suppressMessages(nethist:::hnethist(A, h = h))

test_that("class is c('hnethist', 'nethist')", {
  expect_equal(class(result), c("hnethist", "nethist"))
})

test_that("all required fields exist", {
  expect_named(
    result,
    c("cluster", "thetahat", "rho_hat", "normalized_LL", "MSE", "method",
      "blockcluster", "BIC", "s", "details", "initial"),
    ignore.order = TRUE
  )
})

test_that("thetahat is symmetric", {
  expect_true(isSymmetric(result$thetahat))
})

test_that("thetahat values are in [0, 1]", {
  expect_true(all(result$thetahat >= 0 & result$thetahat <= 1))
})

test_that("cluster length equals n", {
  expect_length(result$cluster, n)
})

test_that("details length equals s_max", {
  k     <- length(unique(result$initial$cluster))
  th    <- result$initial$thetahat
  s_max <- min(k * (k + 1L) / 2L,
               length(unique(th[upper.tri(th, diag = TRUE)])))
  expect_length(result$details, s_max)
})

test_that("initial field is a nethist object", {
  expect_s3_class(result$initial, "nethist")
})

# -- Model selection ---------------------------------------------------------
test_that("selected model minimizes BIC over all candidates", {
  BICs <- sapply(result$details, function(x) x$BIC)
  expect_equal(result$BIC, min(BICs))
  expect_equal(result$MSE, result$details[[which.min(BICs)]]$MSE)
})

test_that("s is a positive integer-valued scalar", {
  expect_length(result$s, 1L)
  expect_gte(result$s, 1L)
  expect_equal(result$s, as.integer(result$s))
})

# -- plot.hnethist -----------------------------------------------------------
test_that("plot(result) draws heatmap without error", {
  expect_no_error(plot(result))
})

test_that("plot(result, type='BIC') runs without error", {
  expect_no_error(plot(result, type = "BIC"))
})

test_that("plot(result, type='bic') accepts lowercase", {
  expect_no_error(plot(result, type = "bic"))
})

test_that("plot(result, type='BIC') returns NULL invisibly", {
  expect_null(plot(result, type = "BIC"))
})

test_that("plot(result, type='invalid') raises error", {
  expect_error(plot(result, type = "invalid"))
})

# -- print -------------------------------------------------------------------
test_that("print() dispatches via nethist inheritance without error", {
  expect_no_error(print(result))
})

# -- Helper: hnethist_BIC ----------------------------------------------------
test_that("hnethist_BIC returns correct Bernoulli BIC value", {
  nh       <- list(n = 100L, s = 3L)
  ll       <- -500.0
  N        <- 100L * 99L / 2L
  expected <- -2 * ll + 3L * log(N)
  expect_equal(nethist:::hnethist_BIC(nh, ll), expected)
})

test_that("hnethist_BIC penalty increases with s given same ll", {
  make_nh <- function(s) list(n = 100L, s = s)
  ll <- -500.0
  expect_lt(nethist:::hnethist_BIC(make_nh(1L), ll),
            nethist:::hnethist_BIC(make_nh(5L), ll))
})

test_that("hnethist_BIC decreases with better ll given same s", {
  nh <- list(n = 100L, s = 3L)
  expect_gt(nethist:::hnethist_BIC(nh, -1000.0),
            nethist:::hnethist_BIC(nh, -100.0))
})

# -- Helper: hnethist_LSE ----------------------------------------------------
test_that("hnethist_LSE is non-negative", {
  lse <- nethist:::hnethist_LSE(result$thetahat, result$cluster, A)
  expect_gte(lse, 0)
})

test_that("hnethist_LSE is 0 for perfect fit on zero matrix", {
  n_s  <- 6L
  A_z  <- matrix(0, n_s, n_s)
  th_z <- matrix(0.0, 1L, 1L)
  cl_z <- rep(1L, n_s)
  expect_equal(nethist:::hnethist_LSE(th_z, cl_z, A_z), 0)
})

test_that("hnethist_LSE equals raw sum divided by n^2", {
  n_s  <- 4L
  A_s  <- matrix(c(0, 1, 0, 1,
                   1, 0, 1, 0,
                   0, 1, 0, 1,
                   1, 0, 1, 0), n_s, n_s)
  th_s <- matrix(c(0.2, 0.6, 0.6, 0.2), 2L, 2L)
  cl_s <- c(1L, 2L, 1L, 2L)
  expected <- sum((A_s - th_s[cl_s, cl_s])^2) / n_s^2
  expect_equal(nethist:::hnethist_LSE(th_s, cl_s, A_s), expected)
})

# -- Helper: hnethist_normalized_LL ------------------------------------------
test_that("hnethist_normalized_LL is non-positive", {
  nll <- nethist:::hnethist_normalized_LL(result$thetahat, result$cluster, A)
  expect_lte(nll, 0)
})

test_that("hnethist_normalized_LL equals Bernoulli LL / sum(A)", {
  n_s  <- 4L
  A_s  <- matrix(c(0, 1, 0, 1,
                   1, 0, 1, 0,
                   0, 1, 0, 1,
                   1, 0, 1, 0), n_s, n_s)
  th_s <- matrix(c(0.2, 0.6, 0.6, 0.2), 2L, 2L)
  cl_s <- c(1L, 2L, 1L, 2L)
  theta_mat <- th_s[cl_s, cl_s]
  upper     <- upper.tri(A_s)
  expected  <- sum(A_s[upper] * log(theta_mat[upper]) +
                   (1 - A_s[upper]) * log(1 - theta_mat[upper])) / sum(A_s)
  expect_equal(nethist:::hnethist_normalized_LL(th_s, cl_s, A_s), expected)
})

# -- normalized_LL from best model -------------------------------------------
test_that("normalized_LL equals best details LL", {
  BICs <- sapply(result$details, function(x) x$BIC)
  expect_equal(result$normalized_LL,
               result$details[[which.min(BICs)]]$normalized_LL)
})

test_that("normalized_LL is a length-1 negative scalar", {
  expect_length(result$normalized_LL, 1L)
  expect_lt(result$normalized_LL, 0)
})

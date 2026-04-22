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
    c("cluster", "thetahat", "rho_hat", "normalized_LL", "LSE", "method",
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
  s_max <- min(k * (k - 1L) / 2L,
               length(unique(c(result$initial$thetahat))))
  expect_length(result$details, s_max)
})

test_that("initial field is a nethist object", {
  expect_s3_class(result$initial, "nethist")
})

# -- Model selection ---------------------------------------------------------
test_that("selected model minimizes BIC over all candidates", {
  BICs <- sapply(result$details, function(x) x$BIC)
  expect_equal(result$BIC, min(BICs))
  expect_equal(result$LSE, result$details[[which.min(BICs)]]$LSE)
})

test_that("s is a positive integer-valued scalar", {
  expect_length(result$s, 1L)
  expect_gte(result$s, 1L)
  expect_equal(result$s, as.integer(result$s))
})

# -- S3 method inheritance ---------------------------------------------------
test_that("plot() dispatches via nethist inheritance without error", {
  expect_no_error(plot(result))
})

test_that("print() dispatches via nethist inheritance without error", {
  expect_no_error(print(result))
})

# -- Helper: hnethist_BIC ----------------------------------------------------
test_that("hnethist_BIC returns correct Gaussian BIC value", {
  nh       <- list(n = 100L, s = 3L, LSE = 50.0)
  N        <- 100L * 99L / 2L
  expected <- N * log(50.0 / (2 * N)) + 3L * log(N)
  expect_equal(nethist:::hnethist_BIC(nh), expected)
})

test_that("hnethist_BIC penalty increases with s given same LSE", {
  make_nh <- function(s) list(n = 100L, s = s, LSE = 50.0)
  expect_lt(nethist:::hnethist_BIC(make_nh(1L)),
            nethist:::hnethist_BIC(make_nh(5L)))
})

test_that("hnethist_BIC decreases with smaller LSE given same s", {
  make_nh <- function(lse) list(n = 100L, s = 3L, LSE = lse)
  expect_gt(nethist:::hnethist_BIC(make_nh(100.0)),
            nethist:::hnethist_BIC(make_nh(10.0)))
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

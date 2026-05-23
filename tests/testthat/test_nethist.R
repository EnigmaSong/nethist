# Tests for nethist()

set.seed(42)
G        <- igraph::sample_gnp(100, p = 0.1)
A_dense  <- igraph::as_adjacency_matrix(G, sparse = FALSE)
A_sparse <- igraph::as_adjacency_matrix(G)
h_used   <- 10

## Input validation ----
test_that("self-loop input errors", {
  A <- A_dense
  A[1, 1] <- 1
  expect_error(nethist(A))
})

test_that("asymmetric matrix input errors", {
  A <- A_dense
  A[1, 2] <- 1; A[2, 1] <- 0
  expect_error(nethist(A))
})

test_that("non-matrix input errors", {
  expect_error(nethist(list()))
})

test_that("invalid method errors", {
  expect_error(nethist(G, h = h_used, method = "foo"))
})

## Valid input types ----
test_that("igraph input works", {
  expect_no_error(nethist(G, h = h_used))
})

test_that("dense matrix input works", {
  expect_no_error(nethist(A_dense, h = h_used))
})

test_that("sparse matrix input works", {
  expect_no_error(nethist(A_sparse, h = h_used))
})

## Output structure ----
hist_G <- nethist(G, h = h_used)

test_that("return class is nethist", {
  expect_s3_class(hist_G, "nethist")
})

test_that("required fields present", {
  expect_named(hist_G,
               c("cluster", "thetahat", "rho_hat", "normalized_LL",
                 "MSE", "method"),
               ignore.order = TRUE)
})

test_that("thetahat is symmetric", {
  expect_equal(hist_G$thetahat, t(hist_G$thetahat))
})

test_that("thetahat entries are in [0, 1]", {
  expect_true(all(hist_G$thetahat >= 0))
  expect_true(all(hist_G$thetahat <= 1))
})

test_that("rho_hat is in [0, 1]", {
  expect_true(hist_G$rho_hat >= 0)
  expect_true(hist_G$rho_hat <= 1)
})

test_that("cluster vector length matches number of nodes", {
  expect_equal(length(hist_G$cluster), igraph::vcount(G))
})

test_that("cluster labels are positive integers up to k", {
  k <- max(hist_G$cluster)
  expect_true(all(hist_G$cluster %in% seq_len(k)))
})

test_that("all cluster labels appear at least once", {
  k <- max(hist_G$cluster)
  expect_equal(sort(unique(hist_G$cluster)), seq_len(k))
})

test_that("thetahat dimensions equal number of clusters", {
  k <- max(hist_G$cluster)
  expect_equal(dim(hist_G$thetahat), c(k, k))
})

test_that("method field is 'PLL' by default", {
  expect_equal(hist_G$method, "PLL")
})

test_that("method field is 'LSE' when specified", {
  hist_lse <- nethist(G, h = h_used, method = "LSE")
  expect_equal(hist_lse$method, "LSE")
})

## Bandwidth ----
test_that("auto bandwidth selection (h = NA) works", {
  expect_no_error(nethist(G))
})

test_that("bin sizes do not exceed h (except the last bin)", {
  bin_size <- as.integer(table(hist_G$cluster))
  bin_size_excl_last <- sort(bin_size)[-length(bin_size)]
  expect_true(all(bin_size_excl_last <= h_used))
})

## Methods ----
test_that("default method equals explicit PLL", {
  expect_equal({set.seed(42); nethist(G, h = h_used)},
               {set.seed(42); nethist(G, h = h_used, method = "PLL")})
})

## Reproducibility ----
test_that("same seed gives identical results", {
  expect_equal({set.seed(1); nethist(G, h = h_used)},
               {set.seed(1); nethist(G, h = h_used)})
})

## Edge case: n == h ----
test_that("n == h (single-bin) returns without error", {
  n <- igraph::vcount(G)
  suppressMessages(expect_no_error(nethist(G, h = n)))
})

test_that("n == h (single-bin) returns k == 1", {
  n <- igraph::vcount(G)
  result <- suppressMessages(nethist(G, h = n))
  expect_equal(max(result$cluster), 1L)
})

## Bandwidth warnings ----
test_that("h < 2 triggers warning and is clamped to 2", {
  expect_warning(nethist(A_dense, h = 1), regexp = "adjusted")
  expect_warning(nethist(A_dense, h = 0), regexp = "adjusted")
})

test_that("h > n triggers warning and is clamped to n", {
  n <- nrow(A_dense)
  expect_warning(nethist(A_dense, h = n + 1), regexp = "adjusted")
})

test_that("non-integer h triggers warning", {
  expect_warning(nethist(A_dense, h = 2.9), regexp = "adjusted")
  expect_warning(nethist(A_dense, h = 3.1), regexp = "adjusted")
})

test_that("valid integer h in range triggers no warning", {
  expect_no_warning(nethist(A_dense, h = 4))
})

test_that("h = NA (auto-select) triggers no warning", {
  # auto-selected h may be non-integer internally, but should not warn the user
  expect_no_warning(nethist(A_dense, h = NA))
})

## Edge cases: zero and complete graph ----
test_that("zero matrix (no edges) errors gracefully in auto-bandwidth", {
  # rho_hat = 0 causes numerical issues in bandwidth selection
  n <- 20
  A_zero <- matrix(0L, n, n)
  expect_error(nethist(A_zero))
})

test_that("complete graph runs without error and rho_hat is close to 1", {
  n <- 20
  A_full <- matrix(1L, n, n)
  diag(A_full) <- 0L
  result_full <- nethist(A_full)
  expect_true(result_full$rho_hat > 0.9)
})

## Invalid method type ----
test_that("numeric method argument triggers error", {
  expect_error(nethist(G, h = h_used, method = 123))
})

## network input ----
test_that("network input gives identical result to matrix", {
  skip_if_not_installed("network")
  A_net <- network::as.network(A_dense, directed = FALSE)
  expect_equal({set.seed(1L); nethist(A_net,   h = h_used)},
               {set.seed(1L); nethist(A_dense, h = h_used)})
})

test_that("combined_networks input errors with redirect to multinethist", {
  skip_if_not_installed("network")
  skip_if_not_installed("ergm.multi")
  net  <- network::as.network(A_dense, directed = FALSE)
  nets <- ergm.multi::Networks(list(net, net))
  expect_error(nethist(nets, h = h_used), regexp = "multinethist")
})

# Tests for plot.nethist() and plot.multinethist()

set.seed(42)
G      <- igraph::sample_gnp(100, p = 0.1)
hist_G <- nethist(G, h = 10)
k      <- max(hist_G$cluster)

## plot.nethist ----
test_that("plot.nethist with default arguments", {
  expect_no_error(plot(hist_G))
})

test_that("plot.nethist with type = 'prob'", {
  expect_no_error(plot(hist_G, type = "prob"))
})

test_that("plot.nethist with type = 'pmat' (alias for prob)", {
  expect_no_error(plot(hist_G, type = "pmat"))
})

test_that("plot.nethist with prob = TRUE overlays values", {
  expect_no_error(plot(hist_G, type = "prob", prob = TRUE))
})

test_that("plot.nethist with valid idx_order", {
  idx <- unique(hist_G$cluster)
  expect_no_error(plot(hist_G, idx_order = idx))
})

test_that("plot.nethist with invalid idx_order produces warning then plots", {
  expect_warning(plot(hist_G, idx_order = seq_len(k + 1)))
})

test_that("plot.nethist with invalid type errors", {
  expect_error(plot(hist_G, type = "invalid"))
})

## plot.multinethist ----
set.seed(42)
data(kite, package = "igraphdata")
kite_mat <- igraph::as_adjacency_matrix(igraph::upgrade_graph(kite), sparse = FALSE)
kite_arr <- array(kite_mat, dim = c(10, 10, 1))
mnhist   <- multinethist(kite_arr, h = 5, max_itr = 1000)

test_that("plot.multinethist with default arguments", {
  expect_no_error(plot(mnhist))
})

test_that("plot.multinethist with explicit idx_order", {
  idx <- unique(mnhist$cluster)
  expect_no_error(plot(mnhist, idx_order = idx))
})

test_that("plot.multinethist with invalid idx_order warns", {
  k_mn <- max(mnhist$cluster)
  expect_warning(plot(mnhist, idx_order = seq_len(k_mn + 1)))
})

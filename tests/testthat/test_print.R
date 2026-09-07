# Tests for print.nethist() and print.multinethist()

set.seed(42)
G      <- igraph::sample_gnp(100, p = 0.1)
hist_G <- nethist(G, h = 10)

data(kite, package = "igraphdata")
kite_mat <- igraph::as_adjacency_matrix(igraph::upgrade_graph(kite), sparse = FALSE)
kite_arr <- array(kite_mat, dim = c(10, 10, 1))
set.seed(42)
# A one-layer array is returned as a nethist, so a genuine multilayer fixture
# is needed to reach print.multinethist().
kite_nh <- multinethist(kite_arr, h = 5, control = nethist_control(max_itr = 1000))

set.seed(42)
mn_arr <- array(0, dim = c(60L, 60L, 2L))
for (l in 1:2) {
  A_l <- igraph::as_adjacency_matrix(
    igraph::sample_gnp(60L, 0.3, directed = FALSE),
    sparse = FALSE
  )
  diag(A_l) <- 0
  mn_arr[, , l] <- A_l
}
kite_mn <- suppressMessages(multinethist(mn_arr, h = 12))

## print.nethist ----
test_that("print.nethist does not error", {
  expect_no_error(capture.output(print(hist_G)))
})

test_that("print.nethist returns object invisibly", {
  capture.output(result <- withVisible(print(hist_G)))
  expect_false(result$visible)
  expect_identical(result$value, hist_G)
})

test_that("print.nethist reports thetahat by its dimensions only", {
  output <- capture.output(print(hist_G))
  k      <- length(unique(hist_G$cluster))
  expect_true(any(grepl(paste0("thetahat: ", k, " x ", k), output)))
  # The estimated probabilities belong in summary(), not in print().
  expect_false(any(grepl("Available components", output)))
  expect_lt(length(output), 12L)
})

test_that("print.nethist reports the fit identity", {
  output <- capture.output(print(hist_G))
  expect_true(any(grepl("Network histogram fit", output)))
  expect_true(any(grepl(paste("Nodes:", length(hist_G$cluster)), output)))
  expect_true(any(grepl("Method:", output)))
  expect_true(any(grepl("summary\\(\\)", output)))
})

test_that("print.nethist respects digits", {
  short <- capture.output(print(hist_G, digits = 3))
  long  <- capture.output(print(hist_G, digits = 10))
  expect_false(identical(short, long))
})

## print.multinethist ----
test_that("print.multinethist does not error", {
  expect_no_error(capture.output(print(kite_mn)))
})

test_that("print.multinethist returns object invisibly", {
  capture.output(result <- withVisible(print(kite_mn)))
  expect_false(result$visible)
  expect_identical(result$value, kite_mn)
})

test_that("print.multinethist reports thetahat by its dimensions only", {
  output <- capture.output(print(kite_mn))
  expect_true(any(grepl("thetahat: .* probability", output)))
  expect_false(any(grepl("Available components", output)))
  expect_lt(length(output), 12L)
})

test_that("print.multinethist reports layers and the graphon assumption", {
  output <- capture.output(print(kite_mn))
  expect_true(any(grepl("Layers: 2", output)))
  expect_true(any(grepl("Common graphon across layers:", output)))
})

## print.hnethist ----
test_that("print.hnethist reports the selected shapes", {
  set.seed(42)
  hn     <- suppressMessages(hnethist(
    igraph::as_adjacency_matrix(G, sparse = FALSE), h = 10))
  output <- capture.output(print(hn))
  expect_true(any(grepl("Hybrid network histogram fit", output)))
  expect_true(any(grepl(paste("Selected shapes:", hn$s), output)))
  expect_lt(length(output), 12L)
})

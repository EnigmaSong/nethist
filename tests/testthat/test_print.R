# Tests for print.nethist() and print.multinethist()

set.seed(42)
G      <- igraph::sample_gnp(100, p = 0.1)
hist_G <- nethist(G, h = 10)

data(kite, package = "igraphdata")
kite_mat <- igraph::as_adjacency_matrix(igraph::upgrade_graph(kite), sparse = FALSE)
kite_arr <- array(kite_mat, dim = c(10, 10, 1))
set.seed(42)
kite_mn <- multinethist(kite_arr, h = 5, max_itr = 1000)

## print.nethist ----
test_that("print.nethist does not error", {
  expect_no_error(capture.output(print(hist_G)))
})

test_that("print.nethist returns object invisibly", {
  capture.output(result <- withVisible(print(hist_G)))
  expect_false(result$visible)
  expect_identical(result$value, hist_G)
})

test_that("print.nethist output mentions thetahat", {
  output <- capture.output(print(hist_G))
  expect_true(any(grepl("thetahat", output, ignore.case = TRUE)))
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

test_that("print.multinethist output mentions Theta_hat", {
  output <- capture.output(print(kite_mn))
  expect_true(any(grepl("Theta_hat|thetahat", output, ignore.case = TRUE)))
})

# Tests for the accessor functions and coef() methods.

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

set.seed(42)
fit_nh  <- suppressMessages(nethist(A, h = h_v))
set.seed(42)
fit_mnh <- suppressMessages(multinethist(A_multi, h = h_v))
set.seed(42)
fit_hnh <- suppressMessages(hnethist(A, h = h_v))

# -- block_id() ----------------------------------------------------------------

test_that("block_id() returns the integer block label of each vertex", {
  for (fit in list(fit_nh, fit_mnh, fit_hnh)) {
    ids <- block_id(fit)
    expect_type(ids, "integer")
    expect_length(ids, n_v)
    expect_equal(ids, as.integer(fit$cluster))
  }
})

test_that("block_id() rejects objects that are not network histogram fits", {
  expect_error(block_id(1:3), "no applicable method")
  expect_error(block_id(list(cluster = 1:3)), "no applicable method")
})

# -- bandwidth() ---------------------------------------------------------------

test_that("bandwidth() returns the bandwidth used for the fit", {
  for (fit in list(fit_nh, fit_mnh, fit_hnh)) {
    expect_equal(bandwidth(fit), fit$h)
  }
  expect_equal(bandwidth(fit_nh), h_v)
})

test_that("bandwidth() rejects objects that are not network histogram fits", {
  expect_error(bandwidth(1:3), "no applicable method")
  expect_error(bandwidth(list(h = 5)), "no applicable method")
})

# -- coef() --------------------------------------------------------------------

test_that("coef() on nethist returns the k x k block matrix", {
  k <- length(unique(fit_nh$cluster))
  expect_equal(dim(coef(fit_nh)), c(k, k))
  expect_equal(coef(fit_nh, type = "prob"), fit_nh$thetahat)
  expect_equal(coef(fit_nh), fit_nh$thetahat / fit_nh$rho_hat)
})

test_that("coef() default type matches fitted() default type", {
  expect_equal(coef(fit_nh)[fit_nh$cluster, fit_nh$cluster],
               fitted(fit_nh))
})

test_that("coef() on hnethist returns the merged block probabilities", {
  expect_equal(coef(fit_hnh, type = "prob"), fit_hnh$thetahat)
  expect_equal(coef(fit_hnh), fit_hnh$thetahat / fit_hnh$rho_hat)
})

test_that("coef() on multinethist returns a k x k x L array", {
  k  <- length(unique(fit_mnh$cluster))
  th <- coef(fit_mnh)
  expect_equal(dim(th), c(k, k, n_layers))
  expect_equal(coef(fit_mnh, type = "prob"), fit_mnh$thetahat)
  for (l in seq_len(n_layers)) {
    expect_equal(th[, , l], fit_mnh$thetahat[, , l] / fit_mnh$rho_hat[l])
  }
})

test_that("coef() handles fits that collapse to a single block", {
  fit_one <- suppressWarnings(suppressMessages(nethist(A, h = n_v)))
  expect_equal(dim(coef(fit_one)), c(1L, 1L))
  expect_equal(coef(fit_one), matrix(1))

  fit_one_m <- suppressWarnings(suppressMessages(
    multinethist(A_multi, h = n_v)))
  expect_equal(dim(coef(fit_one_m)), c(1L, 1L, n_layers))
  expect_equal(as.vector(coef(fit_one_m, type = "prob")), fit_one_m$rho_hat)
  expect_equal(as.vector(coef(fit_one_m)), rep(1, n_layers))
})

test_that("coef() rejects an unknown type", {
  expect_error(coef(fit_nh, type = "unknown"))
})

# -- shape_id() ----------------------------------------------------------------

test_that("shape_id() returns a symmetric k x k matrix of shape labels", {
  ids <- shape_id(fit_hnh)
  k   <- length(unique(fit_hnh$cluster))
  expect_equal(dim(ids), c(k, k))
  expect_type(ids, "integer")
  expect_true(isSymmetric(ids))
  expect_setequal(unique(as.vector(ids)), seq_len(fit_hnh$s))
})

test_that("shape_id() numbers shapes in increasing order of probability", {
  ids  <- shape_id(fit_hnh)
  prob <- coef(fit_hnh, type = "prob")
  keep <- upper.tri(ids, diag = TRUE)
  shape_prob <- tapply(prob[keep], ids[keep], mean)
  expect_false(is.unsorted(shape_prob))
})

test_that("blocks with the same shape share one probability", {
  ids  <- shape_id(fit_hnh)
  prob <- coef(fit_hnh, type = "prob")
  for (s in unique(as.vector(ids))) {
    expect_equal(length(unique(prob[ids == s])), 1L)
  }
})

test_that("shape_id() keeps the block clustering and only relabels it", {
  ids  <- shape_id(fit_hnh)
  keep <- upper.tri(ids, diag = TRUE)
  bc   <- fit_hnh$blockcluster
  # The relabeling is a bijection from the original labels to 1, ..., s.
  pairs <- unique(cbind(bc$cluster, ids[keep]))
  expect_equal(nrow(pairs), fit_hnh$s)
  expect_equal(length(unique(pairs[, 1])), fit_hnh$s)
  expect_equal(length(unique(pairs[, 2])), fit_hnh$s)
  expected <- match(bc$cluster, order(bc$centers[, 1]))
  expect_equal(ids[keep], expected)
})

test_that("shape_id() handles a merging step without kmeans", {
  # hnethist skips kmeans when every block keeps its own shape.
  fit_triv <- fit_hnh
  vec      <- fit_triv$initial$thetahat[upper.tri(fit_triv$initial$thetahat,
                                                  diag = TRUE)]
  fit_triv$blockcluster <- list(cluster = seq_along(vec),
                                centers = matrix(vec, ncol = 1))
  ids  <- shape_id(fit_triv)
  keep <- upper.tri(ids, diag = TRUE)
  expect_equal(ids[keep], match(seq_along(vec), order(vec)))
  expect_true(isSymmetric(ids))
})

test_that("shape_id() is only defined for hnethist fits", {
  expect_error(shape_id(fit_nh), "only defined for hnethist")
  expect_error(shape_id(fit_mnh), "only defined for hnethist")
})

# -- summary.hnethist() --------------------------------------------------------

test_that("summary() of hnethist reports blocks and probability by shape", {
  sm   <- summary(fit_hnh)
  ids  <- shape_id(fit_hnh)
  keep <- upper.tri(ids, diag = TRUE)
  expect_s3_class(sm$shapes, "data.frame")
  expect_equal(nrow(sm$shapes), fit_hnh$s)
  expect_equal(sum(sm$shapes$n_blocks), sum(keep))
  expect_false(is.unsorted(sm$shapes$probability))
  expect_output(print(sm), "Blocks and probability by shape")
})

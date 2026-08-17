##' Network histogram estimation
##'
##' Estimating network histogram for multiplex networks and returning the indices of partitions.
##'
##' @name multinethist
##' @param A Adjacency data for one or more network layers. Accepted formats:
##'   \itemize{
##'     \item \strong{Single-layer}: a \code{matrix}, \code{igraph} object, or
##'       \code{network} object.
##'     \item \strong{Multilayer}: a 3D \code{array} of dimension
##'       \eqn{n \times n \times L}; a \code{list} of \code{igraph} or
##'       \code{network} objects (all sharing a common vertex set); or a
##'       \code{combined_networks} object from \code{ergm.multi::Networks()}.
##'   }
##'   Plain \code{matrix} or sparse \code{Matrix} elements inside a list are
##'   not accepted because vertex ordering across layers cannot be verified.
##'   Use a 3D array instead.
##'
##'   When a list of \code{igraph} or \code{network} objects is supplied,
##'   vertex names are used to align layers if present (via
##'   \code{igraph::vertex_attr(g, "name")} or
##'   \code{network::network.vertex.names()}). If vertex names are absent,
##'   positional correspondence is assumed: vertex \eqn{i} in each layer is
##'   treated as the same vertex. All layers must have the same number of
##'   vertices, and when names are present they must form the same set.
##' @param h A bandwidth parameter. If `NA`, the bandwidth is selected by Olhede and Wolfe (2014). If specified, the user-supplied value is used.
##' @param common_f A logical; if `TRUE`, assumes a common network histogram function for all layers.
##' @param method Type of loss function for network histogram. Must be one of `PLL` (default) or `LSE`. `LSE` is implemented only for single-layer networks. See Details.
##' @param control A control object from \code{\link{nethist_control}}. Governs
##'   \code{max_itr}, \code{greedy_swap_rule}, \code{greedy_stop_threshold},
##'   and \code{verbose}.
##' @param ... Currently unused.
##' @returns
##' If the number of layers is greater than 1, it returns an object of class ``multinethist``:
##'
##' \itemize{
##' \item `cluster` a vector of partition indices.
##' \item `thetahat` a probability array from multinetwork histogram ordered by group labels.
##' \item `rho_hat` a vector of estimated sparsity parameters.
##' \item `normalized_LL` a normalized likelihood from the algorithm.
##' \item `homogeneous` a logical variable indicating homogeneous multinetwork histogram.
##' \item `h` bandwidth used for estimation.
##' }
##' @usage multinethist(A, h = NA, common_f = FALSE, method = "PLL",
##'   control = nethist_control(), ...)
##' @details {
##' The \eqn{l}th layer's multi-network histogram is defined by thetahat/rho_hat. The multinetwork histogram can be plotted using [plot()] and [plot3d()].
##'
##' If the number of layers is 1, multinethist() fits a single-layer network histogram, equivalent to nethist().
##'
##' `method` is only used for single-layer networks. `method = "PLL"` is for Olhede and Wolfe (2014), and `method = "LSE"` is for Gao et al. (2015).
##'
##' Note that `cluster` only shows a partition of vertices, and the index labels are not ordered. For example, vertices in cluster 1 do not have to be more similar to vertices in cluster 2 than to vertices in cluster 10. Hence, users may specify a custom order in [plot.multinethist()].
##' }
##' @examples
##' \donttest{
##'    #single-layer network histogram
##'    set.seed(42)
##'    data(polblog)
##'    nethist(polblog)
##'    nethist_polblog <- nethist(polblog)
##'    nethist_polblog_with_h <- nethist(polblog, h = 72)
##'
##'    #multi-network histogram
##'    set.seed(42)
##'    data(IndianVil)
##'    IndianVil
##'    mnethist_Ind_Vil <- multinethist(IndianVil)
##' }
##' @seealso [plot.multinethist()], [plot.nethist()], [nethist_control()]
##' @references Song, Y. & Olhede, S. C. (2026). Joint Estimation of Sparse Multilayer Networks via Graph Limits. https://arxiv.org/abs/2608.14536
##' @references Olhede, S. C. & Wolfe, P. J. (2014). Network Histograms and Universality of Blockmodel Approximation. Proceedings of the National Academy of Sciences, 111(41), 14722-14727. doi:10.1073/pnas.1400374111
##' @references Gao, C., Lu, Y., & Zhou, H. H. (2015). Rate-Optimal Graphon Estimation. The Annals of Statistics, 43(6), 2624-2652. doi:10.1214/15-AOS1354
##' @import Rcpp
##' @importFrom methods is
##' @importFrom stats .lm.fit dist pnorm weighted.mean
##' @importFrom RSpectra eigs
##' @export
multinethist <- function(A, h = NA, common_f = FALSE,
                         method = "PLL",
                         control = nethist_control(), ...) {
  UseMethod("multinethist")
}

##' @exportS3Method
multinethist.igraph <- function(A, h = NA, common_f = FALSE,
                                method = "PLL",
                                control = nethist_control(), ...) {
  multinethist.matrix(igraph::as_adjacency_matrix(A, sparse = FALSE),
                      h, common_f, method, control, ...)
}

##' @exportS3Method
multinethist.network <- function(A, h = NA, common_f = FALSE,
                                 method = "PLL",
                                 control = nethist_control(), ...) {
  multinethist.matrix(as.matrix(A, matrix.type = "adjacency"),
                      h, common_f, method, control, ...)
}

##' @exportS3Method
multinethist.combined_networks <- function(A, h = NA, common_f = FALSE,
                                           method = "PLL",
                                           control = nethist_control(), ...) {
  layer_id <- network::get.vertex.attribute(A, ".NetworkID")
  layers   <- sort(unique(layer_id))
  M        <- as.matrix(A, matrix.type = "adjacency")
  mats     <- lapply(layers, function(l) M[layer_id == l, layer_id == l])
  arr      <- .adj_list_to_array(mats)
  return(multinethist.array(arr, h, common_f, method, control, ...))
}

##' @exportS3Method
multinethist.list <- function(A, h = NA, common_f = FALSE,
                              method = "PLL",
                              control = nethist_control(), ...) {
  if (length(A) == 0L) stop("A must be a non-empty list.")
  mats <- lapply(seq_along(A), function(l) {
    el <- A[[l]]
    if (inherits(el, "combined_networks")) {
      stop(sprintf(
        "Layer %d is a combined_networks object. Unnest it first.", l))
    } else if (is.matrix(el) || (isS4(el) && methods::is(el, "Matrix"))) {
      stop(sprintf(paste0(
        "Layer %d is a plain matrix. Vertex ordering cannot be verified ",
        "across layers. Use a 3D array instead."), l))
    } else if (inherits(el, "igraph")) {
      igraph::as_adjacency_matrix(el, sparse = FALSE)
    } else if (inherits(el, "network")) {
      as.matrix(el, matrix.type = "adjacency")
    } else {
      stop(sprintf("Layer %d has unsupported class: %s.", l, class(el)[1L]))
    }
  })
  arr <- .adj_list_to_array(mats)
  return(multinethist.array(arr, h, common_f, method, control, ...))
}

# Helper: align a list of adjacency matrices by rownames and stack into array.
# Matrices with identical rownames (including positional "1","2",...) are
# accepted as-is; differing orders are permuted to match the first layer.
# Mismatched vertex sets raise an error.
.adj_list_to_array <- function(mats) {
  L    <- length(mats)
  dims <- vapply(mats, nrow, integer(1L))
  if (length(unique(dims)) > 1L)
    stop("All layers must have the same number of vertices.")
  n   <- dims[1L]
  ref <- rownames(mats[[1L]])
  for (l in seq_len(L)[-1L]) {
    cur <- rownames(mats[[l]])
    if (is.null(ref) || is.null(cur)) next   # no names: positional
    if (!setequal(cur, ref))
      stop(sprintf("Layer %d has different vertex names from layer 1.", l))
    if (!identical(cur, ref)) {
      perm       <- match(ref, cur)
      mats[[l]]  <- mats[[l]][perm, perm]
    }
  }
  arr <- array(0L, dim = c(n, n, L))
  for (l in seq_len(L)) arr[, , l] <- mats[[l]]
  return(arr)
}

##' @exportS3Method
multinethist.matrix <- function(A, h = NA, common_f = FALSE,
                                method = "PLL",
                                control = nethist_control(), ...) {
  multinethist.array(array(A, dim = c(nrow(A), ncol(A), 1)),
                     h, common_f, method, control, ...)
}

##' @exportS3Method
multinethist.array <- function(A, h = NA, common_f = FALSE,
                               method = "PLL",
                               control = nethist_control(), ...) {
  if (!is.array(A)) stop(paste0("A is not supported object:", class(A)))
  control <- .resolve_control(control, list(...))

  n_nodes  <- dim(A)[1]
  n_layers <- dim(A)[3]
  for (l in seq_len(n_layers)) {
    A_l <- A[, , l]
    if (any(diag(A_l) != 0))
      stop(paste0("Layer ", l, ": A has self-loops. Network A must be a simple graph (no self-loops)."))
    if (!isSymmetric(unname(A_l)))
      stop(paste0("Layer ", l, ": A is not symmetric. Network A must be an undirected graph."))
    if (any(A_l[A_l != 0] != 1))
      stop(paste0("Layer ", l, ": A is not a simple graph. All non-zero entries must be 1 (binary adjacency matrix)."))
  }
  method_char <- method
  method <- pmatch(method, c("PLL", "LSE")) # PLL = 1, LSE = 2
  if (is.na(method)) stop("method must be one of the followings: PLL, LSE")
  swap_rule_int <- pmatch(control$greedy_swap_rule, c("single_random"))
  if (is.na(swap_rule_int)) stop("greedy_swap_rule must be one of: single_random")

  # Compute necessary summaries from A
  rhoHat <- apply(A, 3, function(A_l) sum(A_l) / (n_nodes * (n_nodes - 1)))

  # Pick a bandwidth
  if (is.na(h)) {
    h <- .oracbwplugin(A, min(4, sqrt(n_nodes) / 8),
                       'degs', 1, rhoHat, common_f, control$verbose)$h
    if (control$verbose) message(paste("Determining bandwidth from data:", round(h)))
    h <- max(2, min(n_nodes, round(h)))
  } else {
    if (control$verbose) message(paste("Determining bandwidth from user input:", round(h)))
    h_clamped <- max(2, min(n_nodes, round(h)))
    if (h_clamped != h) {
      warning(paste0("h was adjusted to ", h_clamped,
                     " (original: ", h, ", valid range: [2, ", n_nodes, "])."))
    }
    h <- h_clamped
  }
  if (control$verbose) {
    message(paste("Final bandwidth:", h))
    message(paste0('Adjacency matrix has ', n_nodes, ' rows/cols'))
  }

  if (n_nodes < 2 * h) {
    warning(paste0("Bandwidth h=", h, " exceeds half the number of ",
                   "vertices n=", n_nodes, "; only one block is possible. ",
                   "Returning theta_hat = rho_hat."))
    if (n_layers == 1) {
      Log_Likelihood <- log(1 - rhoHat) * (n_nodes * (n_nodes - 1) / 2) +
        log(rhoHat / (1 - rhoHat)) * (sum(A) / 2)
      LSE <- sum((A - rhoHat)^2)
    } else {
      Log_Likelihood <- 0
      LSE <- 0
      for (l in 1:n_layers) {
        Log_Likelihood <- Log_Likelihood +
          log(1 - rhoHat[l]) * (n_nodes * (n_nodes - 1) / 2) +
          log(rhoHat[l] / (1 - rhoHat[l])) * (sum(A[, , l]) / 2)
        LSE <- LSE + sum((A[, , l] - rhoHat[l])^2)
      }
    }
    result <- list(cluster       = rep(1, n_nodes),
                   thetahat      = rhoHat,
                   rho_hat       = rhoHat,
                   normalized_LL = Log_Likelihood / (sum(A) / 2),
                   MSE           = LSE / n_nodes^2,
                   method        = method_char,
                   homogeneous   = common_f,
                   h             = h)
    result <- structure(result,
                        class = ifelse(n_layers > 1, "multinethist", "nethist"))
    return(result)
  }

  # Initialize using regularized spectral clustering (densest layer)
  tstart  <- Sys.time()
  idxInit <- initialize_index(A[, , which.max(rhoHat)], n_nodes, h,
                              control$verbose)
  k <- ceiling(n_nodes / h)

  if (control$verbose)
    message(paste0('Initial label vector assigned from row-similarity ordering; time ',
                   round(difftime(Sys.time(), tstart), 4), ' sec'))

  if (n_layers == 1) {
    res <- .nethist_fastgreedy(A[, , 1], h, Rind_to_Cind(idxInit), method,
                               control$max_itr, swap_rule_int,
                               control$greedy_stop_threshold,
                               control$verbose)
  } else if (common_f) {
    res <- .mnhistCommon_fastgreedy(A, h, Rind_to_Cind(idxInit),
                                    control$max_itr, swap_rule_int,
                                    control$greedy_stop_threshold,
                                    control$verbose)
  } else {
    res <- .multinethist_fastgreedy(A, h, Rind_to_Cind(idxInit),
                                    control$max_itr, swap_rule_int,
                                    control$greedy_stop_threshold,
                                    control$verbose)
  }

  result <- list(cluster       = as.vector(res$node_labels),
                 thetahat      = res$ThetaHat,
                 rho_hat       = rhoHat,
                 normalized_LL = res$norm_LL,
                 MSE           = res$LSE,
                 method        = method_char,
                 homogeneous   = common_f,
                 h             = h)
  result <- structure(result,
                      class = ifelse(n_layers > 1, "multinethist", "nethist"))
  return(result)
}

##' @exportS3Method
multinethist.default <- function(A, h = NA, common_f = FALSE,
                                 method = "PLL",
                                 control = nethist_control(), ...) {
  multinethist.array(A, h, common_f, method, control, ...)
}

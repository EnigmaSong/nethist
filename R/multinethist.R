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
##'   treated as the same node. All layers must have the same number of
##'   vertices, and when names are present they must form the same set.
##' @param h A bandwidth parameter. If `NA`, selecting bandwidth by Olhede and Wolfe (2014). If specified, the user input value is used.
##' @param common_f a logical variable indicating assume the common network histogram function for all layers.
##' @param method Type of loss functions for network histogram. It must be one of `PLL` (default) or `LSE`. `LSE` is implemented only for single-layer network histogram. See details
##' @param max_itr an integer for number of maximum iteration for greedy search.
##' @param swap_rule string of swap node selection methods. only "random" implemented. See details.
##' @param consecutive_iter_threshold an integer for stopping criterion. If the log-likelihood does not improve for the last `consecutive_iter_threshold` iterations, stop the algorithm.
##' @param verbose logical value indicating whether verbose output is generated.
##' @returns 
##' If number of layer is greater than 1, it returns an object of class ``multinethist``:
##' 
##' \itemize{
##' \item `cluster` a vector of partition indices.
##' \item `thetahat` a probability array from multinetwork histogram ordered by cluster labels. 
##' \item `rho_hat` a vector of estimated sparsity parameters. 
##' \item `normalized_LL` a normalized likelihood from the algorithm.
##' \item `homogeneous` a logical variable indicating homogeneous multinetwork histogram.
##' }
##' @usage nethist(A, h = NA, method = "PLL", max_itr = 5e6,
##'   swap_rule = "random", consecutive_iter_threshold = 2e4,
##'   verbose = FALSE)
##' @usage multinethist(A, h = NA, common_f = FALSE, method = "PLL",
##'   max_itr = 5e6, swap_rule = "random",
##'   consecutive_iter_threshold = 2e4, verbose = FALSE)
##' @details {
##' lth layer's multi-network histogram is defined by thetahat/rho_hat. We can plot multinetwork histogram using [plot()] and [plot3d()].
##' 
##' If the number of layer is 1, then it calls single-layer network histogram. See. nethist() in nethist package.
##' 
##' `method` is only used for single-layer networks. `method = "PLL"` is for Olhede and Wolfe (2014), and `method = "LSE"` is for Gao et al. (2015).
##' 
##' Note that `cluster` only shows a partition of vertices, and the index labels is not an ordered variable. For example, nodes in cluster 1 do not have to more similar to nodes in cluster 2 than nodes in cluster 10. Hence, users would use a user-specified order in [plot.multinethist()].
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
##' @seealso [plot.multinethist()], [plot.nethist()]
##' @references Song, Y. & Olhede, S. C. (2026+). Graph Limits for Sparse Multilayer Networks.
##' @references Olhede, S. C. & Wolfe, P. J. (2014). Network Histograms and Universality of Blockmodel Approximation. Proceedings of the National Academy of Sciences, 111(41), 14722-14727. doi:10.1073/pnas.1400374111
##' @references Gao, C., Lu, Y., & Zhou, H. H. (2015). Rate-Optimal Graphon Estimation. The Annals of Statistics, 43(6), 2624-2652. doi:10.1214/15-AOS1354
##' @import Rcpp
##' @importFrom stats .lm.fit dist pnorm weighted.mean
##' @importFrom RSpectra eigs
##' @export
multinethist <- function(A, h = NA, common_f = FALSE, 
                         method = "PLL",
                         max_itr = 5e6,
                         swap_rule = "random", 
                         consecutive_iter_threshold = 2e4,
                         verbose = FALSE){
  UseMethod("multinethist")
} 

##' @exportS3Method
multinethist.igraph <-  function(A, h = NA, common_f = FALSE, 
                                 method = "PLL",
                                 max_itr = 5e6,
                                 swap_rule = "random", 
                                 consecutive_iter_threshold = 2e4,
                                 verbose = FALSE){
  return(multinethist.matrix(igraph::as_adjacency_matrix(A, sparse = FALSE),
                             h, common_f, method,
                             max_itr, swap_rule, 
                             consecutive_iter_threshold, verbose)) 
}

##' @exportS3Method
multinethist.network <- function(A, h = NA, common_f = FALSE,
                                 method = "PLL",
                                 max_itr = 5e6,
                                 swap_rule = "random",
                                 consecutive_iter_threshold = 2e4,
                                 verbose = FALSE) {
  return(multinethist.matrix(as.matrix(A, matrix.type = "adjacency"),
                             h, common_f, method,
                             max_itr, swap_rule,
                             consecutive_iter_threshold, verbose))
}

##' @exportS3Method
multinethist.combined_networks <- function(A, h = NA, common_f = FALSE,
                                           method = "PLL",
                                           max_itr = 5e6,
                                           swap_rule = "random",
                                           consecutive_iter_threshold = 2e4,
                                           verbose = FALSE) {
  layer_id <- network::get.vertex.attribute(A, ".NetworkID")
  layers   <- sort(unique(layer_id))
  M        <- as.matrix(A, matrix.type = "adjacency")
  mats     <- lapply(layers, function(l) M[layer_id == l, layer_id == l])
  arr      <- .adj_list_to_array(mats)
  return(multinethist.array(arr, h, common_f, method, max_itr, swap_rule,
                            consecutive_iter_threshold, verbose))
}

##' @exportS3Method
multinethist.list <- function(A, h = NA, common_f = FALSE,
                              method = "PLL",
                              max_itr = 5e6,
                              swap_rule = "random",
                              consecutive_iter_threshold = 2e4,
                              verbose = FALSE) {
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
  return(multinethist.array(arr, h, common_f, method, max_itr, swap_rule,
                            consecutive_iter_threshold, verbose))
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
multinethist.matrix <-  function(A, h = NA, common_f = FALSE,
                                 method = "PLL",
                                 max_itr = 5e6,
                                 swap_rule = "random",
                                 consecutive_iter_threshold = 2e4,
                                 verbose = FALSE){
  return(multinethist.array(array(A, dim=c(nrow(A), ncol(A), 1)), 
                            h, common_f, method,
                            max_itr, swap_rule, 
                            consecutive_iter_threshold, verbose)) #should think about the design of code.
}

##' @exportS3Method
multinethist.array <- function(A, h = NA, common_f = FALSE, 
                               method = "PLL",
                         max_itr = 5e6,
                         swap_rule = "random", 
                         consecutive_iter_threshold = 2e4,
                         verbose = FALSE){
  if(!is.array(A)) stop(paste0("A is not supported object:", class(A)))
  n_nodes <- dim(A)[1]
  n_layers <- dim(A)[3]
  for(l in seq_len(n_layers)) {
    A_l <- A[, , l]
    if(any(diag(A_l) != 0))
      stop(paste0("Layer ", l, ": A has self-loops. Network A must be a simple graph (no self-loops)."))
    if(!isSymmetric(unname(A_l)))
      stop(paste0("Layer ", l, ": A is not symmetric. Network A must be an undirected graph."))
    if(any(A_l[A_l != 0] != 1))
      stop(paste0("Layer ", l, ": A is not a simple graph. All non-zero entries must be 1 (binary adjacency matrix)."))
  }
  method_char <- method
  method <- pmatch(method, c("PLL","LSE")) # PLL = 1, LSE = 2
  if(is.na(method)) stop("method must be one of the followings: PLL, LSE")
  swap_rule <- pmatch(swap_rule, c("random"))
  if(is.na(swap_rule)) stop("swap_rule must be one of the followings: random")
  
  # Compute necessary summaries from A
  rhoHat <- apply(A, 3, function(A_l) sum(A_l)/(n_nodes*(n_nodes-1)))
  
  # Pick a bandwidth
  if(is.na(h)){
    h <- .oracbwplugin(A, min(4, sqrt(n_nodes)/8),
                       'degs', 1, rhoHat, common_f, verbose)$h
    if(verbose) message(paste("Determining bandwidth from data:", round(h)))
    h <- max(2, min(n_nodes, round(h)))
  }else{
    if(verbose) message(paste("Determining bandwidth from user input:", round(h)))
    h_clamped <- max(2, min(n_nodes, round(h)))
    if (h_clamped != h) {
      warning(paste0("h was adjusted to ", h_clamped,
                     " (original: ", h, ", valid range: [2, ", n_nodes, "])."))
    }
    h <- h_clamped
  }
  if(verbose){
    message(paste("Final bandwidth:", h))
    message(paste0('Adjacency matrix has ', n_nodes, ' rows/cols'))
  }
  
  if(n_nodes == h){
    message(paste("Bandwidth h=", h, "and number of nodes=",n_nodes,"are equal. Return theta_hat = rho_n",sep =" "))
    if(n_layers == 1){
      Log_Likelihood <- log(1-rhoHat)*(n_nodes*(n_nodes-1)/2) + log(rhoHat/(1-rhoHat))*(sum(A)/2)
      LSE <- sum((A-rhoHat)^2)
    }else{
      Log_Likelihood <- 0
      LSE <- 0
      for(l in 1:n_layers){
        Log_Likelihood <- Log_Likelihood + log(1-rhoHat[l])*(n_nodes*(n_nodes-1)/2) + log(rhoHat[l]/(1-rhoHat[l]))*(sum(A[,,l])/2)
        LSE <- LSE + sum((A[,,l]-rhoHat[l])^2)
      }
    }
    result <- list(cluster = rep(1, n_nodes), 
                   thetahat =  rhoHat,
                   rho_hat = rhoHat,
                   normalized_LL = Log_Likelihood/(sum(A)/2),
                   MSE = LSE/n_nodes^2,
                   method = method_char,
                   homogeneous = common_f)
    result <- structure(result, class= ifelse(n_layers > 1, "multinethist", "nethist"))
    return(result)
  }
  # Initialize using regularized spectral clustering based on row similarity (used densest layer)
  tstart <- Sys.time()
  idxInit <- initialize_index(A[,,which.max(rhoHat)], n_nodes, h, verbose)
  k <- ceiling(n_nodes/h)
  
  if(verbose) message(paste0('Initial label vector assigned from row-similarity ordering; time ',
                             round(difftime(Sys.time(),tstart),4), ' sec'))
  
  if(n_layers == 1){
    res<- .nethist_fastgreedy(A[,,1], h, Rind_to_Cind(idxInit), method,
                              max_itr, swap_rule, consecutive_iter_threshold, verbose)
  }else if(common_f){
    res<- .mnhistCommon_fastgreedy(A, h, Rind_to_Cind(idxInit), 
                                   max_itr, swap_rule, consecutive_iter_threshold, verbose)
  }else{
    res <- .multinethist_fastgreedy(A, h, Rind_to_Cind(idxInit), 
                                    max_itr, swap_rule, consecutive_iter_threshold, verbose)
  }
  
  
  result <- list(cluster = as.vector(res$node_labels), 
                 thetahat =  res$ThetaHat,
                 rho_hat = rhoHat,
                 normalized_LL = res$norm_LL,
                 MSE = res$LSE,
                 method = method_char,
                 homogeneous = common_f)
  result <- structure(result, class= ifelse(n_layers > 1, "multinethist", "nethist"))
  return(result)
}

##' @exportS3Method
multinethist.default <- function(A, h = NA, common_f = FALSE, 
                                         method = "PLL",
                                         max_itr = 5e6,
                                         swap_rule = "random", 
                                         consecutive_iter_threshold = 2e4,
                                         verbose = FALSE){
  args <- as.list(environment())
  do.call(multinethist.array, args)
}
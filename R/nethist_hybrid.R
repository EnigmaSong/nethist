##' Hybrid Network histogram estimation
##'
##' Estimating hybrid network histogram for single-layer networks and returning the indices of partitions.
##'
##' @param A An adjacency matrix or an igraph object. It must be an undirected and simple graph.
##' @param h A bandwidth parameter. If `NA`, selecting bandwidth by Olhede and Wolfe (2014). If specified, the user input value is used.
##' @param method Type of loss functions for network histogram. It must be one of `LSE` (default) and `PLL` for single-layer hybrid network histogram. See details
##' @param max_itr an integer for number of maximum iteration for greedy search.
##' @param swap_rule string of swap node selection methods. only "random" implemented. See details.
##' @param consecutive_iter_threshold an integer for stopping criterion. If the log-likelihood does not improve for the last `consecutive_iter_threshold` iterations, stop the algorithm.
##' @param verbose logical value indicating whether verbose output is generated.
##' @returns 
##' If number of layer is greater than 1, it returns an object of class ``multinethist``:
##' 
##' \itemize{
##' \item `cluster` a vector of partition indices.
##' \item `thetahat` a probability matrix from hybrid network histogram ordered by cluster labels. 
##' \item `rho_hat` a vector of estimated sparsity parameters. 
##' \item `LSE` a LSE of the selected model.
##' \item `BIC` BIC value of the selected model
##' \item `details` list of models information: thetahat, LSE, BIC, and the number of shapes.
##' }
##' @usage hnethist(A, h = NA, method = "LSE", max_itr = 5e6, swap_rule = "random", consecutive_iter_threshold = 2e4, verbose = FALSE)
##' @details {
##' lth layer's multi-network histogram is defined by thetahat/rho_hat. We can plot multinetwork histogram using [plot()] and [plot3d()].
##' 
##' If the number of layer is 1, then it calls single-layer network histogram. See. nethist() in nethist package.
##' 
##' `method` is only used for single-layer networks. `method = "PLL"` is for Olhede and Wolfe (2014), and `method = "LSE"` is for Gao et al. (2015).
##' 
##' Note that `cluster` only shows a partition of vertices, and the index labels is not an ordered variable. For example, nodes in cluster 1 do not have to more similar to nodes in cluster 2 than nodes in cluster 10. Hence, users would use a user-specified order in [plot.multinethist()].
##' }
##' @seealso [plot.multinethist()], [plot.nethist()]
##' @references Song, Y. & Olhede. S. C. (2024)
##' @references Verdeyme, A. & Olhede. S. C. (2024). Hybrid of Node and Link Communities for Graphon Estimation. arXiv:2401.05088
##' @import Rcpp 
##' @importFrom stats kmeans
##' @export

hnethist <- function(A, h = NA, 
                         method = "LSE",
                         max_itr = 5e6,
                         swap_rule = "random", 
                         consecutive_iter_threshold = 2e4,
                         verbose = FALSE){
  UseMethod("hnethist")
} 

##' @exportS3Method
hnethist.igraph <- function(A, h = NA, 
                            method = "LSE",
                            max_itr = 5e6,
                            swap_rule = "random", 
                            consecutive_iter_threshold = 2e4,
                            verbose = FALSE){
  A <- igraph::as_adjacency_matrix(A, sparse = FALSE)
  return(hnethist.matrix(A, h, method, max_itr, swap_rule,
                         consecutive_iter_threshold, verbose))
}

##' @exportS3Method
hnethist.matrix <- function(A, h = NA, 
                                   method = "LSE",
                                   max_itr = 5e6,
                                   swap_rule = "random", 
                                   consecutive_iter_threshold = 2e4,
                                   verbose = FALSE){
  result <- list()
  
  result$initial <- nethist(A, h = h, method = method, max_itr = max_itr, 
          swap_rule = swap_rule, consecutive_iter_threshold = consecutive_iter_threshold,
          verbose = verbose)  
  n <- nrow(A)
  k <- length(unique(result$initial$cluster))
  s <- min(k*(k-1)/2.0, length(unique(c(result$initial$thetahat)))) #in case lots of duplicated values in thetahat
  result$details <- vector(mode = "list", length = s)
  
  #Model selection
  for(ss in s:1){
    result$details[[ss]]$n <- n
    result$details[[ss]]$s <- ss
    result$details[[ss]]$cluster <- hnethist_kmeans(result$initial, centers = ss)
    result$details[[ss]]$thetahat <- hnethist_smoothing(result$details[[ss]]$cluster, k)
    result$details[[ss]]$LSE <- hnethist_LSE(result$details[[ss]]$thetahat, result$initial$cluster, A)
    result$details[[ss]]$BIC <- hnethist_BIC(result$details[[ss]])
  }
  
  #Best model based on BIC
  BICs <- sapply(result$details, function(x) return(x$BIC))
  min_BIC_index <- which.min(BICs)
  
  result$cluster <- result$details[[min_BIC_index]]$cluster
  result$thetahat <- result$details[[min_BIC_index]]$thetahat
  
  result <- structure(result, class= "hnethist")
  
  return(result)
}
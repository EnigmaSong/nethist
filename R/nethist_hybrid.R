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
##' 
##' \itemize{
##' \item `initial` a nethist object from nethist(), which serves as the initial point of block clustering.
##' \item `blockcluster` A matrix of block partitions, where two blocks with the same index are part of the same cluster.
##' \item `thetahat` a probability matrix from hybrid network histogram ordered by cluster labels. 
##' \item `rho_hat` a vector of estimated sparsity parameters. 
##' \item `LSE` a LSE of the selected model.
##' \item `BIC` BIC value of the selected model.
##' \item `details` list of models information: kmeans clustering outputs, thetahat, LSE, BIC, and the number of shapes.
##' }
##' @usage hnethist(A, h = NA, method = "LSE", max_itr = 5e6, swap_rule = "random", 
##'                 consecutive_iter_threshold = 2e4, verbose = FALSE)
##' @details {
##' Among the outputs, the best model is selected based on the BIC criterion.
##' }
##' @references Verdeyme, A. & Olhede. S. C. (2024). Hybrid of Node and Link Communities for Graphon Estimation. arXiv:2401.05088
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
  
  result$blockcluster <- result$details[[min_BIC_index]]$cluster
  result$thetahat <- result$details[[min_BIC_index]]$thetahat
  result$LSE <- result$details[[min_BIC_index]]$LSE
  result$rho_hat <- result$initial$rho_hat
  
  result <- structure(result, class= "hnethist")
  
  return(result)
}
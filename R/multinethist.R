##' Network histogram estimation
##'
##' Estimating network histogram for multiplex networks and returning the indices of partitions.
##'
##' @name multinethist
##' @param A An adjacency array or list of igraph object. It must be an undirected and simple graph.
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
##' @usage nethist(A, h = NA, method = "PLL", max_itr = 5e6, swap_rule = "random", 
##' @usage         consecutive_iter_threshold = 2e4, verbose = FALSE)
##' @usage multinethist(A, h = NA, common_f = FALSE, method = "PLL", max_itr = 5e6, 
##' @usage              swap_rule = "random",consecutive_iter_threshold = 2e4, verbose = FALSE)
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
##.    nethist(polblog)
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
##' @references Song, Y. & Olhede. S. C. (2024)
##' @references Olhede, S. C. & Wolfe, P. J. (2014). Network Histograms and Universality of Blockmodel Approximation. Proceedings of the National Academy of Sciences, 111(41), 14722-14727. doi:10.1073/pnas.1400374111
##' @references #' Gao, C., Lu, Y., & Zhou, H. H. (2015). Rate-Optimal Graphon Estimation. The Annals of Statistics, 43(6), 2624-2652. doi:10.1214/15-AOS1354
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
  return(multinethist.matrix(igraph::as_adjacency_matrix(A), 
                             h, common_f, method,
                             max_itr, swap_rule, 
                             consecutive_iter_threshold, verbose)) 
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
  if(!is.array(A)) stop("A is not supported object")
  if(!all(apply(A, 3, .is_undirected_simple))) stop("Network A must be an undirected simple network.")
  method_char <- method
  method <- pmatch(method, c("PLL","LSE")) # PLL = 1, LSE = 2
  if(is.na(method)) stop("method must be one of the followings: PLL, LSE")
  n_nodes <- dim(A)[1]
  n_layers <- dim(A)[3]
  swap_rule <- pmatch(swap_rule, c("random"))
  if(is.na(swap_rule)) stop("swap_rule must be one of the followings: random")
  
  # Compute necessary summaries from A
  rhoHat <- apply(A, 3, function(A_l) sum(A_l)/(n_nodes*(n_nodes-1)))
  
  # Pick a bandwidth
  if(is.na(h)){
    h <- .oracbwplugin(A, min(4, sqrt(n_nodes)/8), 
                       'degs', 1, rhoHat, common_f, verbose)$h
    if(verbose) message(paste("Determining bandwidth from data:", round(h)))
  }else{
    if(verbose) message(paste("Determining bandwidth from user input:", round(h)))
  }
  
  h <- max(2, min(n_nodes, round(h)))
  if(verbose){
    message(paste("Final bandwidth:", h))
    message(paste0('Adjacency matrix has ', n_nodes, ' rows/cols'))
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
                 LSE = res$LSE,
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
  args <- ls()
  do.call(multinethist.array, args)
}
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
##' @returns A list of class `c("hnethist", "nethist")` with the following fields:
##'
##' \itemize{
##' \item `cluster` an integer vector of length n with node-level block assignments from the initial nethist fit.
##' \item `thetahat` a k-by-k probability matrix of the selected hybrid model, where k is the number of nethist blocks.
##' \item `rho_hat` estimated sparsity parameter from the initial nethist fit.
##' \item `normalized_LL` normalized log-likelihood from the initial nethist fit.
##' \item `LSE` least squared error of the selected model.
##' \item `method` loss function used (`"LSE"` or `"PLL"`).
##' \item `blockcluster` a `kmeans` object describing how the k-by-k blocks were merged into `s` shapes.
##' \item `BIC` BIC value of the selected model.
##' \item `s` number of distinct shapes in the selected model.
##' \item `details` list of length `s_max`, one entry per candidate model, each containing `s`, `cluster`, `thetahat`, `LSE`, and `BIC`.
##' \item `initial` the nethist object used as the starting point for block clustering.
##' }
##' @usage hnethist(A, h = NA, method = "LSE", max_itr = 5e6, swap_rule = "random",
##'                 consecutive_iter_threshold = 2e4, verbose = FALSE)
##' @details {
##' Among the outputs, the best model is selected based on the BIC criterion.
##' 
##' The original reference provided theoretical guarantees for LSE, but we also allow PLL for the initial nethist fit. The block clustering and model selection steps are performed on the LSE loss regardless of the initial method, as the theoretical results pertain to LSE.
##' }
##' @references Verdeyme, A. & Olhede. S. C. (2024). Hybrid of Node and Link Communities for Graphon Estimation. arXiv:2401.05088
##' @rdname hnethist
##' @export
hnethist <- function(A, h = NA,
                         method = "LSE",
                         max_itr = 5e6,
                         swap_rule = "random",
                         consecutive_iter_threshold = 2e4,
                         verbose = FALSE){
  if (inherits(A, "igraph")) {
    A <- igraph::as_adjacency_matrix(A, sparse = FALSE)
  } else if (!is.matrix(A)) {
    stop("A must be a matrix or igraph object.")
  }
  return(.hnethist_matrix(A, h, method, max_itr, swap_rule,
                          consecutive_iter_threshold, verbose))
}

.hnethist_matrix <- function(A, h = NA,
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

  # Select the best model by minimum BIC
  BICs <- sapply(result$details, function(x) x$BIC)
  best <- result$details[[which.min(BICs)]]

  return(structure(
    list(
      # nethist-compatible fields (enables S3 method inheritance)
      cluster       = result$initial$cluster,
      thetahat      = best$thetahat,
      rho_hat       = result$initial$rho_hat,
      normalized_LL = result$initial$normalized_LL,
      LSE           = best$LSE,
      method        = result$initial$method,
      # hnethist-specific fields
      blockcluster  = best$cluster,
      BIC           = best$BIC,
      s             = best$s,
      details       = result$details,
      initial       = result$initial
    ),
    class = c("hnethist", "nethist")
  ))
}

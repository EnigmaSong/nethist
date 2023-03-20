##' Network histogram estimation
##'
##' Estimating network histogram and returning the indices of partitions.
##'
##' @param A An adjacency matrix or igraph object. It must be an undirected and simple graph.
##' @param h A bandwidth parameter. If `NA`, selecting bandwidth by Olhede and Wolfe (2014). If specified, the user input value is used.
##' @param outfile A filename for saving cluster indices. If it is missing, the results are not saved in a file.
##' @param verbose logical value indicating whether verbose output is generated.
##' @returns 
##' An object of class ``nethist``:
##' 
##' \itemize{
##' \item `cluster` a vector of partition indices.
##' \item `p_mat` a probability matrix from network histogram ordered by cluster labels. 
##' \item `rho_hat` estimated sparsity parameter. 
##' }
##' @details 
##' Note that `cluster` only shows a partition of vertices, and the index labels is not an ordered variable. For example, nodes in cluster 1 do not have to more similar to nodes in cluster 2 than nodes in cluster 10. Hence, users would use a user-specified order in [plot.nethist()].
##' 
##' Network histogram values are computed by `p_mat/rho_hat` in `nethist` object.
##' 
##' ``nethist()`` runs the following two steps:
##' \itemize{
##' \item Initialization: With spectral clustering, find an initial partitions of vertices. If bandwidth `h` is not specified, automatic bandwidth selection is used.
##' \item Using a greedy search algorithm, find a (local) optima of normalized profile log-likelihood, which is `cluster` in the ``nethist`` object.
##' }
##' @seealso [plot.nethist()]
##' @references Olhede, S. C., & Wolfe, P. J. (2014). Network histograms and universality of blockmodel approximation. Proceedings of the National Academy of Sciences, 111(41), 14722-14727.
##' @references Wolfe, P. J., & Olhede, S. C. (2013). Nonparametric graphon estimation. arXiv preprint arXiv:1309.5936.
##' @examples
##' \dontrun{
##' set.seed(2022)
##' #Generating Erdos-Renyi graph
##' A <- igraph::sample_gnp(100, 0.05)
##' 
##' #With automatic bandwidth selection
##' hist_A <- nethist(A) 
##' 
##' #with user-specified bandwidth
##' hist_A <- nethist(A, h = 20) 
##' 
##' #with adjancency matrix
##' hist_A <- nethist(igraph::as_adj(A)) 
##' }
##' @importFrom stats .lm.fit dist pnorm
##' @importFrom graphics par
##' @importFrom utils write.table
##' @importFrom RSpectra eigs
##' @export
nethist<-function(A, h = NA, outfile, verbose = FALSE){
  UseMethod("nethist")
}
##' @exportS3Method 
nethist.igraph<-function(A, h, outfile, verbose){
  args <- as.list(environment())
  args$A <- igraph::as_adj(args$A, sparse = FALSE)
  do.call("nethist.default", args)
}
##' @exportS3Method 
nethist.matrix<-function(A, h, outfile, verbose){
  args <- as.list(environment())
  do.call("nethist.default", args)
}
##' @exportS3Method 
nethist.dgCMatrix<-function(A, h, outfile, verbose){
  args <- as.list(environment())
  args$A <- as.matrix(args$A)
  do.call("nethist.default", args)
}
##' 
nethist.default <- function(A, h = NA, outfile, verbose = F){
  check_input_error(A, h, verbose)
  
  n <- dim(A)[1L]
  rhoHat <- sum(A)/(n*(n-1))
  if(verbose) message(paste0('Adjacency matrix has ', n, ' rows/cols'))
  
  h <- get_bandwidth(A, n, rhoHat, h, verbose)
  
  # Initialize using regularized spectral clustering based on row similarity
  tstart <- Sys.time()
  idx <- initialize_index(A, n, h, verbose)
  if(verbose) message(paste0('Initial label vector assigned from row-similarity ordering; time ',
                 round(difftime(Sys.time(),tstart),4), ' sec'))
  
  idx <- .graphest_fastgreedy(A, h, idx, verbose)
  
  if(!missing(outfile)){
    write.table(file=outfile, x = idx, row.names = FALSE, col.names = FALSE)
  }
  
  result <- list(cluster = as.vector(idx), 
                 p_mat = .prob_mat_from_adj(A,idx),
                 rho_hat = rhoHat)
  result <- structure(result, class="nethist")
  return(result)
}


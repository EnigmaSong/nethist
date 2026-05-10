##' @rdname multinethist
##' @export
nethist<-function(A, h = NA, 
                  method = "PLL",
                  max_itr = 5e6,
                  swap_rule = "random", 
                  consecutive_iter_threshold = 2e4,
                  verbose = FALSE){
  UseMethod("nethist")
}

##' @exportS3Method 
nethist.igraph<-function(A, h = NA, 
                         method = "PLL",
                         max_itr = 5e6,
                         swap_rule = "random", 
                         consecutive_iter_threshold = 2e4,
                         verbose = FALSE){
  args <- as.list(environment())
  args$A <- igraph::as_adjacency_matrix(args$A, sparse = FALSE)
  do.call("nethist.default", args)
}
##' @exportS3Method
nethist.matrix<-function(A, h = NA,
                         method = "PLL",
                         max_itr = 5e6,
                         swap_rule = "random",
                         consecutive_iter_threshold = 2e4,
                         verbose = FALSE){
  # Explicit dispatch for matrix class: passes through to nethist.default without conversion.
  # Required so that S3 dispatch does not fall through to an unintended method when A is a matrix.
  args <- as.list(environment())
  do.call("nethist.default", args)
}
##' @exportS3Method
nethist.dgCMatrix<-function(A, h = NA,
                            method = "PLL",
                            max_itr = 5e6,
                            swap_rule = "random",
                            consecutive_iter_threshold = 2e4,
                            verbose = FALSE){
  args <- as.list(environment())
  args$A <- as.matrix(args$A)
  do.call("nethist.default", args)
}
##' @exportS3Method
nethist.combined_networks <- function(A, h = NA,
                                      method = "PLL",
                                      max_itr = 5e6,
                                      swap_rule = "random",
                                      consecutive_iter_threshold = 2e4,
                                      verbose = FALSE) {
  stop("A is a combined_networks (multilayer) object. Use multinethist() instead.")
}
##' @exportS3Method
nethist.network <- function(A, h = NA,
                            method = "PLL",
                            max_itr = 5e6,
                            swap_rule = "random",
                            consecutive_iter_threshold = 2e4,
                            verbose = FALSE) {
  args <- as.list(environment())
  args$A <- as.matrix(args$A, matrix.type = "adjacency")
  do.call("nethist.default", args)
}
##' @exportS3Method
nethist.default <- function(A, h = NA, 
                            method = "PLL",
                            max_itr = 5e6,
                            swap_rule = "random", 
                            consecutive_iter_threshold = 2e4,
                            verbose = FALSE){
  if(!is.matrix(A)) stop("A must be a matrix.")
  if(nrow(A) != ncol(A)) stop("A must be a square matrix.")
  if(any(diag(A) != 0)) stop("A has self-loops. Network A must be a simple graph (no self-loops).")
  if(!isSymmetric(unname(A))) stop("A is not symmetric. Network A must be an undirected graph.")
  if(any(A[A != 0] != 1)) stop("A is not a simple graph. All non-zero entries must be 1 (binary adjacency matrix).")
  multinethist.array(array(A, dim=c(nrow(A), ncol(A), 1)), 
                     h, common_f = FALSE, method,
                     max_itr, swap_rule, 
                     consecutive_iter_threshold, verbose)
}


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
nethist.default <- function(A, h = NA, 
                            method = "PLL",
                            max_itr = 5e6,
                            swap_rule = "random", 
                            consecutive_iter_threshold = 2e4,
                            verbose = FALSE){
  if(!is.matrix(A)) stop("A must be a matrix.")
  if(!.is_undirected_simple(A)) stop("Network A must be an undirected simple single-layer network.")
  multinethist.array(array(A, dim=c(nrow(A), ncol(A), 1)), 
                     h, common_f = FALSE, method,
                     max_itr, swap_rule, 
                     consecutive_iter_threshold, verbose)
}


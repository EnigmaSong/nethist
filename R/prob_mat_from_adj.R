##' Probability matrix computation from network histogram cluster indices 
##'
##' From an adjacency matrix and a cluster index vector, it returns a probability matrix
##'
##' @param A A symmetric adjacency matrix. It must have binary entries (0 or 1).
##' @param idx an index vector consists of integers from 1 to K
##' @returns 
##' A matrix whose entries are in between 0 and 1.
##' 
##' @examples
##' set.seed(2022)
##' #Generating Erdos-Renyi graph
##' A <- igraph::sample_gnp(100, 0.05)
##' hist_A <- nethist(A) 
##' 
##' all.equal(.prob_mat_from_adj(igraph::as_adj(A, sparse= FALSE), hist_A$cluster), hist_A$p_mat)
##' @keywords internal
##' @export
.prob_mat_from_adj<-function(A, idx){
  K <- max(idx)  
  p_mat <- matrix(0,K,K)
  
  numer <- rowsum(t(rowsum(A, idx)), idx)
  
  bin_size <- table(idx)
  denom <- outer(bin_size,bin_size)
  diag(denom) <- diag(denom) - bin_size
  p_mat<- numer/denom
  dimnames(p_mat)<- list(1:K, 1:K)
  return(p_mat)
}



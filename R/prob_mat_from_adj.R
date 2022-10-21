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
##' \dontrun{
##' set.seed(2022)
##' #Generating Erdos-Renyi graph
##' A <- igraph::sample_gnp(100, 0.05)
##' A <- igraph::as_adj(A)
##' res <- nethist(A) 
##' 
##' prob_mat_from_adj(A, res$cluster)
##' all.equal(prob_mat_from_adj(A, res$cluster), res$p_mat)
##' }
##' @keywords internal
##' @export
prob_mat_from_adj<-function(A, idx){
  K <- max(idx)  
  p_mat <- matrix(0,K,K)
  
  numer <- rowsum(t(rowsum(A, idx)), idx)
  
  bin_size <- table(idx)
  denom <- outer(bin_size,bin_size)
  diag(denom) <- diag(denom) - bin_size
  p_mat<- numer/denom
  dimnames(p_mat)<-NULL
  return(p_mat)
}



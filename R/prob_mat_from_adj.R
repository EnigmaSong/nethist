##' Probability matrix computation from cluster indices
##'
##' From an adjancency matrix and a cluster index vector, it returns probability matrix
##'
##' @param A A symmetric adjacency matrix. It must have binary entries (0 or 1).
##' @param idx an index vector consists of integers from 1 to K
##' @returns 
##' A matrix whose entires in between 0 and 1.
##' 
##' @examples
##' \dontrun{
##' set.seed(2022)
##' #Generating Erdos-Renyi graph
##' A <- igraph::sample_gnp(100, 0.05)
##' A <- igraph::as_adj(A)
##' res <- nethist(A) #Save the result in idx, do not save it in a csv file.
##' 
##' prob_mat_from_adj(A, res$idx)
##' all.equal(prob_mat_from_adj(A, res$idx), res$p_mat)
##' }
##' @export
prob_mat_from_adj<-function(A, idx){
  K <- max(idx)  
  p_mat <- matrix(0,K,K)
  for(i in 1:K){
    for(j in i:K){
      adj.mat.block = A[idx==i,idx==j]
      if(i != j){
        p_mat[i,j] = sum(adj.mat.block)/prod(dim(adj.mat.block))
        p_mat[j,i] = p_mat[i,j]
      }else{
        dim_block = dim(adj.mat.block)[1]
        p_mat[i,j] = sum(adj.mat.block)/max(1,dim_block*(dim_block-1))
      }
    }
  }
  return(p_mat)
}
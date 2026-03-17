##' @noRd
hnethist_kmeans <- function(nh, centers){
  vec_thetahat <- c(nh$thetahat)
  cluster <- kmeans(vec_thetahat, centers = centers)

  return(cluster)
}
#kmeans++ later?

hnethist_smoothing <- function(nh, k){
  thetahat <- matrix(nh$centers[nh$cluster], k, k)
  return(thetahat)
}
#LSE 
hnethist_LSE <- function(thetahat, cluster_initial, A){
  LSE <- sum((A - thetahat[cluster_initial,cluster_initial])^2)
  return(LSE)
}

## BIC
hnethist_BIC <- function(nh){
  #Need to check the implementation of BIC in Arthur's code.
  return(- 2*log(nh$LSE) + nh$n*log(nh$s))
}
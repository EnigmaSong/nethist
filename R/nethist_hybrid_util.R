##' @importFrom stats kmeans
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
## Gaussian BIC: N = n(n-1)/2 unique pairs, RSS = LSE/2 (hnethist_LSE sums
## both (i,j) and (j,i)), penalty = s * log(N).
hnethist_BIC <- function(nh){
  N <- nh$n * (nh$n - 1) / 2
  return(N * log(nh$LSE / (2 * N)) + nh$s * log(N))
}
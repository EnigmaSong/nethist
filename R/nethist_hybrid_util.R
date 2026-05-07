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
hnethist_LSE <- function(thetahat, cluster_initial, A){
  n <- nrow(A)
  LSE <- sum((A - thetahat[cluster_initial, cluster_initial])^2)
  return(LSE / n^2)
}

hnethist_normalized_LL <- function(thetahat, cluster_initial, A) {
  theta_mat <- thetahat[cluster_initial, cluster_initial]
  eps <- .Machine$double.eps
  theta_mat <- pmax(pmin(theta_mat, 1 - eps), eps)
  LL <- sum(theta_mat * log(theta_mat) + (1 - theta_mat) * log(1 - theta_mat))
  return(LL / sum(A))
}

## BIC
## Gaussian BIC: N = n(n-1)/2 unique pairs, RSS = MSE * n^2 / 2 (hnethist_LSE
## sums both (i,j) and (j,i) then divides by n^2), penalty = s * log(N).
hnethist_BIC <- function(nh){
  N <- nh$n * (nh$n - 1) / 2
  return(N * log(nh$MSE * nh$n^2 / (2 * N)) + nh$s * log(N))
}
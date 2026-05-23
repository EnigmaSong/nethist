##' @importFrom stats kmeans
##' @noRd
hnethist_kmeans <- function(nh, centers){
  vec_thetahat <- nh$thetahat[upper.tri(nh$thetahat, diag = TRUE)]
  # kmeans requires centers < length(vec): handle trivial case separately
  if (centers >= length(vec_thetahat)) {
    return(list(
      cluster = seq_along(vec_thetahat),
      centers = matrix(vec_thetahat, ncol = 1)
    ))
  }
  cluster <- kmeans(vec_thetahat, centers = centers)
  return(cluster)
}
#kmeans++ later?

hnethist_smoothing <- function(nh, k){
  mat <- matrix(0, k, k)
  mat[upper.tri(mat, diag = TRUE)] <- nh$centers[nh$cluster]
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
  return(mat)
}
hnethist_LSE <- function(thetahat, cluster_initial, A){
  n <- nrow(A)
  LSE <- sum((A - thetahat[cluster_initial, cluster_initial])^2)
  return(LSE / n^2)
}

hnethist_LL <- function(thetahat, cluster_initial, A) {
  theta_mat <- thetahat[cluster_initial, cluster_initial]
  eps <- .Machine$double.eps
  theta_mat <- pmax(pmin(theta_mat, 1 - eps), eps)
  upper <- upper.tri(A)
  return(sum(A[upper] * log(theta_mat[upper]) +
             (1 - A[upper]) * log(1 - theta_mat[upper])))
}

hnethist_normalized_LL <- function(thetahat, cluster_initial, A) {
  return(hnethist_LL(thetahat, cluster_initial, A) / sum(A))
}

## BIC
## Bernoulli BIC: N = n(n-1)/2 unique pairs, penalty = s * log(N).
hnethist_BIC <- function(nh, ll){
  N <- nh$n * (nh$n - 1) / 2
  return(-2 * ll + nh$s * log(N))
}
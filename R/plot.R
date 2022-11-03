##' Network histogram plot
##'
##' Drawing [heatmap()] using an `nethist` object with an user-specified order.
##'
##' @param x a nethist object from [nethist()].
##' @param idx_order A numeric vector for index label order, which must be a permutation of x$cluster.
##' @param ... other arguments to pass to [stats::heatmap()].
##' @returns 
##' a heatmap 
##' @examples
##' \dontrun{
##' set.seed(2022)
##' #Generating Erdos-Renyi graph
##' A <- igraph::sample_gnp(200, 0.05)
##' hist_A <- nethist(A)
##' plot(hist_A)
##' }
##' @importFrom stats heatmap
##' @exportS3Method 
##' @export
plot.nethist <- function(x, idx_order = 1:max(x$cluster), ...){
  if(!.is_valid_order(idx_order, 1:max(x$cluster))){
    warning(paste0("idx_order is invalid. Set idx_order = 1:",max(x$cluster)))
    idx_order <- 1:max(x$cluster)
  }
  heatmap(x$p_mat[idx_order, idx_order], Rowv = NA, symm = TRUE, ...)
}
##' Network histogram plot
##'
##' Estimating network histogram and returning the indices of partitions.
##'
##' @param x a nethist object from [nethist()].
##' @param idx_order index label order specifies 
##' @param ... other arguements to pass to [stats::heatmap()].
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
  if(!is_valid_order(idx_order)) 
  heatmap(x$p_mat[idx_order, idx_order], Rowv = NA, symm = TRUE, ...)
}
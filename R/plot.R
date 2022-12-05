##' Network histogram plot
##'
##' Drawing [heatmap()] using an `nethist` object with an user-specified order.
##'
##' @param x a nethist object from [nethist()].
##' @param idx_order A numeric vector for index label order, which must be a permutation of `x$cluster`. If `NA`, it uses `1:max(x$clsuter)`. 
##' @param ... other arguments to pass to [stats::heatmap()]. See details.
##' @details 
##' ... includes various graphical parameters passes to [stats::heatmap()], then [graphics::image()]. 
##' @returns 
##' a heatmap of `p_mat` orderd by `idx_order` in ``nethist`` object.
##' @examples
##' \dontrun{
##' set.seed(2022)
##' #Generating Erdos-Renyi graph
##' A <- igraph::sample_gnp(200, 0.05)
##' hist_A <- nethist(A)
##' plot(hist_A)
##' 
##' #with user-specified order
##' idx<- unique(hist_A$cluster) 
##' plot(hist_A, idx_order = idx)
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
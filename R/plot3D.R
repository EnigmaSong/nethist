##' Network histogram 3D plot
##'
##' Drawing [heatmap()] using an `nethist` object with an user-specified order.
##'
##' @param x a nethist object from [nethist()].
##' @param idx_order A numeric vector for index label order, which must be a permutation of `x$cluster`. If `NA`, it uses `1:max(x$clsuter)`. 
##' @param ... other arguments to pass to [plot3D::hist3D()]. See details.
##' @details 
##' ... includes various [`graphical parameters`] passes to [plot3D::hist3D()], then [graphics::image()]. 
##' @returns 
##' a heatmap of network histogram or `p_mat` ordered by `idx_order` from ``nethist`` object.
##' @examples
##' \dontrun{
##' set.seed(42)
##' hist_polblog <- nethist(polblog)
##' plot3d(hist_polblog)
##' 
##' #with user-specified order
##' idx<- c(17,10,13,7,12,4,8,5,1,2,6,3,9,15,11,14,16)
##' plot3d(hist_polblog, idx_order = idx)
##' 
##' #Rotate 3D plot
##' plot3d(hist_polblog, idx_order = idx, phi = 20, theta = 60)
##' }
##' @importFrom plot3D hist3D
##' @export
plot3d <- function(x, idx_order = 1:max(x$cluster), 
                   ...){
  UseMethod("plot3d")
}
##' @exportS3Method 
plot3d.nethist <- function(x, 
                         idx_order = 1:max(x$cluster), 
                         ...){
  k<-max(x$cluster)

  if(!.is_valid_order(idx_order, 1:k)){
    warning(paste0("idx_order is invalid. Set idx_order = 1:",k))
    idx_order <- 1:k
  }
  
  mat <- x$p_mat[idx_order, idx_order]/x$rho_hat
  
  hist3D(z=mat, ...)
}
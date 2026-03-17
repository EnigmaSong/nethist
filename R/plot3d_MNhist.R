##' Multinetwork histogram 3D plot
##'
##' Drawing [heatmap()] using an `multinethist` object with an user-specified order.
##'
##' @param x a nethist object from [multinethist()].
##' @param idx_order A numeric vector for index label order, which must be a permutation of `x$cluster`. If `NA`, it uses `1:max(x$cluster)`. 
##' @param type One of `MNhist` or `prob`.
##' @param ... other arguments to pass to [plot3D::hist3D()]. See details.
##' @details 
##' ... includes various [`graphical parameters`] passes to [plot3D::hist3D()], then [graphics::image()]. 
##' @returns
##' Called for its side effects (plotting). Returns `NULL` invisibly.
##' @examples
##' \donttest{
##' set.seed(42)
##' data(IndianVil)
##' mnhist_Ind_vil <- multinethist(IndianVil)
##' plot3d(mnhist_Ind_vil)
##' }
##' @importFrom plot3D hist3D
##' @rdname plot3d.multinethist
##' @export
plot3d <- function(x, idx_order = 1:max(x$cluster), 
                   type = "MNhist",
                   ...){
  UseMethod("plot3d")
}
##' @exportS3Method 
plot3d.multinethist <- function(x, 
                           idx_order = 1:max(x$cluster), 
                           type = "MNhist",
                           ...){
  k<-max(x$cluster)
  
  if(!.is_valid_order(idx_order, 1:k)){
    warning(paste0("idx_order is invalid. Set idx_order = 1:",k))
    idx_order <- 1:k
  }
  
  n_loop <- ifelse((x$homogeneous)&(type=="MNhist"), 1, dim(x$thetahat)[3])
  for(l in 1:n_loop){
    mat <- x$thetahat[idx_order, idx_order, l]/x$rho_hat[l]

    hist3D(z=mat, main = paste("Layer", l),...)
  }
  invisible(NULL)
}
##' Network histogram plot
##'
##' Drawing [heatmap()] using an `nethist` object with an user-specified order.
##'
##' @param x a nethist object from [nethist()].
##' @param idx_order A numeric vector for index label order, which must be a permutation of `x$cluster`. If `NA`, it uses `1:max(x$clsuter)`. 
##' @param prob A logical variable indicating block probabilities are printed on the plot. Default is FALSE.
##' @param digits integer indicating the number of decimal places for probability
##' @param prob.cex A numerical value of `cex` of probabilites. 
##' @param prob.col color of probabilities on each bin.
##' @param ... other arguments to pass to [stats::heatmap()]. See details.
##' @details 
##' ... includes various [`graphical parameters`] passes to [stats::heatmap()], then [graphics::image()]. 
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
##' 
##' #User-speicifc bin color pallete (see [graphical parameters])
##' plot(hist_polblog,  idx_order = ind, col = colorRampPalette(colors=c("#FFFFFF","#000000"))(50))
##' 
##' #Users can print p_mat on the plot using user-specific colors
##' plot(hist_polblog,  idx_order = ind, prob= TRUE, prob.col = "blue",
##'       col = colorRampPalette(colors=c("#FFFFFF","#000000"))(50))
##' }
##' @importFrom stats heatmap
##' @exportS3Method 
##' @export
plot.nethist <- function(x, idx_order = 1:max(x$cluster), 
                         prob = FALSE, digits = 2,
                         prob.cex =  0.1 + 0.5/log10(max(x$cluster)),
                         prob.col = "black",
                         ...){
  k<-max(x$cluster)
  if(!.is_valid_order(idx_order, 1:k)){
    warning(paste0("idx_order is invalid. Set idx_order = 1:",k))
    idx_order <- 1:k
  }
  
  if(prob){
    heatmap(x$p_mat[idx_order, idx_order], Rowv = NA, symm = TRUE, 
            add.expr = {text(rep(1:k,each=k), rev(rep(1:k,k)), 
                             round(as.vector(x$p_mat[idx_order, idx_order]), digits),
                             cex = prob.cex, col = prob.col)},
            ...)
  }else{
    heatmap(x$p_mat[idx_order, idx_order], Rowv = NA, symm = TRUE, ...)
  }
}
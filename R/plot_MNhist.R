##' Network histogram plot
##'
##' Drawing [heatmap()] using an `multinethist` object with an user-specified order.
##'
##' @param x a multinethist object from [multinethist()].
##' @param y a dummy variable. It does not affect to the plot.
##' @param idx_order A numeric vector for index label order, which must be a permutation of `x$cluster`. If `NA`, it uses `1:max(x$cluster)`. 
##' @param type One of `MNhist` or `prob`.
##' @param prob A logical variable indicating block probabilities are printed on the plot. Default is FALSE.
##' @param digits integer indicating the number of decimal places for probability
##' @param prob.cex A numerical value of `cex` of probabilities. 
##' @param prob.col A color vector for network histogram values/probabilities on each bin.
##' @param ... other arguments to pass to [stats::heatmap()]. See details.
##' @details 
##' ... includes various [`graphical parameters`] passes to [stats::heatmap()], then [graphics::image()]. 
##' @returns 
##' a heatmap of network histogram or `p_mat` ordered by `idx_order` from ``nethist`` object.
##' @examples
##' \donttest{
##' set.seed(42)
##' data(IndianVil)
##' mnhist_Ind_vil <- multinethist(IndianVil)
##' 
##' #with user-specified order
##' idx<- unique(mnhist_Ind_vil$cluster) 
##' plot(mnhist_Ind_vil, idx_order = idx)
##' 
##' #User-specific bin color pallete (see [graphical parameters])
##' plot(mnhist_Ind_vil,  idx_order = idx, col = colorRampPalette(colors=c("#FFFFFF","#000000"))(50))
##' 
##' #Users can print p_mat on the plot using user-specific colors
##' plot(mnhist_Ind_vil,  idx_order = idx, prob= TRUE, prob.col = "blue",
##'       col = colorRampPalette(colors=c("#FFFFFF","#000000"))(50))
##' }
##' @importFrom stats heatmap
##' @importFrom graphics text
##' @rdname plot.multinethist
##' @exportS3Method
plot.multinethist <- function(x, y= NA, type = "MNhist",
                         idx_order = 1:max(x$cluster), 
                         prob = FALSE, digits = 2,
                         prob.cex =  0.1 + 0.5/log10(max(x$cluster)),
                         prob.col = "black",
                         ...){
  n_layers <- ifelse(length(dim(x$thetahat)) == 2, 1, dim(x$thetahat)[3])
  
  if(n_layers == 1){
    invisible(plot.nethist(x, type, idx_order, prob, digits,
                           prob.cex, prob.col, ...))
  }else{
    invisible(plot_MNhist_Mlayers(x, type, idx_order, prob, digits,
                                  prob.cex, prob.col, ...))
  }
}
plot_MNhist_Mlayers <- function(x, type = "MNhist",
                              idx_order = 1:max(x$cluster), 
                              prob = FALSE, digits = 2,
                              prob.cex =  0.1 + 0.5/log10(max(x$cluster)),
                              prob.col = "black",
                              ...){
  k<-dim(x$thetahat)[1]
  if(!(type %in% c("MNhist", "prob"))){
    stop("type must be one of MNhist or prob.")
  }
  if(!.is_valid_order(idx_order, 1:k)){
    warning(paste0("idx_order is invalid. Set idx_order = 1:",k))
    idx_order <- 1:k
  }
  
  #If we want to draw homogeneous multinetwork histogram 
  #(that is, x$homegenous == TRUE, type =="MNhist"),
  #just draw 1st layer MNhist.
  n_loop <- ifelse(x$homogeneous & (type=="MNhist"), 
                   1, 
                   dim(x$thetahat)[3])
  for(l in 1:n_loop){
    mat <- switch(type,
                  "MNhist" = (x$thetahat[idx_order, idx_order, l])/x$rho_hat[l],
                  "prob" = x$thetahat[idx_order, idx_order, l])
    rownames(mat) <- idx_order
    colnames(mat) <- idx_order
    if(prob & (type=="prob")){
      tryCatch(
        heatmap(mat, Rowv = NA, symm = TRUE,
                main = paste("Layer ", l),
                add.expr = {text(rep(1:k,each=k), rev(rep(1:k,k)),
                                 round(as.vector(mat), digits),
                                 cex = prob.cex, col = prob.col)},
                ...),
        error = function(e) {
          if (!grepl("pin", conditionMessage(e))) stop(e)
        }
      )
    }else{
      tryCatch(
        heatmap(mat, Rowv = NA, symm = TRUE, main = paste("Layer ", l), ...),
        error = function(e) {
          if (!grepl("pin", conditionMessage(e))) stop(e)
        }
      )
    }
  }
}
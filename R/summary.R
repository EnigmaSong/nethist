##' Bin summary by covariate
##'
##' Drawing a bin summary plot of covariates given an `multinethist` object with an user-specified order.
##'
##' @param object a `multinethist` object from [multinethist()].
##' @param covariate a vector for univariate covariate. If it is factor, draw a stacked barchart. If it is numeric, draw a violin plot.
##' @param idx_order A numeric vector for index label order, which must be a permutation of `object$cluster`. If `NA`, it uses `1:max(object$clsuter)`. 
##' @param main title of summary plot. If NA, the plot has no title.
##' @param ylab label of y-axis. If NA, y-axis label is "covariate"
##' @param legend_title title of legend. If NA, the legend title is "covariate"
##' @param stat variables pass to [ggplot2::geom_bar()]. Only used for a factor covariate.
##' @param position variables pass to [ggplot2::geom_bar()]. Only used for a factor covariate.
##' @details 
##' ... includes various [`graphical parameters`] passes to [stats::heatmap()], then [graphics::image()]. 
##' @returns 
##' a heatmap of network histogram or `p_mat` ordered by `idx_order` from ``nethist`` object.
##' @import ggplot2
##' @examples
##' \dontrun{
##' set.seed(42)
##' mnets <- rnets_graphon(3, 200, function(x,y) 0.2*pmin(x,y))
##' 
##' mnhist <- multinethist(mnets)
##' x <- rep(1:2, 100)
##' summary_plot(mnets, x)
##' }
##' @export
summary_plot <- function(object, covariate,
                          idx_order = 1:max(object$cluster),
                          main = NA,
                          ylab = NA,
                          legend_title = NA,
                          stat = "count", position = "stack"){
  UseMethod("summary_plot")
}
##' @exportS3Method 
summary_plot.multinethist <- function(object, covariate,
                            idx_order = 1:max(object$cluster),
                            main = NA,
                            ylab = NA,
                            legend_title = NA,
                            stat = "count", position = "stack"){
  k = max(object$cluster)
  if(!.is_valid_order(idx_order, 1:k)){
    warning(paste0("idx_order is invalid. Set idx_order = 1:",k))
    idx_order <- 1:k
  }
  
  df <- data.frame(cluster = factor(object$cluster, levels = idx_order),
                   covariate = covariate)
  if(is.factor(covariate)){
    p <- summary_covariate_factor(df, legend_title, stat, position)
  }else if(is.vector(covariate) & is.numeric(covariate)){
    p <- summary_covariate_numeric(df)
  }else{
    stop("covariate must be either a numeric vector or a factor.")
  }
  
  if(!is.na(main)) p <- p + ggplot2::ggtitle(main)
  if(!is.na(ylab)) p <- p + ggplot2::ylab(ylab)
  p <- p + ggplot2::theme(
    panel.background = ggplot2::element_rect(fill='transparent'),
    plot.background = ggplot2::element_rect(fill='transparent', color=NA),
    panel.grid.major = ggplot2::element_blank(), 
    panel.grid.minor = ggplot2::element_blank(), 
    legend.background = ggplot2::element_rect(fill='transparent'), 
    legend.box.background = ggplot2::element_rect(fill='transparent') 
  )
  print(p)  
}

summary_covariate_factor <- function(df, legend_title, stat, position){
  p <- ggplot2::ggplot(data=df, mapping=ggplot2::aes(x = cluster, fill = covariate)) + ggplot2::geom_bar(position = position, stat = stat)
  if(!is.na(legend_title)) p <- p + ggplot2::labs(fill = legend_title)
  p
}

summary_covariate_numeric <- function(df){
  p <- ggplot2::ggplot(df, ggplot2::aes(x = cluster, y = covariate)) + ggplot2::geom_violin()
}

globalVariables(c("cluster", "covariate")) ## to remove warning message (https://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when)
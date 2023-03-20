##' Bin summary by covariate
##'
##' Drawing a bin summary plot of covariates given an `nethist` object with an user-specified order.
##'
##' @param covariate a vector for univariate covariate. If it is factor, draw a stacked barchart. If it is numeric, draw a violin plot.
##' @param x a nethist object from [nethist()].
##' @param idx_order A numeric vector for index label order, which must be a permutation of `x$cluster`. If `NA`, it uses `1:max(x$clsuter)`. 
##' @param main title of summary plot. If NA, the plot has no title.
##' @param ylab label of y-axis. If NA, y-axis label is "covariate"
##' @param legend_title title of legend. If NA, the legend title is "covariate"
##' @param stat, position variables pass to [ggplot2::geom_bar()]. Only used for a factor covariate.
##' @details 
##' ... includes various [`graphical parameters`] passes to [stats::heatmap()], then [graphics::image()]. 
##' @returns 
##' a heatmap of network histogram or `p_mat` ordered by `idx_order` from ``nethist`` object.
##' @examples
##' \dontrun{
##' set.seed(42)
##' nethist_polblog <- nethist(polblog, h = 72)
##' 
##' #Add factor covariate
##' politic <- factor(c(rep("Liberal", 586), rep("Conservative", 1224-586)))
##' 
##' summary(nethist_polblog, politic)
##' 
##' #with user-specified order
##' idx <-  c(17,10,13,7,12,4,8,5,1,2,6,3,9,15,11,14,16)
##' summary(nethist_polblog, politic, idx_order = idx)
##' 
##' #Add title and change ylab
##' summary(nethist_polblog, politic, idx_order = idx, main = "Political View", ylab = "Number of blogs")
##' 
##' #Add numeric covariates
##' numeric_covariate <- nethist_polblog$cluster + rnorm(1224, sd = 3)
##' 
##' summary(nethist_polblog, numeric_covariate) 
##' summary(nethist_polblog, numeric_covariate, idx_order = idx) #With-user specified order
##' 
##' }
##' @exportS3Method 
##' @export
summary.nethist <- function(x, covariate,
                            idx_order = 1:max(x$cluster),
                            main = NA,
                            ylab = NA,
                            legend_title = NA,
                            stat = "count", position = "stack"){
  df <- data.frame(cluster = factor(x$cluster, level = idx_order),
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

  print(p)  
}

summary_covariate_factor <- function(df, legend_title, stat, position){
  p <- ggplot2::ggplot(df, ggplot2::aes(x = cluster, fill = covariate)) + ggplot2::geom_bar(position = position, stat = stat)
  if(!is.na(legend_title)) p <- p + ggplot2::labs(fill = legend_title)
  p
}

summary_covariate_numeric <- function(df){
  p <- ggplot2::ggplot(df, ggplot2::aes(x = cluster, y = covariate)) + ggplot2::geom_violin()
}
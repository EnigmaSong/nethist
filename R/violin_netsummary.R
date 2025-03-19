##' Network summary plots
##'
##' Draw a network summary plot proposed by Maugis et al. (2017). To count k-cycles, Alon et al. (1997) is used.
##'
##' @param A an adjacency matrix or igraph object to draw a network summary plot. It must be an undirected and simple graph.
##' @param subsample_sizes a numeric vector of node subsample sizes. If `NA`, a length `Ns` vector is obtained from the automatic subsample size selection.
##' @param max_cycle_order an integer value of the maximum cycle size. Must be `>=3` and `<=7`.
##' @param R an integer value of subsampling replication. If `NA`, `R` is automatically selected by `alpha`.
##' @param Ns number of different subsample sizes. It is only used when `subsample_sizes = NA`, that is, when automatic subsample size selection is used.
##' @param alpha a pre-specified level used in determining `R` and `subsample_sizes` when they are not specified. It must be in (0,1). Default is 0.05. Smaller `alpha` gives larger `R` and `subsample_sizes`.
##' @param y.max Upper limit of y-axis of the plot. Must be 0 < `y.max` <= 1. If `NA`, the upper limit is automatically selected.
##' @param save.plot logical variable whether save the generated figure or not. If `TRUE`, the plot is saved by [ggplot2::ggsave()] in the specified file name. Otherwise, display the generated plot.
##' @param filename file name to save the generated figure. 
##' @param width a numeric value of the width of the generated figure in inch. It is only used when `save.plot = TRUE`.
##' @param height a numeric value of the height of the generated figure in inch. It is only used when `save.plot = TRUE`.
##' @return 
##' A network summary plot, and a data.frame about the networks summaries.
##' @details 
##' Vertex sampling is done by simple random sampling without replacement.
##' 
##' Following matrix classes are supported: [base::matrix], [Matrix::dgCMatrix-class]
##' @references Maugis et al. (2017). Topology reveals universal features for network comparison. arXiv: 1705.05677 
##' @references Alon et al. (1997). Finding and counting given length cycles. Algorithmica 17, 209–223 (1997). https://doi.org/10.1007/BF02523189
##' @examples
##' {
##' set.seed(2022)
##' #Generating Erdos-Renyi graph
##' n <- 400
##' #igraph object
##' A <- igraph::sample_gnp(n, 0.05)
##' violin_netsummary(A)
##' 
##' #sparse adjacency matrix
##' A2 <- igraph::as_adj(A)
##' violin_netsummary(A2)
##' 
##' #dense adjacency matrix
##' A2 <- igraph::as_adj(A, sparse = FALSE)
##' violin_netsummary(A2)
##' 
##' #user-specified R and subsample_sizes
##' violin_netsummary(A, R = 500, subsample_sizes = 150)
##'
##' #user-specified alpha
##' violin_netsummary(A, alpha = 0.1)
##'
##' #saving the plot with user-specified file name
##' \dontrun{
##' violin_netsummary(A, save.plot = TRUE, filename = "myfig.pdf")
##' }
##' }
##' @importFrom ggtext element_markdown
##' @import png 
##' @export
##' 
violin_netsummary <- function(A,
                              subsample_sizes = NA, 
                              max_cycle_order = 4, 
                              R=NA, 
                              Ns = 11, alpha = 0.05,
                              y.max=NA, save.plot = FALSE, 
                              filename = "myplot.pdf", 
                              width = 7, height = 5){
  UseMethod("violin_netsummary")
}
##' @exportS3Method
violin_netsummary.igraph<- function(A, 
                                     subsample_sizes, 
                                     max_cycle_order, 
                                     R, Ns, alpha,
                                     y.max, save.plot, 
                                     filename, width, height){
  args <- as.list(environment())
  args$A<- igraph::as_adjacency_matrix(args$A, sparse = FALSE)
  
  do.call("violin_netsummary.default", args = args)
}

##' @exportS3Method
violin_netsummary.matrix<- function(A, 
                                    subsample_sizes, 
                                    max_cycle_order, 
                                    R, Ns, alpha,
                                    y.max, save.plot, 
                                    filename, width, height){
  args <- as.list(environment())
  do.call("violin_netsummary.default", args = args)
}

##' @exportS3Method
violin_netsummary.dgCMatrix<- function(A, 
                                    subsample_sizes, 
                                    max_cycle_order, 
                                    R, Ns, alpha,
                                    y.max, save.plot, 
                                    filename, width, height){
  args <- as.list(environment())
  args$A <- as.matrix(args$A)
  do.call("violin_netsummary.default", args = args)
}

##' @exportS3Method 
violin_netsummary.default<- function(A, 
                                     subsample_sizes = NA, 
                                     max_cycle_order = 4, 
                                     R=NA, Ns = 11, alpha = 0.05,
                                     y.max=NA, save.plot = FALSE, 
                                     filename = "myplot.pdf", width = 7, height = 5){
  if(!.is_undirected_simple(A)) stop("Network A must be an undirected simple network.")
  
  if((max_cycle_order < 3)|(max_cycle_order%%1 != 0)){
    stop("order_cycle must be >= 3 integers.")
    max_cycle_order <- 7
  }else if(max_cycle_order >7){
    message("order_cycle >7 is not implemented. Use an integer between 3 and 7.")
    max_cycle_order <- 7
  }
  if(is.na(R)){
    R <- ceiling((1/(2*alpha)*stats::qnorm(1-alpha/(2*(max_cycle_order-1))))^2)
    message(paste("Use R=",R))
  }
  if(is.na(subsample_sizes)){
    subsample_sizes <- auto_select_subsample_sizes(A, Ns, k_max = max_cycle_order, R, alpha=0.05, delta = 0.05)
  }
  
  result <- .net_summary_subsample_adj(A, subsample_sizes, max_cycle_order, R)
  colnames(result) <- c("v-shape","triangle","square","pentagon", "hexagon", 'septagon')[1:(max_cycle_order-1)]
  if(!is.na(y.max) & ((y.max > 1)|(y.max < 0))){
    warning("y.max: Use a number between 0 and 1")
    y.max = NA
  }
  if(is.na(y.max)){
    y.max <- max(result)
  }
  
  result <- data.frame(result)
  melt_result <- suppressMessages(reshape2::melt(result, 
                           value.name= "summaries", variable.name= "subgraphs"))
  xlabels <- get_netsummary_xlables(max_cycle_order)
  p <- ggplot2::ggplot(melt_result, ggplot2::aes(.data$subgraphs, .data$summaries))
  p <- p + ggplot2::geom_violin() + ggplot2::ylim(0,y.max) 
  p <- p + ggplot2::ylab("Prevalence and local variability") 
  p <- p + ggplot2::scale_x_discrete(name=NULL, labels=xlabels)
  p <- p + ggplot2::theme(axis.text.x=ggtext::element_markdown(color="black",size=rel(1)))
  
  if(save.plot){
    ggplot2::ggsave(filename,width = width, height = height, unit = "in")
  }else{
    print(p)
  }
  return(invisible(result))
}

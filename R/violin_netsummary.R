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
##' @return 
##' A network summary plot, and a data.frame about the networks summaries.
##' @details 
##' Vertex sampling is done by simple random sampling without replacement.
##' 
##' Following matrix classes are supported: [base::matrix], [Matrix::dgCMatrix-class]
##' @references Maugis et al. (2017). Topology reveals universal features for network comparison. arXiv: 1705.05677 
##' @references Alon et al. (1997). Finding and counting given length cycles. Algorithmica 17, 209–223 (1997). https://doi.org/10.1007/BF02523189
##' @examples
##' \dontrun{
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
##' violin_netsummary(A, save.plot = TRUE, filename = "myfig.pdf")
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
                              filename = "myplot.pdf"){
  UseMethod("violin_netsummary")
}
##' @exportS3Method
violin_netsummary.igraph<- function(A, 
                                     subsample_sizes, 
                                     max_cycle_order, 
                                     R, Ns, alpha,
                                     y.max, save.plot, 
                                     filename){
  args <- as.list(environment())
  args$A<- igraph::as_adj(args$A, sparse = FALSE)
  
  do.call("violin_netsummary.default", args = args)
}

##' @exportS3Method
violin_netsummary.matrix<- function(A, 
                                    subsample_sizes, 
                                    max_cycle_order, 
                                    R, Ns, alpha,
                                    y.max, save.plot, 
                                    filename){
  args <- as.list(environment())
  do.call("violin_netsummary.default", args = args)
}

##' @exportS3Method
violin_netsummary.dgCMatrix<- function(A, 
                                    subsample_sizes, 
                                    max_cycle_order, 
                                    R, Ns, alpha,
                                    y.max, save.plot, 
                                    filename){
  args <- as.list(environment())
  args$A <- as.matrix(args$A)
  do.call("violin_netsummary.default", args = args)
}

violin_netsummary.default<- function(A, 
                                     subsample_sizes = NA, 
                                     max_cycle_order = 4, 
                                     R=NA, Ns = 11, alpha = 0.05,
                                     y.max=NA, save.plot = FALSE, 
                                     filename = "myplot.pdf"){
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
  colnames(result) <- c("tree","triangle","square","pentagon", "hexagon", 'septagon')[1:(max_cycle_order-1)]
  if(!is.na(y.max) & ((y.max > 1)|(y.max < 0))){
    warning("y.max: Use a number between 0 and 1")
    y.max = NA
  }
  if(is.na(y.max)){
    y.max <- max(result)
  }

  suppressMessages(result <- reshape2::melt(data.frame(result)))
  p <- ggplot2::ggplot(result, ggplot2::aes(variable, value))
  p <- p + ggplot2::geom_violin() + ggplot2::ylim(0,y.max) + ggplot2::ylab("Prevalence and local variability") + ggplot2::xlab("")
  if(save.plot){
    ggplot2::ggsave(filename,width = 7, height = 5, unit = "in")
  }else{
    print(p)
  }
  return(invisible(result))
}

auto_select_subsample_sizes <- function(A, Ns, k_max, R, alpha=0.05, delta){
  n <- dim(A)[1]

  s_star <- min(max(k_max + 1, min(floor(n/4), 3*(k_max+1))), n)/(1+delta)
  K_set <- 2:k_max
  s_max <- n

  for(i in 1:ceiling(log(n)/log(1+delta))){
    s_star <- min(ceiling((1+delta)*s_star), n)
    t_k <- .net_summary_subsample_adj(A = A, subsample_sizes = s_star,
                                     max_cycle_order = k_max, R = R)
    #Check summary separated from 0
    R_cols <- apply(t_k,2,function(x) sum(!is.na(x)))
    colSums_t_k <- colSums(t_k, na.rm=TRUE)
    colSums_t_k_square <- colSums(t_k^2, na.rm=TRUE)
    p_k_numer <- 1/R*colSums_t_k
    p_k_denom <- sqrt(1/(R_cols-1)*colSums_t_k_square - 1/(R_cols*(R_cols-1))*(colSums_t_k)^2)
    p_k_denom[which(p_k_denom==0)] <- 1 #denominator is 0 iff all t_k are zero.
    p_k <- pnorm(p_k_numer/p_k_denom, lower.tail = FALSE)
    p <- max(p_k[K_set-1]) #quantify least-separated summary

    if((p <= alpha/(k_max-1))|(s_star >= s_max)){
      if(s_max == floor(0.8*n)){
        subsample_sizes <- round(seq(0.9,1.1, length.out = Ns)*s_star)
        #Reset and halt
        return(subsample_sizes)
      }
      s_star <- min(max(k_max + 1, min(floor(n/4), 3*(k_max+1))), n)/(1+delta)
      K_set <- K_set[p_k < 1/2] #ignore all-zero summaries
      s_max <- floor(0.8*n) #restrict maximum subgraph size
    }
  }

  return(subsample_sizes)
}

globalVariables(c("variable","value"))

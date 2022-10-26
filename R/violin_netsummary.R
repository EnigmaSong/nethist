##' Network summary plots
##'
##' Draw a network summary plot proposed by Maugis et al. (2017). Counting k-cycles by Alon et al. (1997).
##'
##' @param A an adjacency matrix to draw a network summary plot.
##' @param Ns number of different subsample size
##' @param subsample_sizes an integer value of node subsample size. 
##' @param max_cycle_order an integer value of the maximum cycle size. Must be >=3 and <=7.
##' @param R an integer value of subsampling replication. If NA, automatically selected using normal CDF
##' @param alpha level used 
##' @param y.max maximum value of y. Must be 0 < y.max <= 1
##' @param save.plot logical variable whether save the generated figure or not. If TRUE, save the generated plot in the specified file name. Otherwise, display the generated plot.
##' @param filename file name to save the generated figure
##' @return A network summary plot as a violin plot 
##' @details Vertex sampling is done by simple random sampling without replacement.
##' @references Maugis et al. (2017). Topology reveals universal features for network comparison. arXiv: 1705.05677 
##' @references Alon et al. (1997). Finding and counting given length cycles. Algorithmica 17, 209–223 (1997). https://doi.org/10.1007/BF02523189
##' @examples
##' \dontrun{
##' set.seed(2022)
##' #Generating Erdos-Renyi graph
##' n <- 400
##' A <- igraph::sample_gnp(n, 0.05)
##' A <- igraph::as_adj(A)
##' violin_netsummary(A, save.plot = FALSE)
##' }
##' @export
violin_netsummary <- function(A, 
                              Ns = 11, subsample_sizes = NA, 
                              max_cycle_order = 4, 
                              R=NA, alpha = 0.05,
                              y.max=NA, save.plot = FALSE, 
                              filename = "myplot.pdf"){
  cl <- match.call()
  if(!is.matrix(A)){
    stop(paste("A must be a (dense) matrix. Use as.matrix() to convert the object:", class(A)))
  }
  if(any((A!=0)&(A!=1))){
    message("There are entries neither 0 nor 1. Convert A into a binary matrix.")
    A <- (A!=0)
  }
  if(FALSE){
    
  }
  if((max_cycle_order < 3)|(max_cycle_order%%1 != 0)){
    stop("order_cycle must be >= 3 integers.")
    max_cycle_order <- 7
  }else if(max_cycle_order >7){
    message("order_cycle >7 is not implemented. Use an integer between 3 and 7.")
    max_cycle_order <- 7
  }
  if(is.na(R)){
    R <- ceiling((1/(2*alpha)*stats::qnorm(1-alpha/(2*(max_cycle_order-1))))^2)
    message(paste("R=",R))
  }
  if(is.na(subsample_sizes)){
    subsample_sizes <- auto_select_subsample_sizes(A, Ns, k_max = max_cycle_order, R, alpha=0.05, delta = 0.05)
  }
  
  result <- net_summary_subsample_adj(A, subsample_sizes, max_cycle_order, R)

  if(!is.na(y.max) & ((y.max > 1)|(y.max < 0))){
    warning("y.max: Use a number between 0 and 1")
    y.max = NA
  }
  if(is.na(y.max)){
    y.max <- max(result)
  }

  result <- reshape2::melt(data.frame(result))  
  p <- ggplot2::ggplot(result, ggplot2::aes(variable, value))
  p <- p + ggplot2::geom_violin() + ggplot2::ylim(0,y.max) + ggplot2::ylab("Prevalence and local variability") + ggplot2::xlab("")
  if(save.plot){
    ggplot2::ggsave(filename,width = 4, height = 4, unit = "in")
  }else{
    print(p)
  }
  return(result)
}

net_summary_subsample_adj <- function(A, subsample_sizes, max_cycle_order, R){
  n <- dim(A)[1]
  Ns <- length(subsample_sizes)
  
  
  result <- matrix(0, R, max_cycle_order-1)
  colname_summary <- c("Trees", 'Triangles', 'Squares', "Pentagons", "Hexagons", "Septagons")
  colnames(result) <- colname_summary[1:(max_cycle_order-1)]
  
  ind_cycles <- 2:(max_cycle_order-1)
  powers <- 1/seq(3,max_cycle_order,1)

  for(r in 1:R){
    subsample_size <- subsample_sizes[ceiling(Ns*r/R)]
    denoms <- sapply(3:max_cycle_order, function(k, n) .ffct(subsample_size,k)/(2*k), n=subsample_size)
    sample_ind <- sample.int(n, subsample_size)
    sub_A <- A[sample_ind, sample_ind]
    deg_A <- colSums(sub_A)
    result[r,1] <- (1/2*sum(deg_A*(deg_A-1))) / ((subsample_size-2)/2*sum(deg_A))# (S.57) in Maugis et al. (2017)
    # Compute k-cycle ratio defined in (S.58)
    numers <- .count_k_cycle(sub_A, max_cycle_order)
    result[r,ind_cycles] <- (numers/denoms)^powers
  }
  return(result)
}

auto_select_subsample_sizes <- function(A, Ns, k_max, R, alpha=0.05, delta){
  n <- dim(A)[1]
  
  s_star <- min(max(k_max + 1, min(floor(n/4), 3*(k_max+1))), n)/(1+delta)
  K_set <- 2:k_max
  s_max <- n
  
  for(i in 1:ceiling(log(n)/log(1+delta))){
    s_star <- min(ceiling((1+delta)*s_star), n)
    message(paste(i,"th iter, s_star=", s_star))
    t_k <- net_summary_subsample_adj(A = A, subsample_sizes = s_star, 
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
      
      print(K_set)
    }
  }
  
  return(subsample_sizes)
}

#' Automatically Selects Subsample Sizes for Network Analysis
#'
#' This function determines optimal subsample sizes for network analysis based on 
#' statistical summaries, ensuring that key network properties remain distinguishable.
#'
#' @param A A square adjacency matrix representing the network.
#' @param Ns Integer. The number of subsample sizes to generate around the optimal size.
#' @param k_max Integer. The maximum cycle order to consider in network summaries.
#' @param R Integer. The number of random subsampling iterations for statistical estimation.
#' @param alpha Numeric. Significance level for hypothesis testing (default is 0.05).
#' @param delta Numeric. Growth factor for subsample sizes, controlling step increments.
#'
#' @return A numeric vector of `Ns` subsample sizes, determined based on statistical properties.
#'
#' @details
#' The function iteratively increases subsample sizes while checking whether key 
#' network summaries significantly deviate from zero. If the computed p-value 
#' (quantifying the least-separated summary) is below `alpha / (k_max - 1)`, 
#' or if the subsample size exceeds a threshold, the function returns the selected sizes.
#'
#'##' @noRd
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

##' Get the figure lables
##' @noRd
get_netsummary_xlables <- function(max_cycle_order){
  fig_path <- system.file("violin_summary",package="nethist")
  subgraphs <- c("v_shape","triangle","square","pentagon", "hexagon", 'septagon')
  
  
  labels <- paste0("<img src='",fig_path,"/",
                   subgraphs[ 1:(max_cycle_order-1)], ".png' width='20' />")
  
  return(labels)
}
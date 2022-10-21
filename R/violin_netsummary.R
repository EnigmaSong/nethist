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
    deg_A <- colSums(A[sample_ind, sample_ind])
    result[r,1] <- (1/2*sum(deg_A*(deg_A-1))) / ((subsample_size-2)/2*sum(deg_A))# (S.57) in Maugis et al. (2017)
    # Compute k-cycle ratio defined in (S.58)
    numers <- count_k_cycle(A[sample_ind, sample_ind], max_cycle_order)
    result[r,ind_cycles] <- (numers/denoms)^powers
  }
  return(result)
}

#Counting k-cycles
count_k_cycle<- function(A, max_cycle_order){
  n <- nrow(A)
  deg_A <- rowSums(A)
  nonzero_deg <- which(deg_A > 0)
  A <- A[nonzero_deg,nonzero_deg]
  deg_A <- deg_A[nonzero_deg]
  edge_num_A <- sum(A)/2
  count <- rep(0, max_cycle_order-2)

  A2 <- A%*%A
  dA2 <- diag(A2)
  Ak <- A2
  
  # Counting k-cycles is from Alon et al. (1997)
  for(k in 3:max_cycle_order){
    Ak <- Ak %*% A
    if(k==3){
      A3 <- Ak
      dA3 <- diag(A3)
      tA3 <- sum(dA3)
      sA3 <- sum(A3)
      count[1] <- tA3/6
    }else if(k==4){
      A4 <- Ak
      dA4 <- diag(A4)
      tA4 <- sum(dA4)
      
      count[2]<- tA4 + 2*edge_num_A - 2*sum(deg_A^2)
      count[2]<- count[2]/8
    }else if(k==5){
      A5 <- Ak
      dA5 <- diag(A5)
      tA5 <- sum(dA5)
      
      count[3]<- tA5 - 5*sum((deg_A-1)*dA3)
      count[3]<- count[3]/10
    }else if(k==6){
      A6 <- Ak
      dA6 <- diag(A6)
      tA6 <- sum(dA6)
      
      count[4]<- tA6 - 3*sum(dA3^2)
      count[4]<- count[4] + 9*sum((A2^2)*A)
      count[4]<- count[4] - 6*sum(dA4 * (deg_A-1))
      count[4]<- count[4] - 4*sum(dA3 - deg_A^3)
      count[4]<- count[4] + 3*sA3
      count[4]<- count[4] - 12*sum(deg_A^2)
      count[4]<- count[4] + 4*sum(deg_A)
      count[4]<- count[4]/12
    }else if(k==7){
      dA7 <- diag(Ak)
      tA7 <- sum(dA7)
      count[5] <- tA7 - 7*sum(dA3*dA4) + 7*sum((A2^3) * A)
      count[5] <- count[5] - 7*sum(dA5*deg_A)
      count[5] <- count[5] + 21*sum(A3*A2*A)
      count[5] <- count[5] + 7*tA5
      count[5] <- count[5] - 28*sum((A2^2)*A)
      count[5] <- count[5] + 7*sum(A2*A*outer(deg_A,deg_A))
      count[5] <- count[5] + 14*sum(dA3*deg_A^2)
      count[5] <- count[5] + 7 *sum(dA3*rowSums(A2))
      count[5] <- count[5] - 77*sum(dA3*deg_A)
      count[5] <- count[5] + 56*tA3
      count[5] <- count[5]/14
    }else if(k==8){
      ### Not implemented completely as of 2022/09/10
      count[6] <-  sum(diag(Ak)) - 4*sum(dA4^2) - 8*sum(dA3*dA5) - 8*sum(dA2*dA6)
      count[6] <- count[6] + 16*sum(dA2*(dA3^2)) + 8*sum(dA6) + 16*sum(dA4*(dA2^2))
      count[6] <- count[6] - 72*sum(dA3^2) - 96*sum(dA4*dA2) - 12*sum(dA2^4)
      count[6] <- count[6] + 64*sum(dA2*dA3) + 73*sum(dA4) + 72*sum(dA2^3) - 112*sum(dA3) + 36*sum(dA2)
      # inter = (
      #   np.sum(np.diag(Ak)) - 4 * np.sum(dA4 * dA4) - 8 * np.sum(dA3 * dA5) - 8 * np.sum(dA2 * dA6)
      #   + 16 * np.sum(dA2 * dA3 * dA3) + 8 * np.sum(dA6) + 16 * np.sum(dA4 * dA2 * dA2)
      #   - 72 * np.sum(dA3 * dA3) - 96 * np.sum(dA4 * dA2) - 12 * np.sum(dA2 * dA2 * dA2 * dA2)
      #   + 64 * np.sum(dA2 * dA3) + 73 * np.sum(dA4) + 72 * np.sum(dA2 * dA2 * dA2)
      #   - 112 * np.sum(dA3) + 36 * np.sum(dA2)
      # )
      count[6] <- count[6] + 2*sum(A2^4) + 24*sum((A2^2)*A*A3) + 4*sum(outer(dA3,dA3)*A)
      count[6] <- count[6] + 16*sum(outer(dA2,dA3)*A*A2) + 12*sum(A*(A3^2)) + 24*sum(A*A4*A2)
      count[6] <- count[6] + 4*sum(outer(dA2,dA2)*(A2^2)) + 8*sum(A2%*%diag(dA4)) + 8*sum(outer(dA2,dA2)*A*A3)#
      count[6] <- count[6] - 16*sum(A2^3) - 32*sum(A*A2*A3) - 96*sum(dA2*diag((A2^2)*A))
      count[6] <- count[6] - 4*sum(A4) - 16*sum((dA2^2)*A2)#
      count[6] <- count[6] + 272*sum(A2*A2*A) + 48*sum(A3) - 132*sum(A2)
      # inter += (
      #   2 * np.sum(A2 * A2 * A2 * A2) + 24 * np.sum(A2 * A2 * A1 * A3) + 4 * np.sum(np.outer(dA3, dA3) * A1)
      #   + 16 * np.sum(np.outer(dA2, dA3) * A1 * A2) + 12 * np.sum(A1 * A3 * A3) + 24 * np.sum(A1 * A4 * A2)
      #   + 4 * np.sum(np.outer(dA2, dA2) * A2 * A2) + 8 * np.sum(A2 @ np.diag(dA4)) + 8 * np.sum(np.outer(dA2, dA2) * A1 * A3)
      #   - 16 * np.sum(A2 * A2 * A2) - 32 * np.sum(A1 * A2 * A3) - 96 * np.sum(np.diag(dA2) * (A2 * A2 * A1))
      #   - 4 * np.sum(A4) - 16 * np.sum(np.diag(dA2 * dA2) @ A2)
      #   + 272 * np.sum(A2 * A2 * A1) + 48 * np.sum(A3) - 132 * np.sum(A2)
      # )
      count[6] <- count[6] - 64*sum(A*((A*A2)^2)) - 24*sum(A*(A%*%diag(dA2)%*%A)*A2)  #ADA = (sum_{k=1}^n d_k a_{ik}a_{kj})_{i,j}
      xk4 <- 0
      for(i in 1:length(nonzero_deg)){
        xk4 <- xk4 + A[i,]%*%(A*(A%*%diag(A[i,])%*%A))%*%A[,i]
      }
      count[6] <- count[6] + 22*xk4
      count[6] <- count[6]/16
      # inter += -64 * np.sum(A1 * ((A1 * A2) ** 2)) - 24 * np.sum(
      #   A1 * (A1 @ np.diag(dA2) @ A1) * A2
      # )
      # xk4 = 0
      # for i in range(A1.shape[0]):
      #   xk4 += A1[i, :] @ (A1 * (A1 @ np.diag(A1[i, :]) @ A1)) @ (A1[i, :].T)
      # inter += 22 * xk4
    }else if(k==9){
      
    }
  }
  return(count)
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

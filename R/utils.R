#' @noRd

# Falling factorial
.ffct <- function(x, k) {
  sapply(x, function(y,k) {
    prod(y:(y - k + 1))
  }, k= k)
}
.count_k_cycle_R<-function(A, max_cycle_order){
  n <- dim(A)[1]
  deg_A <- rowSums(A)
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
    }
  }
  return(count)
}

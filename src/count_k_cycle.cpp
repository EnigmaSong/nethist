#include <RcppArmadillo.h>
using namespace Rcpp;
const double eps = arma::datum::eps;

// Rcpp implementation for counting k-cycles

arma::vec count_k_cycle_cpp(arma::mat &A, int max_cycle_order) {
  arma::vec count(max_cycle_order-1);
  arma::vec deg_A;
  double edge_num_A;
  
  arma::mat A2 = A*A;
  arma::mat A3 = A2*A;
  
  // Count C_3
  count[0] = trace(A3)/6.0;
  // Count C_4
  if(max_cycle_order >= 4){
    arma::mat A4 = A3*A;
    arma::vec dA4 = A4.diag();
    tA4 <- A4.sum();
    
  }
  
  // Count C_5
  if(max_cycle_order >= 5){
    
  }
  
  // Count C_6
  if(max_cycle_order >= 6){
    
  }
  
  // Count C_6
  if(max_cycle_order >= 7){
    
  }
  
  
  return count;
}

// #Counting k-cycles
// count_k_cycle<- function(A, max_cycle_order){
//   n <- nrow(A)
//   deg_A <- rowSums(A)
//   nonzero_deg <- which(deg_A > 0)
//   A <- A[nonzero_deg,nonzero_deg]
//   deg_A <- deg_A[nonzero_deg]
//   edge_num_A <- sum(A)/2
//   count <- rep(0, max_cycle_order-2)
//   
//   A2 <- A%*%A
//   dA2 <- diag(A2)
//   Ak <- A2
//   
// # Counting k-cycles is from Alon et al. (1997)
//   for(k in 3:max_cycle_order){
//     Ak <- Ak %*% A
//     if(k==3){
//       A3 <- Ak
//       dA3 <- diag(A3)
//       tA3 <- sum(dA3)
//       sA3 <- sum(A3)
//       count[1] <- tA3/6
//     }else if(k==4){
//       A4 <- Ak
//       dA4 <- diag(A4)
//       tA4 <- sum(dA4)
//       
//       count[2]<- tA4 + 2*edge_num_A - 2*sum(deg_A^2)
//       count[2]<- count[2]/8
//     }else if(k==5){
//       A5 <- Ak
//       dA5 <- diag(A5)
//       tA5 <- sum(dA5)
//       
//       count[3]<- tA5 - 5*sum((deg_A-1)*dA3)
//       count[3]<- count[3]/10
//     }else if(k==6){
//       A6 <- Ak
//       dA6 <- diag(A6)
//       tA6 <- sum(dA6)
//       
//       count[4]<- tA6 - 3*sum(dA3^2)
//       count[4]<- count[4] + 9*sum((A2^2)*A)
//       count[4]<- count[4] - 6*sum(dA4 * (deg_A-1))
//       count[4]<- count[4] - 4*sum(dA3 - deg_A^3)
//       count[4]<- count[4] + 3*sA3
//       count[4]<- count[4] - 12*sum(deg_A^2)
//       count[4]<- count[4] + 4*sum(deg_A)
//       count[4]<- count[4]/12
//     }else if(k==7){
//       dA7 <- diag(Ak)
//       tA7 <- sum(dA7)
//       count[5] <- tA7 - 7*sum(dA3*dA4) + 7*sum((A2^3) * A)
//       count[5] <- count[5] - 7*sum(dA5*deg_A)
//       count[5] <- count[5] + 21*sum(A3*A2*A)
//       count[5] <- count[5] + 7*tA5
//       count[5] <- count[5] - 28*sum((A2^2)*A)
//       count[5] <- count[5] + 7*sum(A2*A*outer(deg_A,deg_A))
//       count[5] <- count[5] + 14*sum(dA3*deg_A^2)
//       count[5] <- count[5] + 7 *sum(dA3*rowSums(A2))
//       count[5] <- count[5] - 77*sum(dA3*deg_A)
//       count[5] <- count[5] + 56*tA3
//       count[5] <- count[5]/14
//     }else if(k==8){
// ### Not implemented completely as of 2022/09/10
//       count[6] <-  sum(diag(Ak)) - 4*sum(dA4^2) - 8*sum(dA3*dA5) - 8*sum(dA2*dA6)
//       count[6] <- count[6] + 16*sum(dA2*(dA3^2)) + 8*sum(dA6) + 16*sum(dA4*(dA2^2))
//       count[6] <- count[6] - 72*sum(dA3^2) - 96*sum(dA4*dA2) - 12*sum(dA2^4)
//       count[6] <- count[6] + 64*sum(dA2*dA3) + 73*sum(dA4) + 72*sum(dA2^3) - 112*sum(dA3) + 36*sum(dA2)
// # inter = (
// #   np.sum(np.diag(Ak)) - 4 * np.sum(dA4 * dA4) - 8 * np.sum(dA3 * dA5) - 8 * np.sum(dA2 * dA6)
// #   + 16 * np.sum(dA2 * dA3 * dA3) + 8 * np.sum(dA6) + 16 * np.sum(dA4 * dA2 * dA2)
// #   - 72 * np.sum(dA3 * dA3) - 96 * np.sum(dA4 * dA2) - 12 * np.sum(dA2 * dA2 * dA2 * dA2)
// #   + 64 * np.sum(dA2 * dA3) + 73 * np.sum(dA4) + 72 * np.sum(dA2 * dA2 * dA2)
// #   - 112 * np.sum(dA3) + 36 * np.sum(dA2)
// # )
//       count[6] <- count[6] + 2*sum(A2^4) + 24*sum((A2^2)*A*A3) + 4*sum(outer(dA3,dA3)*A)
//       count[6] <- count[6] + 16*sum(outer(dA2,dA3)*A*A2) + 12*sum(A*(A3^2)) + 24*sum(A*A4*A2)
//       count[6] <- count[6] + 4*sum(outer(dA2,dA2)*(A2^2)) + 8*sum(A2%*%diag(dA4)) + 8*sum(outer(dA2,dA2)*A*A3)#
//       count[6] <- count[6] - 16*sum(A2^3) - 32*sum(A*A2*A3) - 96*sum(dA2*diag((A2^2)*A))
//       count[6] <- count[6] - 4*sum(A4) - 16*sum((dA2^2)*A2)#
//       count[6] <- count[6] + 272*sum(A2*A2*A) + 48*sum(A3) - 132*sum(A2)
// # inter += (
// #   2 * np.sum(A2 * A2 * A2 * A2) + 24 * np.sum(A2 * A2 * A1 * A3) + 4 * np.sum(np.outer(dA3, dA3) * A1)
// #   + 16 * np.sum(np.outer(dA2, dA3) * A1 * A2) + 12 * np.sum(A1 * A3 * A3) + 24 * np.sum(A1 * A4 * A2)
// #   + 4 * np.sum(np.outer(dA2, dA2) * A2 * A2) + 8 * np.sum(A2 @ np.diag(dA4)) + 8 * np.sum(np.outer(dA2, dA2) * A1 * A3)
// #   - 16 * np.sum(A2 * A2 * A2) - 32 * np.sum(A1 * A2 * A3) - 96 * np.sum(np.diag(dA2) * (A2 * A2 * A1))
// #   - 4 * np.sum(A4) - 16 * np.sum(np.diag(dA2 * dA2) @ A2)
// #   + 272 * np.sum(A2 * A2 * A1) + 48 * np.sum(A3) - 132 * np.sum(A2)
// # )
//       count[6] <- count[6] - 64*sum(A*((A*A2)^2)) - 24*sum(A*(A%*%diag(dA2)%*%A)*A2)  #ADA = (sum_{k=1}^n d_k a_{ik}a_{kj})_{i,j}
//       xk4 <- 0
//       for(i in 1:length(nonzero_deg)){
//         xk4 <- xk4 + A[i,]%*%(A*(A%*%diag(A[i,])%*%A))%*%A[,i]
//       }
//       count[6] <- count[6] + 22*xk4
//         count[6] <- count[6]/16
// # inter += -64 * np.sum(A1 * ((A1 * A2) ** 2)) - 24 * np.sum(
// #   A1 * (A1 @ np.diag(dA2) @ A1) * A2
// # )
// # xk4 = 0
// # for i in range(A1.shape[0]):
// #   xk4 += A1[i, :] @ (A1 * (A1 @ np.diag(A1[i, :]) @ A1)) @ (A1[i, :].T)
// # inter += 22 * xk4
//     }else if(k==9){
//       
//     }
//   }
//   return(count)
// }


/*** R
A <- igraph::sample_gnp(100,0.2)
count_k_cycle_cpp(A, 7)
*/

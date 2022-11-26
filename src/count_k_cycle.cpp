#include <RcppArmadillo.h>
using namespace Rcpp;

// Rcpp implementation for counting k-cycles
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.count_k_cycle)]]
arma::vec count_k_cycle(const arma::mat &A, const int &max_cycle_order) {
  arma::vec count(max_cycle_order-2);
  arma::vec deg_A = sum(A,1);
  
  arma::uvec nonzero_v = find(deg_A>0);
  arma::mat A_nonzero = A.submat(nonzero_v,nonzero_v);
  deg_A = deg_A(nonzero_v);
  
  double edge_num_A = accu(A_nonzero)/2.0;
  
  arma::mat A2 = A_nonzero*A_nonzero;
  arma::mat A3 = A2*A_nonzero;
  arma::mat A4, A5, A6, A7;
  arma::vec dA3, dA4, dA5;
  double tA5;
  
  // Count C_3
  double tA3 = trace(A3);
  count[0] = tA3/6.0;
  // Count C_4
  if(max_cycle_order >= 4){
    A4 = A3*A_nonzero;
    
    count[1] = trace(A4) + 2.0*edge_num_A - 2*sum(deg_A%deg_A);
    count[1] = count[1]/8.0;
  }
  
  // Count C_5
  if(max_cycle_order >= 5){
    A5 = A4*A_nonzero;
    dA3 = A3.diag();
    tA5 = trace(A5);
    
    count[2] = tA5 - 5*sum((deg_A-1)%dA3);
    count[2] = count[2]/10.0;
  }
  
  // Count C_6
  if(max_cycle_order >= 6){
    A6 = A5*A_nonzero;
    dA4 = A4.diag();
    
    count[3] = trace(A6) - 3*sum(dA3%dA3);
    count[3] = count[3] + 9*accu(A2%A2%A_nonzero);
    count[3] = count[3] - 6*sum(dA4%(deg_A-1));
    count[3] = count[3] - 4*sum(dA3 - pow(deg_A,3));
    count[3] = count[3] + 3*accu(A3);
    count[3] = count[3] - 12*sum(deg_A%deg_A);
    count[3] = count[3] + 4*sum(deg_A);
    count[3] = count[3]/12.0;
  }
  
  // Count C_7
  if(max_cycle_order >= 7){
    A7 = A6*A_nonzero;
    dA5 = A5.diag();
    
    count[4] = trace(A7) - 7*sum(dA3%dA4) + 7*accu(pow(A2,3) % A_nonzero);
    count[4] = count[4] - 7*sum(dA5%deg_A);
    count[4] = count[4] + 21*accu(A3%A2%A_nonzero);
    count[4] = count[4] + 7*tA5;
    count[4] = count[4] - 28*accu(pow(A2,2)%A_nonzero);
    count[4] = count[4] + 7*accu(A2%A_nonzero%(deg_A*deg_A.t()));
    count[4] = count[4] + 14*sum(dA3%pow(deg_A,2));
    count[4] = count[4] + 7*sum(dA3%sum(A2,1));
    count[4] = count[4] - 77*sum(dA3%deg_A);
    count[4] = count[4] + 56*tA3;
    count[4] = count[4]/14.0;
  }
  
  return count;
}
#include <RcppArmadillo.h>
using namespace Rcpp;

arma::vec count_k_cycle(arma::mat &A, int max_cycle_order);

// Rcpp implementation for counting k-cycles
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.net_summary_subsample_adj)]]
arma::mat net_summary_subsample_adj(arma::mat &A, arma::vec subsample_sizes, int max_cycle_order, int R) {
  int n = A.n_cols;
  int Ns = subsample_sizes.n_elem;
  int subsample_size;
  arma::mat result(R, max_cycle_order-1);
  arma::vec deg_A;

  arma::vec powers  = 1.0/arma::regspace(3, max_cycle_order);
  arma::uvec sample_ind;
  arma::uvec ind_cycles = arma::regspace<arma::uvec>(1, max_cycle_order-2);
  arma::vec numers;
  arma::vec denoms;
  
  for(int r = 0; r<R; r++){
    // subsample_size = subsample_sizes(ceil((double)Ns*r/R));
    // 
    // sample_ind =  <arma::uvec>sample(n, subsample_size, false);
    // deg_A = sum(A.submat(sample_ind,sample_ind), 1);
    // result(r,0) = (1/2*sum(deg_A*(deg_A-1.0))) / ((subsample_size-2)/2*sum(deg_A)); // (S.57) in Maugis et al. (2017)
    // //  Compute k-cycle ratio defined in (S.58)
    // numers = count_k_cycle(A.submat(sample_ind,sample_ind), max_cycle_order);
    // denoms = numers;
    // result(r, ind_cycles) = arma::pow(numers/denoms, powers);
  }
  
  return result;
}

// 
// net_summary_subsample_adj <- function(A, subsample_sizes, max_cycle_order, R){
//   n <- dim(A)[1]
//   Ns <- length(subsample_sizes)
//   
//   
//   result <- matrix(0, R, max_cycle_order-1)
//   colname_summary <- c("Trees", 'Triangles', 'Squares', "Pentagons", "Hexagons", "Septagons")
//   colnames(result) <- colname_summary[1:(max_cycle_order-1)]
//   
//   ind_cycles <- 2:(max_cycle_order-1)
//   powers <- 1/seq(3,max_cycle_order,1)
//   
//   for(r in 1:R){
//     subsample_size <- subsample_sizes[ceiling(Ns*r/R)]
//     denoms <- sapply(3:max_cycle_order, function(k, n) .ffct(subsample_size,k)/(2*k), n=subsample_size)
//     sample_ind <- sample.int(n, subsample_size)
//     deg_A <- colSums(A[sample_ind, sample_ind])
//     result[r,1] <- (1/2*sum(deg_A*(deg_A-1))) / ((subsample_size-2)/2*sum(deg_A))# (S.57) in Maugis et al. (2017)
// # Compute k-cycle ratio defined in (S.58)
//     numers <- .count_k_cycle(A[sample_ind, sample_ind], max_cycle_order)
//     result[r,ind_cycles] <- (numers/denoms)^powers
//   }
//   return(result)
// }

/*** R
timesTwo(42)
*/

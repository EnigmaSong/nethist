#include <RcppArmadillo.h>
#include "count_k_cycle.h"
#include "nethist_utils.h"

using namespace Rcpp;

// Rcpp implementation for counting k-cycles
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.net_summary_subsample_adj)]]
arma::mat net_summary_subsample_adj(const arma::mat &A, const arma::vec &subsample_sizes, int max_cycle_order, int R) {
  const int n = A.n_cols;
  const int Ns = subsample_sizes.n_elem;
  int subsample_size;  
  int temp_ind;
  arma::mat result(R, max_cycle_order-1);
  arma::vec deg_A_sub;

  const arma::vec powers = 1.0/arma::regspace<arma::vec>(3, max_cycle_order);
  const arma::uvec total_ind = arma::regspace<arma::uvec>(0, n-1);
  arma::uvec sample_ind;
  arma::vec cycle_ratio;
  arma::mat denoms(max_cycle_order-2,Ns);
  
  for(int k = 3; k <= max_cycle_order; k++){
    for(int j = 0; j < Ns; j++){
      denoms(k-3,j) = ffct((int)subsample_sizes(j), k)/(2.0*k);
    }
  }
  
  for(int r = 0; r<R; r++){
    temp_ind = (int)floor((double)Ns*r/R);
    subsample_size = subsample_sizes(temp_ind);

    sample_ind = sample(total_ind, (int)subsample_size, false);
    deg_A_sub = sum(A(sample_ind,sample_ind), 1);
    
    result(r,0) = (1/2.0*sum(deg_A_sub%(deg_A_sub-1.0))) / ((subsample_size-2)/2*sum(deg_A_sub)); // (S.57) in Maugis et al. (2017)
    //  Compute k-cycle ratio defined in (S.58)
    cycle_ratio = count_k_cycle(A(sample_ind,sample_ind), max_cycle_order);
    cycle_ratio = cycle_ratio/denoms.col(temp_ind);
    result(r, arma::span(1, max_cycle_order-2)) = arma::pow(cycle_ratio, powers).t();
  }
  
  return result;
}

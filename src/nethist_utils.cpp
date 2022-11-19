#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.ffct)]]
double ffct(int n, int k){
  double result = 1;
  for(int i=n; i > n-k; i--){
    result = result*(double)i;
  }
  
  return (double) result;
}

// [[Rcpp::export(.hamming_dist_adj_mat)]]
arma::mat hamming_dist_adj_mat(const arma::mat &A){
  //To compute hamming distance matrix for adjacency matrix from pairs of column vectors
  //It is faster than dist(method="manhattan"), which uses primitive C function, because 
  //i) it does not have error handling in the C primitive (= no if statement for error handling like dist())
  //ii) it runs smaller number of arithmetic operations.
  int n = A.n_cols;
  arma::mat result(n,n);
  
  for(int i = 0; i < n; i++){
    for(int j = i+1; j <n; j++){
      result(i,j) += sum(A.col(i)!= A.col(j));
    }
  }
  result = symmatu(result);
  
  return result;
}

arma::vec sample(const arma::vec& x, const int& size, const bool& replace){
  return Rcpp::RcppArmadillo::sample(x, size, replace);
}
arma::uvec sample(const arma::uvec& x, const int& size, const bool& replace){
  return Rcpp::RcppArmadillo::sample(x, size, replace);
}
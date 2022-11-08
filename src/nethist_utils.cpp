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

arma::vec sample(const arma::vec& x, const int& size, const bool& replace){
  return Rcpp::RcppArmadillo::sample(x, size, replace);
}
arma::uvec sample(const arma::uvec& x, const int& size, const bool& replace){
  return Rcpp::RcppArmadillo::sample(x, size, replace);
}
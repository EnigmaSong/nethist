#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include "multinethist.h"

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
arma::mat hamming_dist_adj_mat(const arma::imat &A){
  const arma::uword n = A.n_cols;
  arma::mat result(n,n);
  
  for(arma::uword j = 0; j < n; j++){
    for(arma::uword i = j+1; i < n; i++){
      result.at(i,j) += sum(A.col(i)!= A.col(j));
    }
  }
  result = symmatl(result);
  
  return result;
}

// [[Rcpp::export(.is_undirected_simple)]]
bool is_undirected_simple(const arma::mat& A){
  // Checking simple & undirected graph
  // If return is true, then assume A be imat or icube.
  const arma::uword n = A.n_cols;
  
  if(n != A.n_rows){
    warning("A is not a square matrix.\n"); 
    return false;
  } 
  
  for(arma::uword j = 0; j < n; j++){
    if(A.at(j,j)!=0){
      warning("A has self-loops.\n");
      return false;
    }
    for(arma::uword i = j+1; i < n; i++){
      if((A.at(i,j)!=0)&&(A.at(i,j)!=1)){
        warning("A is not simple.\n");
        return false;
      }
      if(A.at(i,j)!=A.at(j,i)){
        warning("A is not symmetric.\n");
        return false;
      }
    }
  }
  return true;
}

arma::vec sample(const arma::vec& x, const int& size, const bool& replace){
  return Rcpp::RcppArmadillo::sample(x, size, replace);
}
arma::uvec sample(const arma::uvec& x, const int& size, const bool& replace){
  return Rcpp::RcppArmadillo::sample(x, size, replace);
}

arma::uvec Cind_to_Rind(arma::uvec Cind){
  return Cind + 1;
}

arma::mat NegBinomEnt(const arma::mat& x){
  return x%log(x) + (1-x)%log(1-x);
}
arma::vec NegBinomEnt(const arma::vec& x){
  return x%log(x) + (1-x)%log(1-x);
}

double NegBinomEnt(const double& x){
  return x*log(x) + (1-x)*log(1-x);
}
double clamp_eps(const double &x){
  return x < eps ? eps : (x > 1 - eps ? 1 - eps : x);
}
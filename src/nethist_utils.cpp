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
  int n = A.n_cols;
  arma::mat result(n,n);
  
  for(int i = 0; i < n; i++){
    for(int j = i+1; j < n; j++){
      result(i,j) += sum(A.col(i)!= A.col(j));
    }
  }
  result = symmatu(result);
  
  return result;
}

// [[Rcpp::export(.is_undirected_simple)]]
bool is_undirected_simple(const arma::mat& A){
  // Checking simple & undirected graph
  int n = A.n_cols;
  
  if(n != A.n_rows){
    Rcout<< "A is not a square matrix.\n"; 
    return false;
  } 
  
  for(int i = 0; i < n; i++){
    for(int j = i+1; j <n; j++){
      if(((A.at(i,j)!=0)&&(A.at(i,j)!=1))||(A.at(i,j)!=A.at(j,i))){
        Rcout<< "A is not simple or symmetric.\n";
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
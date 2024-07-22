#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include "multinethist.h"

using namespace Rcpp;

arma::uvec select_swap(const Assignment& assignment, 
                       const int& n_node,
                       const int& swap_rule){
  //Current: only uniform random sample is implemented.
  const arma::uword selected_node1 = sample_node_by_swap_rule(n_node, swap_rule); 
  const int label1 = assignment.node_labels.at(selected_node1); 
  arma::uword selected_node2 = sample_node_by_swap_rule(n_node, swap_rule); 
  arma::uvec selected_nodes(2);
  
  //if nodes have the same label, re-draw node2. 
  int i = 0;
  while((label1 == assignment.node_labels.at(selected_node2)) && (i < 10)){
    selected_node2 = sample_node_by_swap_rule(n_node, swap_rule);
    i++;
  }
  
  selected_nodes.at(0) = selected_node1;
  selected_nodes.at(1) = selected_node2;
  
  return selected_nodes;
}
arma::uword sample_node_by_swap_rule(const int& n_node, const int& swap_rule){
  switch(swap_rule){ 
  case 1: // uniform random sample (equivalent to sample.int() in R)
    return as<arma::uword>(wrap(sample(n_node, 1, false) - 1));
  }
  return 0;
}

bool BestInfo::check_stop_rule_LL(const Assignment& current, 
                                  const int& consecutive_iter_threshold,
                                  const double& normalizeC, 
                                  const int& i){
  if(normalized_bestLL < current.likelihood*normalizeC){
    normalized_bestLL = current.likelihood*normalizeC;
    best_iter = i;
  }
  
  return (i > best_iter + consecutive_iter_threshold);
}

bool BestInfo::check_stop_rule_LL(const AssignSingleLayer& current, 
                                  const int& consecutive_iter_threshold,
                                  const double& normalizeC, 
                                  const int& i){
  if(normalized_bestLL < current.likelihood*normalizeC){
    normalized_bestLL = current.likelihood*normalizeC;
    best_iter = i;
  }
  
  return (i > best_iter + consecutive_iter_threshold);
}

bool BestInfo::check_stop_rule_LSE(const Assignment& current, 
                                   const int& consecutive_iter_threshold,
                                   const double& normalizeC, //dummy variable
                                   const int& i){
  if(bestLSE > current.LSE){
    bestLSE = current.LSE;
    best_iter = i;
  }
  
  return (i > best_iter + consecutive_iter_threshold);
}

bool BestInfo::check_stop_rule_LSE(const AssignSingleLayer& current, 
                                   const int& consecutive_iter_threshold,
                                   const double& normalizeC, //dummy variable
                                   const int& i){
  if(bestLSE > current.LSE){
    bestLSE = current.LSE;
    best_iter = i;
  }
  
  return (i > best_iter + consecutive_iter_threshold);
}

// micellenous

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
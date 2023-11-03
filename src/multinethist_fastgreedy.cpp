#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>
#include "multinethist.h"
#include "nethist_utils.h"

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.multinethist_fastgreedy)]]
List multinethist_fastgreedy(const arma::icube &A, const int &hbar, 
                             const arma::uvec &inputLabelVec, 
                             const int &max_itr,
                             const int &swap_rule,
                             const int &consecutive_iter_threshold,
                             const bool &verbose){
  const double normalizeC = 2.0/accu(A);
  const GroupSize group_size(A.n_rows, hbar);
  Assignment current(A, inputLabelVec, group_size);
  Assignment proposal = current;
  BestInfo best_iter_LL(current.likelihood * normalizeC);
  
  if(verbose) Rcout<< "Initial LL="<<current.likelihood*normalizeC<<"\n";
  
  for(int i = 1; i < max_itr; i++){
    proposal.create_proposal_from(current, A, swap_rule);
    if(proposal > current) current.copy_labels_theta_LL(proposal);
    proposal.copy_labels_theta_LL(current);
    
    if(best_iter_LL.check_stop_rule(current, consecutive_iter_threshold, normalizeC, i)){
      if(verbose) Rcout<<"total num iter = "<<i<<"\n";
      break;
    }
  }
  
  if(verbose) Rcout<< "best_iter at the end="<<best_iter_LL.best_iter<<" w/ LL="<<best_iter_LL.normalized_bestLL<<"\n";

  return List::create(Named("ThetaHat") = current.estimated_theta,
                      Named("node_labels") = Cind_to_Rind(current.node_labels),
                      Named("norm_LL") = current.likelihood*normalizeC
  );
}
// [[Rcpp::export(.mnhistCommon_fastgreedy)]]
List mnhistCommon_fastgreedy(const arma::icube &A, const int &hbar, 
                             const arma::uvec &inputLabelVec, 
                             const int &max_itr,
                             const int &swap_rule,
                             const int &consecutive_iter_threshold,
                             const bool &verbose){
  const double normalizeC = 2.0/accu(A);
  const GroupSize group_size(A.n_rows, hbar);
  AssignCommonF current(A, inputLabelVec, group_size);
  AssignCommonF proposal = current;
  BestInfo best_iter_LL;
  
  if(verbose) Rcout<< "Initial LL="<<current.likelihood*normalizeC<<"\n";
  
  for(int i = 1; i < max_itr; i++){
    proposal.create_proposal_from(current, A, swap_rule);
    if(proposal > current) current.copy_labels_theta_LL(proposal);
    proposal.copy_labels_theta_LL(current);
    
    if(best_iter_LL.check_stop_rule(current, consecutive_iter_threshold, normalizeC, i)){
      if(verbose) Rcout<<"total num iter = "<<i<<"\n";
      break;
    }
  }
  
  if(verbose) Rcout<< "best_iter at the end="<<best_iter_LL.best_iter<<" w/ LL="<<best_iter_LL.normalized_bestLL<<"\n";
  
  return List::create(Named("ThetaHat") = current.estimated_theta,
                      Named("node_labels") = Cind_to_Rind(current.node_labels),
                      Named("norm_LL") = current.likelihood*normalizeC
                      );
}

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

bool BestInfo::check_stop_rule(const Assignment& best, 
                              const int& stop_rule,
                              const double& normalizeC, 
                              const int& i){
  if(normalized_bestLL < best.likelihood*normalizeC){
    normalized_bestLL = best.likelihood*normalizeC;
    best_iter = i;
  }
  
  return (i > best_iter + stop_rule);
}


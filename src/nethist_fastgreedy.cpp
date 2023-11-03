#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>
#include "multinethist.h"
#include "nethist_utils.h"

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.nethist_fastgreedy)]]
List nethist_fastgreedy(const arma::imat &A, const int &hbar, 
                       const arma::uvec &inputLabelVec, 
                       const int &max_itr,
                       const int &swap_rule,
                       const int &consecutive_iter_threshold,
                       const bool &verbose){
  const double normalizeC = 2.0/accu(A);
  GroupSize group_size(A.n_rows, hbar);
  AssignSingleLayer current(A, inputLabelVec, group_size);
  AssignSingleLayer proposal = current;
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

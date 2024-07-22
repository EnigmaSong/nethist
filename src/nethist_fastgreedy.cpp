#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>
#include "multinethist.h"
#include "nethist_utils.h"

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.nethist_fastgreedy)]]
List nethist_fastgreedy(const arma::imat &A, const int &hbar, 
                       const arma::uvec &inputLabelVec, 
                       const int &method,
                       const int &max_itr,
                       const int &swap_rule,
                       const int &consecutive_iter_threshold,
                       const bool &verbose){
  const double normalizeC = 2.0/accu(A);
  GroupSize group_size(A.n_rows, hbar);
  AssignSingleLayer current(A, inputLabelVec, group_size);
  AssignSingleLayer proposal = current;
  BestInfo best_iter;
  
  // Select the appropriate operator
  Comparator comp;
  CheckStopRuleFunction check_stop_rule;
  
  if (method == 1) {
    comp = &AssignSingleLayer::operator>;
    check_stop_rule = &BestInfo::check_stop_rule_LL;
  } else if (method == 2) {
    comp = &AssignSingleLayer::operator<<;
    check_stop_rule = &BestInfo::check_stop_rule_LSE;
  } else {
    Rcpp::stop("Invalid method.");
  }
  
  if(verbose && (method == 1)) Rcout<< "Initial LL="<<current.likelihood*normalizeC<<"\n";
  if(verbose && (method == 2)) Rcout<< "Initial LSE="<<current.LSE<<"\n";
  
  
  for(int i = 1; i < max_itr; i++){
    proposal.create_proposal_from(current, A, swap_rule);
    if((proposal.*comp)(current)) current.copy_labels_theta(proposal);
    proposal.copy_labels_theta(current);
    
    if((best_iter.*check_stop_rule)(current, consecutive_iter_threshold, normalizeC, i)){
      Rcout << "Best LSE=" << current.LSE << ", proposed LSE = " << proposal.LSE <<"\n";
      if(verbose) Rcout<<"total num iter = "<<i<<"\n";
      break;
    }
  }
  
  if(verbose && (method == 1)) Rcout<< "best_iter at the end="<<best_iter.best_iter<<" w/ LL="<<best_iter.normalized_bestLL<<"\n";
  if(verbose && (method == 2)) Rcout<< "best_iter at the end="<<best_iter.best_iter<<" w/ LSE="<<best_iter.bestLSE<<"\n";

  return List::create(Named("ThetaHat") = current.estimated_theta,
                      Named("node_labels") = Cind_to_Rind(current.node_labels),
                      Named("norm_LL") = current.likelihood*normalizeC,
                      Named("LSE") = current.LSE
  );
}

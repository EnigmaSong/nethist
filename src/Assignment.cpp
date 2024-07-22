#include <RcppArmadillo.h>
#include "multinethist.h"
#include "nethist_utils.h"
using namespace Rcpp;

//Assignment class methods
Assignment::Assignment(const arma::icube& A, 
           const arma::uvec& v_labels,
           const GroupSize& g_size
) : n_nodes(A.n_rows), n_layers(A.n_slices),  node_labels(v_labels), 
group_size(g_size),
bin_cell_counts(arma::zeros(group_size.number_group, group_size.number_group)),
bin_edge_counts(arma::zeros(group_size.number_group, group_size.number_group, n_layers)),
estimated_theta(arma::zeros(group_size.number_group, group_size.number_group, n_layers)),
likelihood(0)
{
  arma::uvec node_labels_j;

  for(arma::uword l = 0; l < n_layers; l++){  
    for(arma::uword j = 0; j < group_size.number_group; j++){
      node_labels_j = find(node_labels==j);
      bin_edge_counts.at(j, j, l) = accu(A.slice(l)(node_labels_j,node_labels_j)) / 2.0;
      for(arma::uword i = j + 1; i < group_size.number_group; i++){
        bin_edge_counts.at(i,j,l) = accu(A.slice(l)(find(node_labels==i),node_labels_j));
        bin_edge_counts.at(j,i,l) = bin_edge_counts.at(i,j,l);
      }
    }
  }
  bin_cell_counts = group_size.group_number*group_size.group_number.t();
  bin_cell_counts.diag() = (group_size.group_number-1)%group_size.group_number/2;
  
  estimated_theta = bin_edge_counts.each_slice()/bin_cell_counts;
  likelihood = compute_normalize_log_likelihood();
  LSE = compute_LSE();
}

void Assignment::create_proposal_from(const Assignment& current, 
                           const arma::icube& A, 
                           const int& swap_rule){
  const arma::uvec swap = select_swap(current, (int)A.n_rows, swap_rule);

  update_thetahat(swap, A);
  swap_nodes(current, swap);
  updateLL();
  updateLSE();
}
double Assignment::compute_normalize_log_likelihood(){
  double loglik = 0.0;
  double theta_ijl;
  
  for(arma::uword l = 0; l < n_layers; l++){
    for(arma::uword j = 0; j < num_groups(); j++){
      for(arma::uword i = j; i< num_groups() ; i++){
        theta_ijl = clamp_eps(estimated_theta.at(i,j,l));
        loglik += (theta_ijl * log(theta_ijl) + (1 - theta_ijl)*log(1 - theta_ijl))*bin_cell_counts.at(i,j);
      }
    }
  }
  return loglik;
}
double Assignment::compute_LSE(){
  double LSE = 0.0;
  double theta_ijl;
  
  for(arma::uword l = 0; l < n_layers; l++){
    for(arma::uword j = 0; j < num_groups(); j++){
      for(arma::uword i = j; i< num_groups() ; i++){
        theta_ijl = clamp_eps(estimated_theta.at(i,j,l));
        LSE += (pow(theta_ijl - 1, 2) * bin_edge_counts.at(i,j,l)) + (pow(theta_ijl,2)*(bin_cell_counts.at(i,j) - bin_edge_counts.at(i,j,l)));
      }
    }
  }
  return LSE;
}
void Assignment::copy_labels_theta(const Assignment& other){
  //assume no change on bin size
  node_labels = other.node_labels;
  bin_edge_counts = other.bin_edge_counts;
  estimated_theta = other.estimated_theta;
  likelihood = other.likelihood;
  LSE = other.LSE;
}

//Update node lables by swap rule
void Assignment::swap_nodes(const Assignment& assignment, const arma::uvec& swap){
  node_labels.at(swap.at(0)) = assignment.node_labels.at(swap.at(1));
  node_labels.at(swap.at(1)) = assignment.node_labels.at(swap.at(0));
}
void Assignment::update_thetahat(const arma::uvec& swap, const arma::icube& A){
  //Assume node labels are not swapped yet.
  const arma::uword group_node1 = node_labels(swap.at(0));
  const arma::uword group_node2 = node_labels(swap.at(1));
  arma::uword group_i;

  for(arma::uword i = 0; i < n_nodes; i++){
    if(any(swap == i)){
      continue;
    }
    group_i = node_labels.at(i);
    for(arma::uword l = 0; l < n_layers; l++){  
      if(A.at(i, swap.at(0), l)==A.at(i, swap.at(1), l)){
        continue;
      }
      //Only one of A(i, swap.at(0), l) and A(i, swap.at(1), l) is 1.
      if(A.at(i, swap.at(0), l) == 1){
        bin_edge_counts.at(group_node1, group_i, l) -= 1;
        bin_edge_counts.at(group_i, group_node1, l) = bin_edge_counts.at(group_node1, group_i, l);
        
        bin_edge_counts.at(group_node2, group_i, l) += 1;
        bin_edge_counts.at(group_i, group_node2, l) = bin_edge_counts.at(group_node2, group_i, l);
        continue;
      }
      bin_edge_counts.at(group_node1, group_i, l) += 1;
      bin_edge_counts.at(group_i, group_node1, l) = bin_edge_counts.at(group_node1, group_i, l);
      
      bin_edge_counts.at(group_node2, group_i, l) -= 1;
      bin_edge_counts.at(group_i, group_node2, l) = bin_edge_counts.at(group_node2, group_i, l);
    }
  }

  estimated_theta = bin_edge_counts.each_slice()/bin_cell_counts;
}
void Assignment::updateLL(){
  likelihood = compute_normalize_log_likelihood();
}
void Assignment::updateLSE(){
  LSE = compute_LSE();
}

// AssignCommonF class methods
AssignCommonF::AssignCommonF(const arma::icube& A, 
                const arma::uvec& v_labels,
                const GroupSize& g_size
) : Assignment(A, v_labels, g_size),
    fhat_common(arma::zeros(group_size.number_group, group_size.number_group)),
    rho_hat(arma::zeros(n_layers))
{
  arma::uvec node_labels_j;
  
  for(arma::uword l = 0; l < n_layers; l++){
    rho_hat.at(l) = accu(A.slice(l))/((double)n_nodes*(n_nodes-1));
    for(arma::uword j = 0; j < group_size.number_group; j++){
      node_labels_j = find(node_labels==j);
      bin_edge_counts.at(j, j, l) = accu(A.slice(l)(node_labels_j,node_labels_j)) / 2.0;
      for(arma::uword i = j+1; i < group_size.number_group; i++){
        bin_edge_counts.at(i,j,l) = accu(A.slice(l)(find(node_labels==i), node_labels_j));
        bin_edge_counts.at(j,i,l) = bin_edge_counts.at(i,j,l);
      }
    }
  }
  bin_cell_counts = group_size.group_number*group_size.group_number.t();
  bin_cell_counts.diag() = (group_size.group_number-1)%group_size.group_number/2;
  
  get_thetahat_common_f();
  
  likelihood = compute_normalize_log_likelihood();
  LSE = compute_LSE();
}
void AssignCommonF::copy_labels_theta(const AssignCommonF& other){
  node_labels = other.node_labels;
  bin_edge_counts = other.bin_edge_counts;
  estimated_theta = other.estimated_theta;
  likelihood = other.likelihood;
  fhat_common = other.fhat_common;
  
}
void AssignCommonF::update_thetahat(const arma::uvec& swap, const arma::icube& A){
  const arma::uword group_node1 = node_labels(swap.at(0));
  const arma::uword group_node2 = node_labels(swap.at(1));
  arma::uword group_i;

  for(arma::uword i = 0; i < n_nodes; i++){
    if(any(swap == i)){
      continue;
    }
    group_i = node_labels.at(i);
    for(arma::uword l = 0; l < n_layers; l++){  
      if(A.at(i, swap.at(0), l)==A.at(i, swap.at(1), l)){
        continue;
      }
      //Only one of A(i, swap.at(0), l) and A(i, swap.at(1), l) is 1.
      if(A.at(i, swap.at(0), l) == 1){
        bin_edge_counts.at(group_node1, group_i, l) -= 1;
        bin_edge_counts.at(group_i, group_node1, l) = bin_edge_counts.at(group_node1, group_i, l);
        
        bin_edge_counts.at(group_node2, group_i, l) += 1;
        bin_edge_counts.at(group_i, group_node2, l) = bin_edge_counts.at(group_node2, group_i, l);
        continue;
      }
      bin_edge_counts.at(group_node1, group_i, l) += 1;
      bin_edge_counts.at(group_i, group_node1, l) = bin_edge_counts.at(group_node1, group_i, l);
      
      bin_edge_counts.at(group_node2, group_i, l) -= 1;
      bin_edge_counts.at(group_i, group_node2, l) = bin_edge_counts.at(group_node2, group_i, l);
    }
  }
  
  estimated_theta = bin_edge_counts.each_slice()/bin_cell_counts;
  get_thetahat_common_f();
  
}

void AssignCommonF::get_thetahat_common_f(){
  fhat_common = sum(estimated_theta, 2)/accu(rho_hat);// weighted average of \hat{f}_l with weight rho_hat(l)/sum(rho_hat).
  for(arma::uword l = 0; l < n_layers; l++){
    estimated_theta.slice(l) = fhat_common * rho_hat(l);
  }
}

// AssignSingleLayer class
AssignSingleLayer::AssignSingleLayer(const arma::icube& A, 
                                     const arma::uvec& v_labels,
                                     const GroupSize& g_size
) : AssignSingleLayer(A.slice(0), v_labels, g_size) {
  
}
AssignSingleLayer::AssignSingleLayer(const arma::imat& A, 
                    const arma::uvec& v_labels,
                    const GroupSize& g_size
) : Assignment(arma::zeros<arma::icube>(A.n_rows, A.n_cols, 1), v_labels, g_size),
    bin_edge_counts(arma::zeros(Assignment::num_groups(), Assignment::num_groups())),
    estimated_theta(arma::zeros(Assignment::num_groups(), Assignment::num_groups()))
{
  arma::uword label_node_j;
  
  for(arma::uword j = 0; j < n_nodes; j++){
    label_node_j = node_labels(j);
    for(arma::uword i = j+1; i < n_nodes; i++){
      bin_edge_counts.at(node_labels.at(i), label_node_j) += 1.0*A.at(i,j);
      bin_edge_counts.at(label_node_j, node_labels.at(i)) += 1.0*A.at(i,j);
    }
  }
  bin_edge_counts.diag() /= 2.0;
  
  bin_cell_counts = group_size.group_number*group_size.group_number.t();
  bin_cell_counts.diag() = (group_size.group_number-1)%group_size.group_number/2;
  
  estimated_theta = bin_edge_counts/bin_cell_counts;
  likelihood = compute_normalize_log_likelihood();
  LSE = compute_LSE();
}
void AssignSingleLayer::create_proposal_from(const AssignSingleLayer& current, 
                          const arma::imat& A, 
                          const int& swap_rule){
  const arma::uvec swap = select_swap(current, (int)A.n_rows, swap_rule);
  
  update_thetahat(swap, A);
  swap_nodes(current, swap);
  updateLL();
  updateLSE();
}
double AssignSingleLayer::compute_normalize_log_likelihood(){
  double loglik = 0.0;
  double theta_ij;
  
  for(arma::uword j= 0; j< num_groups(); j++){
    for(arma::uword i = j; i < num_groups(); i++){
      theta_ij = clamp_eps(estimated_theta.at(i,j));
      loglik += (theta_ij * log(theta_ij) + (1 - theta_ij)*log(1 - theta_ij))*bin_cell_counts.at(i,j);
    }
  }
  
  return loglik;
}
double AssignSingleLayer::compute_LSE(){
  double LSE = 0.0;
  double theta_ij;
  
  for(arma::uword j= 0; j< num_groups(); j++){
    for(arma::uword i = j; i < num_groups(); i++){
      theta_ij = clamp_eps(estimated_theta.at(i,j));
      LSE += (pow(theta_ij - 1, 2) * bin_edge_counts.at(i,j)) + (pow(theta_ij,2)*(bin_cell_counts.at(i,j) - bin_edge_counts.at(i,j)));
    }
  }
  
  return LSE;
}
void AssignSingleLayer::copy_labels_theta(const AssignSingleLayer& other){
  node_labels = other.node_labels;
  bin_edge_counts = other.bin_edge_counts;
  estimated_theta = other.estimated_theta;
  likelihood = other.likelihood;
  LSE = other.LSE;
}

void AssignSingleLayer::update_thetahat(const arma::uvec& swap, const arma::icube& A){
  update_thetahat(swap, A.slice(0));
}
void AssignSingleLayer::update_thetahat(const arma::uvec& swap, const arma::imat& A){
  //Assume node labels are not swapped yet.
  const arma::uword group_node1 = node_labels.at(swap.at(0));
  const arma::uword group_node2 = node_labels.at(swap.at(1));
  arma::uword group_i;
  
  for(arma::uword i = 0; i < n_nodes; i++){
    if(any(swap == i) || (A.at(i, swap.at(0))==A.at(i, swap.at(1)))){
      continue;
    }

    //Only one of A(i, swap.at(0)) and A(i, swap.at(1)) is 1.
    group_i = node_labels.at(i);
    if(A.at(i, swap.at(0)) == 1){
      bin_edge_counts.at(group_node1, group_i) -= 1;
      bin_edge_counts.at(group_i, group_node1) = bin_edge_counts.at(group_node1, group_i);

      bin_edge_counts.at(group_node2, group_i) += 1;
      bin_edge_counts.at(group_i, group_node2) = bin_edge_counts.at(group_node2, group_i);
      continue;
    }
    bin_edge_counts.at(group_node1, group_i) += 1;
    bin_edge_counts.at(group_i, group_node1) = bin_edge_counts.at(group_node1, group_i);

    bin_edge_counts.at(group_node2, group_i) -= 1;
    bin_edge_counts.at(group_i, group_node2) = bin_edge_counts.at(group_node2, group_i);
  }
  
  for(arma::uword i = 0; i < 2; i++){
    group_i = node_labels.at(swap.at(i));
    for(arma::uword j = 0; j < num_groups(); j++){
      estimated_theta.at(j, group_i) = bin_edge_counts.at(j, group_i)/bin_cell_counts(j, group_i);
      estimated_theta.at(group_i, j) = estimated_theta.at(j, group_i);
    }
  }
}

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
  const arma::uword g1 = node_labels(swap(0));
  const arma::uword g2 = node_labels(swap(1));

  const double old_contrib = compute_affected_LL(g1, g2);
  update_thetahat(swap, A);
  const double new_contrib = compute_affected_LL(g1, g2);

  likelihood = current.likelihood + (new_contrib - old_contrib);
  swap_nodes(current, swap);
  // LSE not used in multilayer stopping rule; skip updateLSE()
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
  
  //off-diagonal entries
  for(arma::uword l = 0; l < n_layers; l++){
    for(arma::uword j = 0; j < num_groups(); j++){
      for(arma::uword i = 0; i< num_groups() ; i++){
        theta_ijl = clamp_eps(estimated_theta.at(i,j,l));
        LSE += (pow(theta_ijl - 1, 2) * bin_edge_counts.at(i,j,l)) + (pow(theta_ijl,2)*(bin_cell_counts.at(i,j) - bin_edge_counts.at(i,j,l)));
      }
    }
  }
  LSE *= 2.0;
  //diagonal entries
  for(arma::uword l = 0; l < n_layers; l++){
    for(arma::uword i=0; i < num_groups()-1; i++){
      theta_ijl = clamp_eps(estimated_theta.at(i,i,l));
      LSE += pow(theta_ijl,2)*group_size.group_number(0);
    }
    theta_ijl = clamp_eps(estimated_theta.at(num_groups()-1,num_groups()-1,l));
    LSE += pow(theta_ijl,2)*group_size.group_number(1);
  }
  
  return LSE/(n_nodes*n_nodes);
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

  const arma::uword K = group_size.number_group;
  for(arma::uword l = 0; l < (arma::uword)n_layers; ++l){
    for(arma::uword j = 0; j < K; ++j){
      estimated_theta.at(group_node1, j, l) =
        bin_edge_counts.at(group_node1, j, l) / bin_cell_counts.at(group_node1, j);
      estimated_theta.at(j, group_node1, l) = estimated_theta.at(group_node1, j, l);
      estimated_theta.at(group_node2, j, l) =
        bin_edge_counts.at(group_node2, j, l) / bin_cell_counts.at(group_node2, j);
      estimated_theta.at(j, group_node2, l) = estimated_theta.at(group_node2, j, l);
    }
  }
}
void Assignment::updateLL(){
  likelihood = compute_normalize_log_likelihood();
}
void Assignment::updateLSE(){
  LSE = compute_LSE();
}
double Assignment::compute_affected_LL(arma::uword g1, arma::uword g2){
  const arma::uword K = group_size.number_group;
  double result = 0.0;
  for(arma::uword j = 0; j < K; ++j){
    for(arma::uword i = j; i < K; ++i){
      if(i != g1 && i != g2 && j != g1 && j != g2) continue;
      for(arma::uword l = 0; l < (arma::uword)n_layers; ++l){
        const double t = clamp_eps(estimated_theta.at(i, j, l));
        result += (t * log(t) + (1 - t) * log(1 - t)) * bin_cell_counts.at(i, j);
      }
    }
  }
  return result;
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
  
  get_thetahat_common_f_partial(group_node1, group_node2);

}

void AssignCommonF::get_thetahat_common_f(){
  fhat_common = sum(estimated_theta, 2)/accu(rho_hat);// weighted average of \hat{f}_l with weight rho_hat(l)/sum(rho_hat).
  for(arma::uword l = 0; l < n_layers; l++){
    estimated_theta.slice(l) = fhat_common * rho_hat(l);
  }
}
void AssignCommonF::get_thetahat_common_f_partial(arma::uword g1, arma::uword g2){
  const arma::uword K = group_size.number_group;
  const double rho_sum = accu(rho_hat);
  for(arma::uword j = 0; j < K; ++j){
    // Recompute fhat_common for pairs (g1,j) and (g2,j) from raw bin counts
    double raw_g1 = 0.0, raw_g2 = 0.0;
    for(arma::uword l = 0; l < (arma::uword)n_layers; ++l){
      raw_g1 += bin_edge_counts.at(g1, j, l) / bin_cell_counts.at(g1, j);
      raw_g2 += bin_edge_counts.at(g2, j, l) / bin_cell_counts.at(g2, j);
    }
    const double fc_g1 = raw_g1 / rho_sum;
    const double fc_g2 = raw_g2 / rho_sum;
    fhat_common.at(g1, j) = fc_g1;  fhat_common.at(j, g1) = fc_g1;
    fhat_common.at(g2, j) = fc_g2;  fhat_common.at(j, g2) = fc_g2;
    for(arma::uword l = 0; l < (arma::uword)n_layers; ++l){
      estimated_theta.at(g1, j, l) = fc_g1 * rho_hat(l);
      estimated_theta.at(j, g1, l) = fc_g1 * rho_hat(l);
      estimated_theta.at(g2, j, l) = fc_g2 * rho_hat(l);
      estimated_theta.at(j, g2, l) = fc_g2 * rho_hat(l);
    }
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
  
  //off-diagonal entries of (A - P_hat)^2
  for(arma::uword j= 0; j< num_groups(); j++){
    for(arma::uword i = j; i < num_groups(); i++){
      theta_ij = clamp_eps(estimated_theta.at(i,j));
      LSE += (pow(theta_ij - 1, 2) * bin_edge_counts.at(i,j)) + (pow(theta_ij,2)*(bin_cell_counts.at(i,j) - bin_edge_counts.at(i,j)));
    }
  }
  LSE *= 2.0;
  //diagonal entries
  for(arma::uword i=0; i < num_groups()-1; i++){
    theta_ij = clamp_eps(estimated_theta.at(i,i));
    LSE += pow(theta_ij,2)*group_size.group_number(0);
  }
  theta_ij = clamp_eps(estimated_theta.at(num_groups()-1,num_groups()-1));
  LSE += pow(theta_ij,2)*group_size.group_number(1);
  
  return LSE/(n_nodes*n_nodes);
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

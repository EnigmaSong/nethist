#ifndef GRAPHEST_FASTGREEDY_MULTINET
#define GRAPHEST_FASTGREEDY_MULTINET

//global constant
const double eps = arma::datum::eps;

//class declaration
class Assignment;
class AssignCommonF;
class AssignSingleLayer;
class GroupSize;

//Function declaration
//main function for greedy search
//multilayer (general)
Rcpp::List multinethist_fastgreedy(const arma::icube &A, const int &hbar, 
                               const arma::uvec &inputLabelVec, 
                               const int &max_itr,
                               const int &swap_rule,
                               const int &consecutive_iter_threshold,
                               const bool &verbose);
//Multilayer with common f
Rcpp::List mnhistCommon_fastgreedy(const arma::icube &A, const int &hbar, 
                             const arma::uvec &inputLabelVec, 
                             const int &max_itr,
                             const int &swap_rule,
                             const int &consecutive_iter_threshold,
                             const bool &verbose);
//Single Layer
Rcpp::List nethist_fastgreedy(const arma::imat &A, const int &hbar, 
                        const arma::uvec &inputLabelVec, 
                        const int &method,
                        const int &max_itr,
                        const int &swap_rule,
                        const int &consecutive_iter_threshold,
                        const bool &verbose);


arma::uvec select_swap(const Assignment& assignment, 
                      const int& n_node,
                      const int& swap_rule);
arma::uword sample_node_by_swap_rule(const int& n_node, const int& swap_rule);
double clamp_eps(const double &x);

//class declaration
//This version puts remainders in a smaller separate group. 
//Theory is based on remainders in the last group so that it is larger than others.
class GroupSize {
public:
  int number_group = 1;
  arma::vec group_number; //at(0), standard group size, at(1) remainder group size (if it exists)
  
  GroupSize(const int& number_nodes, const int& hbar
  ) : number_group(number_nodes/hbar), 
      group_number(arma::zeros(number_group)){
    
    group_number.fill(hbar);
    if(number_group*hbar < number_nodes){
      group_number.at(number_group-1) += number_nodes - (number_group)*hbar;
    }else if(number_group*hbar > number_nodes){
      Rcpp::stop("Error in group size setting: n < h*k. Check code implementation");
    }
  }
  
  int length(){
    return number_group;
  }
  bool operator==(const GroupSize& other) const{
    return (number_group == other.number_group) && all(group_number == other.group_number);
  }
};

//node label assignment info
class Assignment {
public:
  int n_nodes;
  int n_layers;
  arma::uvec node_labels;
  GroupSize group_size;

  arma::mat bin_cell_counts;
  arma::cube bin_edge_counts;
  arma::cube estimated_theta;
  double likelihood;
  double LSE;
  
  //constructor
  Assignment(const arma::icube& A, 
             const arma::uvec& v_labels,
             const GroupSize& g_size);
  //operator
  ////for likelihood
  virtual bool operator>(const Assignment& other) const {
    return likelihood > other.likelihood;
  }
  virtual bool operator<(const Assignment& other) const {
    return likelihood < other.likelihood;
  }
  ////for LSE
  virtual bool operator>>(const Assignment& other) const {
    return LSE > other.LSE;
  }
  virtual bool operator<<(const Assignment& other) const {
    return LSE < other.LSE;
  }
  bool operator==(const Assignment& other) const {
    bool res = (n_nodes == other.n_nodes) && (n_layers == other.n_layers);
    res = res && all(node_labels == other.node_labels) && (group_size == other.group_size);
    res = res && (likelihood == other.likelihood);
    res = res && (LSE == other.LSE);
    res = res && approx_equal(bin_cell_counts, other.bin_cell_counts, "absdiff", 2*eps);
    for(arma::uword l = 0; l < n_layers; l++){
      res = res && approx_equal(bin_edge_counts.slice(l), other.bin_edge_counts.slice(l), "absdiff", 2*eps);
      res = res && approx_equal(estimated_theta.slice(l), other.estimated_theta.slice(l), "absdiff", 2*eps);
    }
    return res;
  }
  
  //methods
  int num_groups(){
    return group_size.length();
  }
  virtual void create_proposal_from(const Assignment& current, 
                       const arma::icube& A, 
                       const int& swap_rule);
  virtual double compute_normalize_log_likelihood();
  virtual double compute_LSE();
  void copy_labels_theta(const Assignment& other);
  void swap_nodes(const Assignment& assignment, const arma::uvec& swap);
  virtual void update_thetahat(const arma::uvec& swap, const arma::icube& A);
  void updateLL();
  void updateLSE();
  double compute_affected_LL(arma::uword g1, arma::uword g2);
};
//derived class of Assignment when f are the same for all layers
class AssignCommonF : public Assignment {
public:
  arma::mat fhat_common;
  arma::vec rho_hat;
  
  //For multi-nethist
  AssignCommonF(const arma::icube& A, 
             const arma::uvec& v_labels,
             const GroupSize& g_size
  );
  //operator
  bool operator==(const AssignCommonF& other) const {
    bool res = (n_nodes == other.n_nodes) && (n_layers == other.n_layers);
    res = res && all(rho_hat == other.rho_hat);
    res = res && all(node_labels == other.node_labels) && (group_size == other.group_size);
    res = res && (likelihood == other.likelihood);
    res = res && (LSE == other.LSE);
    res = res && approx_equal(bin_cell_counts, other.bin_cell_counts, "absdiff", 2*eps);
    res = res && approx_equal(fhat_common, other.fhat_common, "absdiff", 2*eps);
    for(arma::uword l = 0; l < n_layers; l++){
      res = res && approx_equal(bin_edge_counts.slice(l), other.bin_edge_counts.slice(l), "absdiff", 2*eps);
      res = res && approx_equal(estimated_theta.slice(l), other.estimated_theta.slice(l), "absdiff", 2*eps);
    }
    return res;
  }
  //methods
  void copy_labels_theta(const AssignCommonF& other);
  void update_thetahat(const arma::uvec& swap, const arma::icube& A) override;
  void get_thetahat_common_f();
  void get_thetahat_common_f_partial(arma::uword g1, arma::uword g2);
};
//derived class of Assignment for single-layer network
class AssignSingleLayer : public Assignment {
public:
  arma::mat bin_edge_counts;
  arma::mat estimated_theta;
  
  //For multi-nethist
  AssignSingleLayer(const arma::icube& A, 
                  const arma::uvec& v_labels,
                  const GroupSize& g_size
  );
  AssignSingleLayer(const arma::imat& A, 
                      const arma::uvec& v_labels,
                      const GroupSize& g_size
  );
  //operator
  ////for likelihood
  bool operator>(const AssignSingleLayer& other) const {
    return likelihood > other.likelihood;
  }
  bool operator<(const AssignSingleLayer& other) const {
    return likelihood < other.likelihood;
  }
  ////for LSE
  bool operator>>(const AssignSingleLayer& other) const {
    return LSE > other.LSE;
  }
  bool operator<<(const AssignSingleLayer& other) const {
    return LSE < other.LSE;
  }
  bool operator==(const AssignSingleLayer& other) const{
    bool res = (n_nodes == other.n_nodes);
    res = res && all(node_labels == other.node_labels) && (group_size == other.group_size);
    res = res && (likelihood == other.likelihood);
    res = res && (LSE == other.LSE);
    res = res && approx_equal(bin_cell_counts, other.bin_cell_counts, "absdiff", 2*eps);
    res = res && approx_equal(bin_edge_counts, other.bin_edge_counts, "absdiff", 2*eps);
    res = res && approx_equal(estimated_theta, other.estimated_theta, "absdiff", 2*eps);

    return res;
  }
  //methods
  void create_proposal_from(const Assignment& current, 
                            const arma::icube& A, 
                            const int& swap_rule) override {};
  void create_proposal_from(const AssignSingleLayer& current, 
                                    const arma::imat& A, 
                                    const int& swap_rule);
  double compute_normalize_log_likelihood() override;
  double compute_LSE() override;
  void copy_labels_theta(const AssignSingleLayer& other);
  void update_thetahat(const arma::uvec& swap, const arma::icube& A) override;
  void update_thetahat(const arma::uvec& swap, const arma::imat& A);
};

struct BestInfo{
  int best_iter;
  double normalized_bestLL;
  double bestLSE;
  BestInfo() : best_iter(1), normalized_bestLL(-1e12), bestLSE(1e7){}
  BestInfo(const double& norm_LL, const double& LSE) : best_iter(1), normalized_bestLL(norm_LL), bestLSE(LSE){}
  bool check_stop_rule_LL(const Assignment& best, const int& consecutive_iter_threshold,
                       const double& normalizedC, 
                       const int& i);
  bool check_stop_rule_LSE(const Assignment& best, const int& consecutive_iter_threshold,
                           const double& normalizedC, 
                          const int& i);
  bool check_stop_rule_LL(const AssignSingleLayer& best, const int& consecutive_iter_threshold,
                          const double& normalizedC, 
                          const int& i);
  bool check_stop_rule_LSE(const AssignSingleLayer& best, const int& consecutive_iter_threshold,
                           const double& normalizedC, 
                           const int& i);
};

// Define function pointer type
typedef bool (AssignSingleLayer::*Comparator)(const AssignSingleLayer&) const;
typedef bool (BestInfo::*CheckStopRuleFunction)(const AssignSingleLayer&, const int&, 
                                      const double&, const int&);

#endif
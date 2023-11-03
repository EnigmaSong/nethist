#ifndef NETHIST_UTIL
#define NETHIST_UTIL

double ffct(int n, int k);
arma::mat hamming_dist_adj_mat(const arma::imat& A);
bool is_undirected_simple(const arma::imat& A);
arma::vec sample(const arma::vec& x, const int& size, const bool& replace);
arma::uvec sample(const arma::uvec& x, const int& size, const bool& replace);
arma::uvec Cind_to_Rind(arma::uvec Cind);
arma::mat NegBinomEnt(const arma::mat& x);
arma::vec NegBinomEnt(const arma::vec& x);
double NegBinomEnt(const double& x);
double clamp_eps(const double &x);

#endif
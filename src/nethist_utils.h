#ifndef NETHIST_UTIL
#define NETHIST_UTIL

double ffct(int n, int k);
arma::mat hamming_dist_adj_mat(const arma::mat& A);
bool is_undirected_simple(const arma::mat& A);
arma::vec sample(const arma::vec& x, const int& size, const bool& replace);
arma::uvec sample(const arma::uvec& x, const int& size, const bool& replace);


#endif
#ifndef NETHIST_UTIL
#define NETHIST_UTIL

#include <RcppArmadillo.h>
double ffct(int n, int k);
double hamming_dist_adj_mat(const arma::mat &A);
arma::vec sample(const arma::vec& x, const int& size, const bool& replace);
arma::uvec sample(const arma::uvec& x, const int& size, const bool& replace);


#endif
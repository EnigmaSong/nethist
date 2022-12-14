#ifndef GRAPHEST_FASTGREEDY
#define GRAPHEST_FASTGREEDY

arma::vec graphest_fastgreedy(const arma::mat &A, const int &hbar, const arma::vec bestLabelVec, const bool &verbose);

arma::mat getSampleCounts(const arma::mat &X, const arma::umat &clusterInds, const arma::vec &h);
double fastNormalizedBMLogLik(const arma::vec &thetaVec, const arma::vec &habSqrdVec, const double &sampleSize);
arma::vec rand_oneTwoVec(const int &n, const double &prob);

void clamp(double *x, const double &min_value, const double &max_value);
void check_bestcount_improvecount(int* bestCount, int* consecZeroImprovements);
double Delta_NegEnt(const arma::vec &habSqrdCola,
                    const arma::vec &habSqrdColb,
                    const double &habSqrdEntryab,
                    const arma::vec &thetaCola,
                    const arma::vec &thetaColb,
                    const double &thetaEntryab);
// double NegEnt(const double &x);
// arma::mat NegEnt(const arma::mat &X);
// double Delta_NegEnt(const arma::mat &habSqrdCols_ab,
//                       const double &habSqrdEntryab,
//                       const arma::mat &thetaCols_ab,
//                       const double &thetaEntryab,
//                       const arma::mat &oldThetaCols_ab,
//                       const double &oldThetaEntryab);
void update_tolCounter(const double &normalizedBestLL, const double &oldNormalizedBestLL, int *tolCounter);
bool check_quit_greedy(const int &tolCounter, 
                       const int &consecZeroImprovements, 
                       const int &tol_ZeroImprovements, const bool &verbose);

#endif
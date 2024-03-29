#include <RcppArmadillo.h>
#include "graphest_fastgreedy.h"
#include "nethist_utils.h"
#include <Rcpp/Benchmark/Timer.h>

using namespace Rcpp;
const double eps = arma::datum::eps;
const arma::uword IMPOSSIBLE_INDEX = 10000000;
const double absTol = 2.5*pow(10,-4);
// Rcpp implementation for graphest_fastgreedy(), 
//
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export(.graphest_fastgreedy)]]
arma::vec graphest_fastgreedy(const arma::mat &A, const int &hbar, arma::vec bestLabelVec, const bool &verbose){
  const int n = A.n_rows;
  const double sampleSize = (double)n*((double)n-1.0)/2.0;
  const double numOnes = accu(A);
  const double normalizeC = 2.0*sampleSize/numOnes;
  
  const int maxNumRestarts = 500;
  const bool allInds = (n<=256);
  const int numGreedySteps = (allInds ? (int)sampleSize: 2*pow(10.0,4));
  const bool smallerLastGroup = !(n%hbar == 0);
  const int k = (n%hbar == 0 ? n/hbar : n/hbar + 1);
  const int numEqualSizeGroup = k-smallerLastGroup;
  const int tol_ZeroImprovements = (allInds ? 2 :ceil((double)(k*n*(n-1))/(2.0*(double)numGreedySteps)));;
  
  //Variables used in greedy search algorithms
  //log-likelihood, label vectors, cluster index, ...
  arma::mat bestACounts(k,k);
  arma::umat bestClusterInds(hbar, k);
  double bestLL,oldbestLL;
  const arma::uvec aLeqb = trimatu_ind(arma::size(bestACounts));
  
  double oldNormalizedBestLL;
  int bestCount = 0;
  int consecZeroImprovements = 0;
  int tolCounter = 0;
  
  arma::mat trialACounts(k,k);
  arma::umat trialClusterInds(hbar, k);
  arma::vec trialLabelVec(n);
  double trialLL;
  
  Timer timer; 
  
  arma::vec oneTwoVec(numGreedySteps);
  arma::uvec iVec(numGreedySteps), jVec(numGreedySteps), kVec(numGreedySteps);
  const arma::uvec integerVec_nminusone = arma::regspace<arma::uvec>(0, n-1);
  
  int i,j,a,b;
  
  arma::vec h(k);
  arma::mat habSqrd(k,k);
  arma::vec habSqrdCola(k), habSqrdColb(k);
  double habSqrdEntryab;

  arma::vec oldThetaCola(k), oldThetaColb(k);
  double oldThetaEntryab;
  arma::vec thetaCola(k), thetaColb(k);
  double thetaEntryab;
  arma::vec AColiMinusColj(n);
  
  arma::vec sumAijc(numEqualSizeGroup);
  double sumAijEnd;
  
  double deltaNegEnt, oldDeltaNegEnt;
  double normalizedBestLL;
  
  //initialization
  //// Initialize cluster assignemnts in order
  h = get_h(hbar, n, k, smallerLastGroup, verbose);
  habSqrd = get_habSqrd(h, n);
  
  bestClusterInds = init_ClusterInds(bestLabelVec, h, k,
                                    numEqualSizeGroup,
                                    smallerLastGroup);
  bestACounts = getSampleCounts(A, bestClusterInds, h);
  bestLL = fastNormalizedBMLogLik(clamp(bestACounts(aLeqb)/habSqrd(aLeqb), eps, 1.0-eps), 
                                  habSqrd(aLeqb), sampleSize);
  oldNormalizedBestLL = bestLL*normalizeC;
  oldbestLL = bestLL;

  timer.step("start");
  
  for(int mm=1; mm <= maxNumRestarts; mm++){
    oneTwoVec = rand_oneTwoVec(numGreedySteps, 1.0/3.0); 
    iVec = sample(integerVec_nminusone, numGreedySteps, true);
    jVec = sample(integerVec_nminusone, numGreedySteps, true);
    kVec = sample(integerVec_nminusone, numGreedySteps, true);
    
    for(int m = 0; m < numGreedySteps; m++){
      trialACounts = bestACounts;
      trialClusterInds = bestClusterInds;
      trialLabelVec = bestLabelVec;
      trialLL = bestLL;
      for(int swapNum=1; swapNum <= oneTwoVec.at(m); swapNum++){
        if(swapNum==1){
          // Step 1 of 2
          i = iVec.at(m);
          j = jVec.at(m);
        }else{
          //Step 2 of 2
          i = jVec.at(m);
          j = kVec.at(m);
        }
        a = trialLabelVec.at(i);
        b = trialLabelVec.at(j);
        // Swap and update trial likelihood only if nodes i and j are in different clusters
        if(a == b) continue;
        
        trialLabelVec.at(i) = b;
        trialLabelVec.at(j) = a;
        
        habSqrdCola = habSqrd.col(a-1);
        habSqrdColb = habSqrd.col(b-1);
        habSqrdEntryab = habSqrd.at(a-1,b-1);

        oldThetaCola = trialACounts.col(a-1)/habSqrdCola;
        oldThetaColb = trialACounts.col(b-1)/habSqrdColb;
        oldThetaEntryab = trialACounts.at(a-1,b-1)/habSqrdEntryab;
        
        oldThetaCola.clamp(eps, 1.0-eps);
        oldThetaColb.clamp(eps, 1.0-eps);
        clamp(&oldThetaEntryab,eps, 1.0-eps);

        //Begin updating
        trialClusterInds.col(a-1).replace(i,j); // update that node j has replaced node i
        trialClusterInds.col(b-1).replace(j,i); // update that node i has replaced node j
        AColiMinusColj = A.col(i)-A.col(j);
        
        for(int kk = 0; kk< numEqualSizeGroup; kk++){
          sumAijc.at(kk) = sum(AColiMinusColj(trialClusterInds.col(kk)));
        }
        trialACounts(arma::span(0,numEqualSizeGroup-1), a-1) -= sumAijc;
        trialACounts(arma::span(0,numEqualSizeGroup-1), b-1) += sumAijc;
        
        if(smallerLastGroup){
          sumAijEnd = sum(AColiMinusColj(trialClusterInds(arma::span(0, h.at(k-1)-1),k-1)));
          trialACounts.at(k-1, a-1) -= sumAijEnd;
          trialACounts.at(k-1, b-1) += sumAijEnd;
        }
        // Update the above for special cases (c==a) or (c==b)
        trialACounts.at(a-1,a-1) += A.at(i,j);
        trialACounts.at(b-1,b-1) += A.at(i,j);
        if(smallerLastGroup && (b==k)){
          trialACounts.at(a-1,b-1) += -sum(AColiMinusColj(trialClusterInds(arma::span(0, h.at(k-1)-1),k-1))) - 2*A.at(i,j);
        }else{
          trialACounts.at(a-1,b-1) += -sum(AColiMinusColj(trialClusterInds.col(b-1))) - 2*A.at(i,j);
        }
        
        // Normalize and respect symmetry of trialAbar matrix
        trialACounts.row(b-1) = trialACounts.col(b-1).t();
        trialACounts.row(a-1) = trialACounts.col(a-1).t();
        
        // Now calculate changed likelihood directly
        thetaCola = trialACounts.col(a-1)/habSqrdCola;
        thetaColb = trialACounts.col(b-1)/habSqrdColb;
        thetaEntryab = trialACounts.at(a-1,b-1)/habSqrdEntryab;
        // Error handling to avoid p*log p and (1-p)*log(1-p) are NaN.
        thetaCola.clamp(eps, 1.0-eps);
        thetaColb.clamp(eps, 1.0-eps);
        clamp(&thetaEntryab, eps, 1.0-eps);
        
        // for this to work, we will have had to subtract out terms prior to updating
        deltaNegEnt = Delta_NegEnt(habSqrdCola,habSqrdColb,habSqrdEntryab,
                                   thetaCola,thetaColb,thetaEntryab);
        oldDeltaNegEnt = Delta_NegEnt(habSqrdCola,habSqrdColb,habSqrdEntryab,
                                      oldThetaCola,oldThetaColb,oldThetaEntryab);
        // Update log-likelihood - O(k)
        trialLL += (deltaNegEnt-oldDeltaNegEnt)/sampleSize;
      }
      // Metroplis or greedy step; if trial clustering accepted, then update current <-- trial
      if(trialLL > bestLL){
        bestACounts = trialACounts;
        bestClusterInds = trialClusterInds;
        bestLabelVec = trialLabelVec;
        bestLL = trialLL;
      }
    }
    // Keep track of best clustering overall
    if (bestLL > oldbestLL){// replace and save if trialLL is an improvement
      oldbestLL = bestLL;
      bestCount++;
    }
    // Keep track of best clustering overall
    if(mm%5==0){
      normalizedBestLL = bestLL*normalizeC;
      if(verbose) Rcout<< normalizedBestLL << " LL.  Iter " << mm<< " of max " << maxNumRestarts << "; "
                       << bestCount << " global improvements; \n";
      
      check_bestcount_improvecount(&bestCount, &consecZeroImprovements);
      update_tolCounter(normalizedBestLL, oldNormalizedBestLL, &tolCounter);
      oldNormalizedBestLL = normalizedBestLL;
      
      if(check_quit_greedy(tolCounter, consecZeroImprovements,
                           tol_ZeroImprovements, verbose)) break;
    }
  }
  timer.step("End");
  NumericVector timer_res(timer);
  if(verbose) Rcout<< "Greedy search: "<<(timer_res(1)/pow(10.0,9))<<" seconds\n";
  
  return bestLabelVec;
}

arma::vec get_h(const int &hbar, const int &n, 
                const int &k, const bool &smallerLastGroup,
                const bool &verbose){
  arma::vec h(k, arma::fill::value(hbar));
  
  if(verbose) Rcout << "Fitting a(n) " << k << " group block model\n";
  if(smallerLastGroup){
    h(k-1) = n%hbar; 
    if(verbose) Rcout << "Final group of size: " << h(k-1) << ", smaller than other of size: " << h(0) <<"\n";
  }else{
    if(verbose) Rcout << "All groups of equal size: " << h(1) <<"\n";
  }
  return h;
}

arma::mat get_habSqrd(const arma::vec &h, const int &n){
  arma::mat habSqrd(h.n_elem, h.n_elem);
  habSqrd = h*h.t()-diagmat(h%(h+1)/2); // (4.3) in Wolfe and Olhede (2013)
  // Error handlings
  if(any(any(habSqrd == 0))){
    stop("All clusters must contain at least 2 nodes.");
  }
  if(!(sum(h)==n)){
    stop("Number of cluster assignments must equal number of nodes.");
  }
  
  return habSqrd;
}
arma::mat getSampleCounts(const arma::mat &X, const arma::umat &clusterInds, const arma::vec &h){
  const int numClusters = clusterInds.n_cols;
  arma::mat Xsums(numClusters, numClusters);
  arma::uvec validIndCola, validIndColb;
  
  for(int a=0; a < numClusters; a++){
    validIndCola = clusterInds(arma::span(0,h.at(a)-1),a);
    Xsums.at(a,a) = accu(X.submat(validIndCola,validIndCola))/2.0;
    
    for(int b=a+1; b < numClusters; b++){
      validIndColb = clusterInds(arma::span(0,h.at(b)-1),b);
      
      Xsums.at(a,b) = accu(X.submat(validIndCola, validIndColb));
    }
  }
  Xsums = symmatu(Xsums);
  
  return Xsums;
}
//Convert bestLabelVec, a vector whose elements are between 1 to k, 
//to bestClusterInds, an (hbar x k) matrix whose elements are between 0 to n-1 and each column is group  
arma::umat init_ClusterInds(const arma::vec &LabelVec, 
                           const arma::vec &h, 
                           const int &k,
                           const int &numEqualSizeGroup,
                           const bool &smallerLastGroup){
  arma::umat ClusterInds(h.at(0), k, arma::fill::value(IMPOSSIBLE_INDEX));
  for(int a = 0; a < numEqualSizeGroup; a++){
    ClusterInds.col(a) = find(LabelVec == a + 1);
  }
  if(smallerLastGroup){
    ClusterInds(arma::span(0, h(k-1)-1), k-1) = find(LabelVec==k);
  }
  
  return ClusterInds;
}

double fastNormalizedBMLogLik(const arma::vec &thetaVec, const arma::vec &habSqrdVec, const double &sampleSize){
  return (dot(habSqrdVec,thetaVec%log(thetaVec) + (1.0-thetaVec)%log(1.0-thetaVec))/sampleSize);
}

void clamp(double *x, const double &min_value, const double &max_value){
  if(*x<=min_value){
    *x = min_value;
  }else if(*x>=max_value){
    *x = max_value;
  }
}

arma::vec rand_oneTwoVec(const int &n, const double &prob){
  return (arma::ones<arma::vec>(n) + (arma::randu<arma::vec>(n) >= prob));
}  


void check_bestcount_improvecount(int* bestCount, int* consecZeroImprovements){
  if(*bestCount == 0){
    *consecZeroImprovements = *consecZeroImprovements + 1;
  }else{
    *bestCount = 0;
    *consecZeroImprovements = 0;
  }
}

double Delta_NegEnt(const arma::vec &habSqrdCola,
                    const arma::vec &habSqrdColb,
                    const double &habSqrdEntryab,
                    const arma::vec &thetaCola,
                    const arma::vec &thetaColb,
                    const double &thetaEntryab){

  return(dot(habSqrdCola,(thetaCola%log(thetaCola) + (1.0-thetaCola)%log(1.0-thetaCola)))+
         dot(habSqrdColb,(thetaColb%log(thetaColb) + (1.0-thetaColb)%log(1.0-thetaColb)))-
         (habSqrdEntryab*(thetaEntryab*log(thetaEntryab) + (1.0-thetaEntryab)*log(1.0-thetaEntryab))));
}

void update_tolCounter(const double &normalizedBestLL, const double &oldNormalizedBestLL, int *tolCounter){
  *tolCounter = (normalizedBestLL - oldNormalizedBestLL < absTol ? *tolCounter+1 : 0);
}

bool check_quit_greedy(const int &tolCounter, 
                       const int &consecZeroImprovements, 
                       const int &tol_ZeroImprovements, const bool &verbose){
  if(tolCounter >=3){
    if(verbose) Rcout<<"3 consecutive likelihood improvements less than specified tolerance; quitting now\n";
    return true;
  }
  if(consecZeroImprovements == tol_ZeroImprovements){
    if(verbose) Rcout<<"Local optimum likely reached in random-ordered greedy likelihood search; quitting now\n";
    return true;
  }
  
  return false;
}

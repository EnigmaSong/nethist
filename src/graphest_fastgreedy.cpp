#include <RcppArmadillo.h>
#include "nethist_utils.h"
#include <Rcpp/Benchmark/Timer.h>
using namespace Rcpp;
const double eps = arma::datum::eps;
const arma::uword IMPOSSIBLE_INDEX = 10000000;
// Rcpp implementation for graphest_fastgreedy(), 
//
// [[Rcpp::depends(RcppArmadillo)]]

arma::mat getSampleCounts(arma::mat X, arma::umat clusterInds){
  int numClusters = size(clusterInds,1);
  int n = size(X,0);
  arma::mat Xsums(numClusters, numClusters);
  arma::uword validIndCola_Max, validIndColb_Max;
  arma::uvec validIndCola, validIndColb;
  
  for(int b=1; b<numClusters; b++){
    for(int a=0; a < b; a++){
      validIndCola_Max = find(clusterInds.col(a)>=0 && clusterInds.col(a)<n).max();
      validIndColb_Max = find(clusterInds.col(b)>=0 && clusterInds.col(b)<n).max();
      
      validIndCola = clusterInds.col(a).rows(0,validIndCola_Max);
      validIndColb = clusterInds.col(b).rows(0,validIndColb_Max);
      
      Xsums(a,b) = accu(X.submat(validIndCola, validIndColb));
    }
  }
  Xsums = Xsums + Xsums.t();
  for(int a=0; a< numClusters; a++){
    validIndCola_Max = find(clusterInds.col(a)>=0 && clusterInds.col(a)<n).max();
    validIndCola = clusterInds.col(a).rows(0,validIndCola_Max);
    Xsums(a,a) = accu(X.submat(validIndCola,validIndCola))/2;
  }
  return Xsums;
}

double fastNormalizedBMLogLik(arma::vec thetaVec, arma::vec habSqrdVec, double sampleSize){
  arma::vec negEntVec;
  double normLogLik;
  
  thetaVec(find(thetaVec<=0)).fill(eps);
  thetaVec(find(thetaVec>=1)).fill(1-eps);
  
  negEntVec = thetaVec%log(thetaVec) + (1-thetaVec)%log(1-thetaVec);
  normLogLik = sum(habSqrdVec%negEntVec)/sampleSize;
  
  return normLogLik;
}

// [[Rcpp::export(.graphest_fastgreedy)]]
arma::vec graphest_fastgreedy(arma::mat A, int hbar, arma::vec inputLabelVec, bool verbose){
  int n = A.n_rows;
  double n_float = (double)n;
  double sampleSize = n_float*(n_float-1)/2;
  
  double absTol = 2.5*pow(10,-4);
  const int maxNumRestarts = 500;
  bool allInds;
  int numGreedySteps;
  const bool smallerLastGroup = !(n%hbar == 0);
  int k;
  if(n%hbar == 0){
    k = n/hbar;
  }else{
    k = n/hbar + 1;
  }
  
  const int numEqualSizeGroup = k-smallerLastGroup;
  arma::vec h(k);
  
  arma::vec equalSizeInds = arma::regspace(0, numEqualSizeGroup - 1);
  arma::mat orderedClusterInds(hbar,k);
  arma::mat habSqrd(k,k);
  
  if(n<= 256){
    allInds = true;
  }else{
    allInds = false;
  }
  if(allInds){
    numGreedySteps = (int)sampleSize;
  }else{
    numGreedySteps = 2*pow(10.0,4);
  }
  // Initialize cluster assignemnts in order
  if(verbose) Rcout << "Fitting a(n) " << k << " group block model\n";
  
  h.fill(hbar);
  orderedClusterInds = arma::regspace(1, n);
  orderedClusterInds.reshape(hbar,k);
  orderedClusterInds= orderedClusterInds.t();
  
  if(smallerLastGroup){
    h(k-1) = n%hbar; 
    if(verbose) Rcout << "Final group of size: " << h(k-1) << ", smaller than other of size: " << h(0) <<"\n";
  }else{
    if(verbose) Rcout << "All groups of equal size: " << h(1) <<"\n";
  }
  habSqrd = h*h.t()-diagmat(h%(h-1));
  // Error handlings
  if(any(any(habSqrd == 0))){
    stop("All clusters must contain at least 2 nodes.");
  }
  if(orderedClusterInds.max() != n){
    stop("All nodes must be assigned to a cluster.");
  }
  if(!(sum(h)==n)){
    stop("Number of cluster assignments must equal number of nodes.");
  }
  
  arma::mat initialClusterCentroids(n,k);
  initialClusterCentroids.zeros();
  
  double bestLL;
  arma::mat bestACounts(k,k);
  arma::vec bestLabelVec(n);
  arma::umat bestClusterInds(hbar, k);
  arma::uvec aLeqb = trimatu_ind(size(bestACounts));
  
  bestLabelVec= inputLabelVec;
  bestClusterInds.fill(IMPOSSIBLE_INDEX);
  for(int a = 0; a < numEqualSizeGroup; a++){
    bestClusterInds.col(a) = find(bestLabelVec == a + 1);
    initialClusterCentroids.col(a) = sum(A.cols(bestClusterInds.col(a)), 1);
  }
  if(smallerLastGroup){
    bestClusterInds.col(k-1).rows(0, h(k-1)-1) = find(bestLabelVec==k);
  }
  bestACounts = getSampleCounts(A, bestClusterInds);
  bestLL = fastNormalizedBMLogLik(bestACounts(aLeqb)/habSqrd(aLeqb), habSqrd(aLeqb), sampleSize);
  
  double numOnes = accu(A);
  double oldNormalizedBestLL = bestLL*2*sampleSize/numOnes;
  int bestCount = 0;
  int consecZeroImprovements = 0;
  int tolCounter = 0;

  Timer timer; 
  
  arma::vec oneTwoVec(numGreedySteps);
  arma::vec iVec(numGreedySteps), jVec(numGreedySteps), kVec(numGreedySteps);
  arma::vec integerVec_nminusone = arma::regspace(0, n-1);
  
  arma::mat currentACounts(k,k);
  arma::umat currentClusterInds(hbar, k);
  arma::vec currentLabelVec(n);
  double currentLL;
  
  arma::mat trialACounts(k,k);
  arma::umat trialClusterInds(hbar, k);
  arma::vec trialLabelVec(n);
  double trialLL;
  
  int i,j,a,b;
  
  arma::vec habSqrdCola(k), habSqrdColb(k);
  arma::vec oldThetaCola(k), oldThetaColb(k);
  double habSqrdEntryab, oldThetaEntryab;
  arma::vec AColiMinusColj(n);
  arma::umat clusterIndMat(hbar, k-1);
  
  arma::vec sumAijc(numEqualSizeGroup);
  double sumAijEnd;
  
  arma::vec thetaCola, thetaColb;
  double thetaEntryab;
  double deltaNegEnt, oldDeltaNegEnt;
  double normalizedBestLL;
  
  timer.step("start");
  
  for(int mm=0; mm< maxNumRestarts; mm++){
    oneTwoVec = arma::ones<arma::vec>(numGreedySteps) + (arma::randu<arma::vec>(numGreedySteps) >= 1.0/3.0);
    iVec = sample(integerVec_nminusone, numGreedySteps, true);
    jVec = sample(integerVec_nminusone, numGreedySteps, true);
    kVec = sample(integerVec_nminusone, numGreedySteps, true);
    
    // Update Partition information
    currentACounts = bestACounts;
    currentClusterInds = bestClusterInds;
    currentLabelVec = bestLabelVec;
    currentLL = bestLL;
    for(int m = 0; m < numGreedySteps; m++){
      // Prepare to update quantities for trial clustering
      trialACounts = currentACounts;
      trialClusterInds = currentClusterInds;
      trialLabelVec = currentLabelVec;
      trialLL = currentLL;
      
      for(int swapNum=1; swapNum <= oneTwoVec(m); swapNum++){
        if(swapNum==1){
          // Step 1 of 2
          i = iVec(m);
          j = jVec(m);
          a = trialLabelVec(i);
          b = trialLabelVec(j);
        }else{
          //Step 2 of 2
          i = jVec(m);
          j = kVec(m);
          a = trialLabelVec(i);
          b = trialLabelVec(j);
        }
        // Swap and update trial likelihood only if nodes i and j are in different clusters
        if(a != b){
          trialLabelVec(i) = b;
          trialLabelVec(j) = a;
          
          habSqrdCola = habSqrd.col(a-1);
          habSqrdColb = habSqrd.col(b-1);
          habSqrdEntryab = habSqrd(a-1,b-1);
          
          oldThetaCola = trialACounts.col(a-1)/habSqrdCola;
          oldThetaColb = trialACounts.col(b-1)/habSqrdColb;
          oldThetaEntryab = trialACounts(a-1,b-1)/habSqrdEntryab;
          
          oldThetaCola(find(oldThetaCola<=0)).fill(eps);
          oldThetaCola(find(oldThetaCola>=1)).fill(1-eps);
          
          oldThetaColb(find(oldThetaColb<=0)).fill(eps);
          oldThetaColb(find(oldThetaColb>=1)).fill(1-eps);
          
          if(oldThetaEntryab <= 0){
            oldThetaEntryab = eps;
          }else if(oldThetaEntryab >= 1){
            oldThetaEntryab = 1-eps;
          }
          //Begin updating
          trialClusterInds.col(a-1).replace(i,j);
          trialClusterInds.col(b-1).replace(j,i);
          AColiMinusColj = A.col(i)-A.col(j);
          clusterIndMat = trialClusterInds.cols(0, numEqualSizeGroup-1);
          for(int kk = 0; kk< numEqualSizeGroup; kk++){
            sumAijc(kk) = sum(AColiMinusColj(clusterIndMat.col(kk)));
          }
          trialACounts.col(a-1).rows(0, numEqualSizeGroup-1) -= sumAijc;
          trialACounts.col(b-1).rows(0, numEqualSizeGroup-1) += sumAijc;
          if(smallerLastGroup){
            sumAijEnd = sum(AColiMinusColj(find(trialClusterInds.col(k-1)<k)));//Issue
            trialACounts(k-1, a-1) -= sumAijEnd;
            trialACounts(k-1, b-1) += sumAijEnd;
          }
          // Update the above for special cases (c==a) or (c==b)
          trialACounts(a-1,a-1) += A(i,j);
          trialACounts(b-1,b-1) += A(i,j);
          if(smallerLastGroup && (b==k)){
            trialACounts(a-1,b-1) += -sum(AColiMinusColj(find(trialClusterInds.col(b-1)))) - 2*A(i,j);//Issue
          }else{
            trialACounts(a-1,b-1) += -sum(AColiMinusColj(trialClusterInds.col(b-1))) - 2*A(i,j);
          }
          
          trialACounts(b-1,a-1) = trialACounts(a-1,b-1);
          // Normalize and respect symmetry of trialAbar matrix
          trialACounts.row(a-1) = trans(trialACounts.col(a-1));
          trialACounts.row(b-1) = trans(trialACounts.col(b-1));
          // Now calculate changed likelihood directly
          thetaCola = trialACounts.col(a-1)/habSqrdCola;
          thetaColb = trialACounts.col(b-1)/habSqrdColb;
          thetaEntryab = trialACounts(a-1,b-1)/habSqrdEntryab;
          // Error handling to avoid p*log p and (1-p)*log(1-p) are NaN.
          thetaCola(find(thetaCola <= 0)).fill(eps);
          thetaCola(find(thetaCola >= 1)).fill(1-eps);
          
          thetaColb(find(thetaColb <= 0)).fill(eps);
          thetaColb(find(thetaColb >= 1)).fill(1-eps);
          
          if(thetaEntryab <= 0){
            thetaEntryab = eps;
          }else if(thetaEntryab >= 1){
            thetaEntryab = 1-eps;
          }
          // for this to work, we will have had to subtract out terms prior to updating
          deltaNegEnt = sum(habSqrdCola%(thetaCola%log(thetaCola) + (1-thetaCola)%log(1-thetaCola)))+
            sum(habSqrdColb%(thetaColb%log(thetaColb) + (1-thetaColb)%log(1-thetaColb)))-
            (habSqrdEntryab*(thetaEntryab*log(thetaEntryab) + (1-thetaEntryab)*log(1-thetaEntryab)));
          oldDeltaNegEnt = sum(habSqrdCola%(oldThetaCola%log(oldThetaCola) + (1-oldThetaCola)%log(1-oldThetaCola)))+
            sum(habSqrdColb%(oldThetaColb%log(oldThetaColb) + (1-oldThetaColb)%log(1-oldThetaColb)))-
            (habSqrdEntryab*(oldThetaEntryab*log(oldThetaEntryab) + (1-oldThetaEntryab)*log(1-oldThetaEntryab)));
          // Update log-likelihood - O(k)
          trialLL += (deltaNegEnt-oldDeltaNegEnt)/sampleSize;
        }
      }
      // Metroplis or greedy step; if trial clustering accepted, then update current <-- trial
      if(trialLL > currentLL){
        currentACounts = trialACounts;
        currentClusterInds = trialClusterInds;
        currentLabelVec = trialLabelVec;
        currentLL = trialLL;
      }
    }
    // Keep track of best clustering overall
    if (currentLL > bestLL){// replace and save if trialLL is an improvement
      bestLL = currentLL; // update globally best visited likelihood
      bestACounts = currentACounts;
      bestClusterInds = currentClusterInds;
      bestLabelVec = currentLabelVec;
      bestCount++;
    }
    // Keep track of best clustering overall
    if(mm%5==0){
      // timer.step();
      // tElapsedOuter = difftime(Sys.time(), tStartOuter)
      normalizedBestLL = bestLL*2*sampleSize/numOnes;
      if(verbose) Rcout<< normalizedBestLL << " LL.  Iter " << mm << " of max " << maxNumRestarts << "; "
           << bestCount << " global improvements;\n";
      // << " took " << tElapsedOuter << " s\n";
      // message(paste0(round(normalizedBestLL,4), " LL.  Iter ", mm, " of max ", maxNumRestarts, "; ",
      //                bestCount, " global improvements; took ",
      //                round(tElapsedOuter,4), " s"))
      // tStartOuter = Sys.time()
      
      if(bestCount == 0){
        consecZeroImprovements++;
      }else{
        bestCount = 0;
        consecZeroImprovements = 0;
      }
      if(normalizedBestLL - oldNormalizedBestLL < absTol){
        tolCounter++;
      }else{
        tolCounter = 0;
      }
      oldNormalizedBestLL = normalizedBestLL;
      if(tolCounter >=3){
        if(verbose) Rcout<<"3 consecutive likelihood improvements less than specified tolerance; quitting now\n";
        break;
      }
      if(allInds == 1){
        if(consecZeroImprovements == 2){
          if(verbose) Rcout<<"Local optimum likely reached in random-ordered greedy likelihood search; quitting now\n";
          break;
        }
      }else{
        if(consecZeroImprovements == ceil((double)(k*n*(n-1))/(2.0*(double)numGreedySteps))){
          if(verbose) Rcout<<"Local optimum likely reached in random-ordered greedy likelihood search; quitting now\n";
          break;
        }
      }
    }
  }
  timer.step("End");
  NumericVector timer_res(timer);
  if(verbose) Rcout<< "Greedy search: "<<(timer_res(1)/pow(10.0,9))<<" seconds\n";
  
  return bestLabelVec;
}
// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/

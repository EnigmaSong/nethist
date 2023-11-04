##' Network histogram estimation
##'
##' Estimating network histogram for multiplex networks and returning the indices of partitions.
##'
##' @param A An adjacency array or list of igraph object. It must be an undirected and simple graph.
##' @param h A bandwidth parameter. If `NA`, selecting bandwidth by Olhede and Wolfe (2014). If specified, the user input value is used.
##' @param common_f a logical variable indicating assume the common network histogram function for all layers.
##' @param max_itr an integer for number of maximum iteration for greedy search.
##' @param swap_rule string of swap node selection methods. only "random" implemented. See details.
##' @param consecutive_iter_threshold an integer for stopping criterion. If the log-likelihood does not improve for the last `consecutive_iter_threshold` iterations, stop the algorithm.
##' @param verbose logical value indicating whether verbose output is generated.
##' @returns 
##' If number of layer is greater than 1, it returns an object of class ``multinethist``:
##' 
##' \itemize{
##' \item `cluster` a vector of partition indices.
##' \item `thetahat` a probability array from multinetwork histogram ordered by cluster labels. 
##' \item `rho_hat` a vector of estimated sparsity parameters. 
##' \item `normalized_LL` a normalized likelihood from the algorithm.
##' \item `homogeneous` a logical variable indicating homogeneous multinetwork histogram.
##' }
##' @details {
##' lth layer's multi-network histogram is defined by thetahat/rho_hat. We can plot multinetwork histogram using [plot()] and [plot3d()].
##' 
##' If number of layer is 1, then it calls single-layer network histogram. See. nethist() in nethist package.
##' 
##' Note that `cluster` only shows a partition of vertices, and the index labels is not an ordered variable. For example, nodes in cluster 1 do not have to more similar to nodes in cluster 2 than nodes in cluster 10. Hence, users would use a user-specified order in [plot.multinethist()].
##' }
##' @seealso [plot.multinethist()], [plot.nethist()]
##' @references Song, Y. & Olhede. S. C. (2023)
##' @import Rcpp
##' @importFrom stats .lm.fit dist pnorm weighted.mean
##' @importFrom RSpectra eigs
##' @export

multinethist <- function(A, h = NA, common_f = FALSE, 
                         max_itr = 5e6,
                         swap_rule = "random", 
                         consecutive_iter_threshold = 2e4,
                         verbose = FALSE){
  UseMethod("multinethist")
} 

##' @exportS3Method
multinethist.matrix <-  function(A, h = NA, common_f = FALSE, 
                                 max_itr = 5e6,
                                 swap_rule = "random", 
                                 consecutive_iter_threshold = 2e4,
                                 verbose = FALSE){
  return(multinethist.array(array(A, dim=c(nrow(A), ncol(A), 1)), 
                            h, common_f,
                            max_itr, swap_rule, 
                            consecutive_iter_threshold, verbose)) #should think about the design of code.
}
##' @exportS3Method
multinethist.array <- function(A, h = NA, common_f = FALSE, 
                         max_itr = 5e6,
                         swap_rule = "random", 
                         consecutive_iter_threshold = 2e4,
                         verbose = FALSE){
  if(!all(apply(A, 3, .is_undirected_simple))) stop("Network A must be an undirected simple network.")
  n_nodes <- dim(A)[1]
  n_layers <- dim(A)[3]
  swap_rule <- pmatch(swap_rule, c("random"))
  if(is.na(swap_rule)) stop("swap_rule must be one of the followings: random")
  
  # Compute necessary summaries from A
  rhoHat <- apply(A, 3, function(A_l) sum(A_l)/(n_nodes*(n_nodes-1)))
  
  ##########################################################################
  # Pick an analysis bandwidth and initialize via regularized spectral clustering
  # Temporary (use the densest layer for h)
  ##########################################################################
  if(is.na(h)){
    h <- .oracbwplugin(A, min(4, sqrt(n_nodes)/8), 
                       'degs', 1, rhoHat, common_f, verbose)$h
    if(verbose) message(paste("Determining bandwidth from data:", round(h)))
  }else{
    if(verbose) message(paste("Determining bandwidth from user input:", round(h)))
  }
  
  h <- max(2, min(n_nodes, round(h)))
  if(verbose){
    message(paste("Final bandwidth:", h))
    message(paste0('Adjacency matrix has ', n_nodes, ' rows/cols'))
  }
  
  # Initialize using regularized spectral clustering based on row similarity
  tstart <- Sys.time()
  idxInit <- initialize_index(A[,,which.max(rhoHat)], n_nodes, h, verbose)
  k <- ceiling(n_nodes/h)
  
  if(verbose) message(paste0('Initial label vector assigned from row-similarity ordering; time ',
                             round(difftime(Sys.time(),tstart),4), ' sec'))
  
  if(n_layers == 1){
    res<- .nethist_fastgreedy(A[,,1], h, Rind_to_Cind(idxInit), max_itr,
                              swap_rule, consecutive_iter_threshold, verbose)
  }else if(common_f){
    res<- .mnhistCommon_fastgreedy(A, h, Rind_to_Cind(idxInit), 
                                   max_itr, swap_rule, consecutive_iter_threshold, verbose)
  }else{
    res <- .multinethist_fastgreedy(A, h, Rind_to_Cind(idxInit), 
                                    max_itr, swap_rule, consecutive_iter_threshold, verbose)
  }
  
  
  result <- list(cluster = as.vector(res$node_labels), 
                 thetahat =  res$ThetaHat,
                 rho_hat = rhoHat,
                 normalized_LL = res$norm_LL,
                 homogeneous = common_f)
  result <- structure(result, class= ifelse(n_layers > 1, "multinethist", "nethist"))
  return(result)
}

.oracbwplugin <- function(A,c,type, alpha,
                          rhoHat, common_f, 
                          verbose){
  #Assume A[,,l] is symmetric, simple, and no self-loop
  if(missing(type)) type <- 'degs'
  if(!(type %in% c("degs","eigs"))) stop(paste("Invalid input type",type))
  if(missing(alpha)) alpha <- 1
  if(alpha != 1) stop("Currently only supports alpha = 1")
  
  n <- dim(A)[1]
  L <- ifelse(length(dim(A)) == 2, 1, dim(A)[3])
  
  midPt <- seq(round(n/2-c*sqrt(n),0), round(n/2+c*sqrt(n),0))
  rhoHat_inv <- ginv(rhoHat)
  sampleSize <- n*(n-1)/2
  estMSqrd <- rep(0,L)
  
  #Rank-1 graphon estimate via fhat(x,y) = mult*u(x)*u(y)*pinv(rhoHat);
  for(l in 1:L){
    estMSqrd[l] <- estim_M(A[,,l], type, 
                           MoreArgs = list(n = n, midPt = midPt, 
                                           rhoHat_inv = rhoHat_inv[l]))
  }

  mean_estMSqrd <- mean(estMSqrd)
  mean_edgenum <- mean(sampleSize*rhoHat)  
  weigths <- rhoHat/sum(rhoHat)
  wmean_estMSqrd <- weighted.mean(estMSqrd, weigths)
  mean_inv_edgenum <- mean(1/(sampleSize*rhoHat))
  h <- ifelse(common_f,
              sqrt(n)*(2*mean_estMSqrd^2*sum(rhoHat))^(-0.25),
              sqrt(n)*(2*mean(estMSqrd*rhoHat))^(-0.25)
  )
  
  MISEfhatBnd <- mean_estMSqrd*((2/sqrt(mean_estMSqrd))*sqrt(mean_inv_edgenum) + 1/n)
  WMISEfhatBnd <- wmean_estMSqrd*((2/sqrt(wmean_estMSqrd))*(mean_edgenum)^(-1/2) + 1/n)
  
  if(verbose){
    message("M^2_hat=")
    message(round(estMSqrd,3))
    message(ifelse(common_f, 
                   paste("MISE bound_hat=", round(WMISEfhatBnd,3)),
                   paste("WMISE bound_hat=", round(WMISEfhatBnd,3)))
            )
  }
  
  return(list(h=h, estMSqrd=estMSqrd))
}

estim_M <- function(A, type, MoreArgs){
  n <- MoreArgs$n
  rhoHat_inv <- MoreArgs$rhoHat_inv
  
  #A is a n by n matrix for undirected networks
  switch(type,
         eigs = {
           eig_res <- RSpectra::eigs(A, 1)
           u <- eig_res$vectors
           mult <- eig_res$values
         },
         degs = {
           u <- rowSums(A)
           mult <- (t(u)%*%A%*%u)/(sum(u*u))^2
         }
  )
  
  #Calculation bandwidth
  u <- sort(u)
  uMid <- u[MoreArgs$midPt]
  lmfit.coef <- .lm.fit(cbind(1,1:length(uMid)), uMid)$coefficient
  estMSqrd <- 2*mult^2*(lmfit.coef[2]*length(uMid)/2+lmfit.coef[1])^2*lmfit.coef[2]^2*rhoHat_inv^2*(n+1)^2
  
  return(estMSqrd)
}

initialize_index <- function(A, n, h, verbose){
  ##########################################################################
  # exponential Taylor approximation to L_ij = exp(-||A_i. - A_j.||^2 / 2) = 1 -||A_i. - A_j.||^2 for small ||.||
  # Temporary (use the densest layer for h)
  ##########################################################################
  
  # exponential Taylor approximation to L_ij = exp(-||A_i. - A_j.||^2 / 2) = 1 -||A_i. - A_j.||^2 for small ||.||
  L <- 1 - (.hamming_dist_adj_mat(A)/n)^2 
  d <- rowSums(L)
  L <- outer(d^(-1/2), d^(-1/2))*L - sqrt(d)%o%sqrt(d)/sqrt(sum(d^2))
  eigen_res <- RSpectra::eigs_sym(L, 1) 
  u <- eigen_res$vectors[,1] * sign(eigen_res$vectors[1,1])
  u <- order(u) #Index vectors from smallest to largest.
  k <- floor(n/h)
  remainder <- n - h*k
  
  idx <- rep(0,n)
  for(i in 1:k){
    idx[u[((i-1)*h+1):(i*h)]] <- i
  }
  if(remainder > 0){
    idx[u[(h*k + 1):n]] <- k
  }
  
  return(idx)
}

.prob_mat_from_adj<-function(A, idx){
  K <- max(idx)  
  p_mat <- matrix(0,K,K)
  
  numer <- rowsum(t(rowsum(A, idx)), idx)
  
  bin_size <- table(idx)
  denom <- outer(bin_size,bin_size)
  diag(denom) <- diag(denom) - bin_size
  p_mat<- numer/denom
  dimnames(p_mat)<- list(1:K, 1:K)
  return(p_mat)
}


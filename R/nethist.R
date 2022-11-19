##' Network histogram estimation
##'
##' Estimating network histogram and returning the indices of partitions.
##'
##' @param A An adjacency matrix or igraph object. It must be an undirected and simple graph.
##' @param h A bandwidth parameter. If `NA`, selecting bandwidth by Olhede and Wolfe (2014). If specified, the user input value is used.
##' @param outfile A filename for saving cluster indices. If it is missing, the results are not saved in a file.
##' @param verbose logical value indicating whether verbose output is generated.
##' @returns 
##' An object of class ``nethist``:
##' 
##' \itemize{
##' \item `cluster` a vector of partition indices.
##' \item `p_mat` a probability matrix from network histogram ordered by cluster labels. 
##' \item `rho_hat` estimated sparsity parameter. 
##' }
##' @details 
##' Note that `cluster` only shows a partition of vertices, and the index labels is not an ordered variable. For example, nodes in cluster 1 do not have to more similar to nodes in cluster 2 than nodes in cluster 10. Hence, users would use a user-specified order in [plot.nethist()].
##' 
##' ``nethist()`` runs the following two steps:
##' \itemize{
##' \item Initialization: With spectral clustering, find an initial partitions of vertices. If bandwidth `h` is not specified, automatic bandwidth selection is used.
##' \item Using a greedy search algorithm, find a (local) optima of normalized profile log-likelihood, which is `cluster` in the ``nethist`` object.
##' }
##' @seealso [plot.nethist()]
##' @references Olhede, S. C., & Wolfe, P. J. (2014). Network histograms and universality of blockmodel approximation. Proceedings of the National Academy of Sciences, 111(41), 14722-14727.
##' @references Wolfe, P. J., & Olhede, S. C. (2013). Nonparametric graphon estimation. arXiv preprint arXiv:1309.5936.
##' @examples
##' \dontrun{
##' set.seed(2022)
##' #Generating Erdos-Renyi graph
##' A <- igraph::sample_gnp(100, 0.05)
##' 
##' #With automatic bandwidth selection
##' hist_A <- nethist(A) 
##' 
##' #with user-specified bandwidth
##' hist_A <- nethist(A, h = 20) 
##' 
##' #with adjancency matrix
##' hist_A <- nethist(igraph::as_adj(A)) 
##' }
##' @importFrom stats .lm.fit dist pnorm
##' @importFrom graphics par
##' @importFrom utils write.table
##' @importFrom RSpectra eigs
##' @export
nethist<-function(A, h = NA, outfile, verbose = FALSE){
  UseMethod("nethist")
}
##' @exportS3Method 
nethist.igraph<-function(A, h, outfile, verbose){
  args <- as.list(environment())
  args$A <- igraph::as_adj(args$A, sparse = FALSE)
  do.call("nethist.default", args)
}
##' @exportS3Method 
nethist.matrix<-function(A, h, outfile, verbose){
  args <- as.list(environment())
  do.call("nethist.default", args)
}
##' @exportS3Method 
nethist.dgCMatrix<-function(A, h, outfile, verbose){
  args <- as.list(environment())
  args$A <- as.matrix(args$A)
  do.call("nethist.default", args)
}
##' 
nethist.default <- function(A, h = NA, outfile, verbose = F){
  if(!.is_undirected_simple(A)) stop("Network A must be an undirected simple network.")
  
  # Compute necessary summaries from A
  n <- dim(A)[1]
  rhoHat <- sum(A)/(n*(n-1))
  
  ##########################################################################
  # Pick an analysis bandwidth and initialize via regularized spectral clustering
  ##########################################################################
  if(is.na(h)){
    c <- min(4, sqrt(n)/8)
    h <- .oracbwplugin(A, c, 'degs', 1)$h
    if(verbose) message(paste("Determining bandwidth from data:", round(h)))
  }else{
    if(verbose) message(paste("Determining bandwidth from user input:", round(h)))
  }
  
  h <- max(2, min(n, round(h)))
  if(verbose) message(paste("Final bandwidth:", h))
  
  lastGroupSize <- n %% h
  # step down h, to avoid singleton final group
  while((lastGroupSize==1) & (h>2)){
    h <- h-1
    lastGroupSize <- n %% h
    if(verbose) message('NB: Bandwidth reduced to avoid singleton group')
  }
  if(verbose) message(paste0('Adjacency matrix has ', n, ' rows/cols'))
  
  # Initialize using regularized spectral clustering based on row similarity
  tstart <- Sys.time()
  
  # exponential Taylor approximation to L_ij = exp(-||A_i. - A_j.||^2 / 2) = 1 -||A_i. - A_j.||^2 for small ||.||
  L <- 1 - (.hamming_dist_adj_mat(A)/n)^2 
  d <- rowSums(L)
  L <- outer(d^(-1/2), d^(-1/2))*L - sqrt(d)%o%sqrt(d)/sqrt(sum(d^2))
  eigen_res <- RSpectra::eigs_sym(L, 1) # 2nd eigenvector of normalized Laplacian
  rm(L)
  u <- eigen_res$vectors[,1] * sign(eigen_res$vectors[1,1])
  ind <- order(u) #Index vectors from smallest to largest.
  k <- ceiling(n/h)
  
  idxInit = rep(0,n)
  for(i in 1:k){
    idxInit[ind[((i-1)*h+1):min(n,i*h)]] = i
  }
  if(verbose) message(paste0('Initial label vector assigned from row-similarity ordering; time ',
                 round(difftime(Sys.time(),tstart),4), ' sec'))
  
  idx <- .graphest_fastgreedy(A,h,idxInit, verbose)
  
  if(!missing(outfile)){
    write.table(file=outfile, x = idx, row.names = FALSE, col.names = FALSE)
  }
  
  p_mat <- .prob_mat_from_adj(A,idx)
  
  result <- list(cluster = as.vector(idx), 
                 p_mat = p_mat,
                 rho_hat = rhoHat)
  result <- structure(result, class="nethist")
  return(result)
}

.oracbwplugin <- function(A,c,type, alpha){
  if(missing(type)) type <- 'degs'
  if(missing(alpha)) alpha <- 1
  
  n <- dim(A)[1]
  midPt <- seq(round(n/2-c*sqrt(n),0), round(n/2+c*sqrt(n),0))
  selfLoops <- any(diag(A)!=0)
  sampleSize <- choose(n + selfLoops, 2)
  rhoHat <- sum(A[upper.tri(A, diag = selfLoops)])/sampleSize
  
  if(rhoHat == 0){
    rhoHat_inv <- 0
  }else{
    rhoHat_inv <- 1/rhoHat
  }
  
  #Rank-1 graphon estimate via fhat(x,y) = mult*u(x)*u(y)*pinv(rhoHat);
  if(type=="eigs"){
    eig_res <- RSpectra::eigs(A, 1)
    u <- eig_res$vectors
    mult <- eig_res$values
  }else if(type=='degs'){
    u <- rowSums(A)
    mult <- (t(u)%*%A%*%u)/(sum(u*u))^2
  }else{
    stop(paste("Invalid input type",type))
  }
  
  #Calculation bandwidth
  u <- sort(u)
  uMid <- u[midPt]
  lmfit.coef <- .lm.fit(cbind(1,1:length(uMid)), uMid)$coefficient
  
  if(alpha != 1) stop("Currently only supports alpha = 1")
  h <- (2^(alpha+1)*alpha*mult^2*(lmfit.coef[2]*length(uMid)/2+lmfit.coef[1])^2*lmfit.coef[2]^2*rhoHat_inv)^(-1/(2*(alpha+1)))
  
  estMSqrd <- 2*mult^2*(lmfit.coef[2]*length(uMid)/2+lmfit.coef[1])^2*lmfit.coef[2]^2*rhoHat_inv^2*(n+1)^2
  MISEfhatBnd <- estMSqrd*((2/sqrt(estMSqrd))*(sampleSize*rhoHat)^(-1/2) + 1/n)
  message(paste("M^2_hat =", round(estMSqrd,3), ", MISE bound_hat=", round(MISEfhatBnd,3)))
  
  #Diagnostic plot (if the code is runned on interactive)
  if(interactive()){
    par(mfrow=c(1,2))
    plot(u, main = "Graphon projection for 
bandwidth estimation", type = 'l')
    plot(uMid, main = "Chosen patch of projection component 
    (adjust using c)",type = 'l')
    par(mfrow=c(1,1))#Reset
  }
  
  return(list(h=h, estMSqrd=estMSqrd))
}
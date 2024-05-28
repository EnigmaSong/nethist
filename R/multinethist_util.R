#' @noRd

# Generalized inverse of rho_hat
ginv <- function(x){
  non_zero_ind <- (x!=0)
  x[non_zero_ind] <- 1/x[non_zero_ind]
  return(x)
}

# checking index order vector in plot.nethist is valid or not
.is_valid_order <- function(ind_order, ind_nethist){
  return(setequal(ind_order, ind_nethist) & (length(ind_order)==length(ind_nethist)))
}

#Convert group labels (1,...,K) to (0, ..., K-1) as C++ index
Rind_to_Cind <- function(x){
  if(any(x < 1)) stop("Invalid index input. Must be integers between 1 to K.")
  return(x-1);
}

#estimate M
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

#Initialize group label index given h
initialize_index <- function(A, n, h, verbose){
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

#Compute block probabilities
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

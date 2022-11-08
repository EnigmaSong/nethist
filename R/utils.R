#' @noRd

# Generalized inverse of rho_hat
.ginv <- function(x){
  non_zero_ind <- (x!=0)
  x[non_zero_ind] <- 1/x[non_zero_ind]
  return(x)
}

# Checking simple & undirected graph
.is_undirected_simple<-function(A){
  if(dim(A)[1]!=dim(A)[2]) {message("A is not a square matrix");return(FALSE)}
  if(!isSymmetric(A)) {message("A is not symmetric.");return(FALSE)}
  if(any((A!=0)&(A!=1))) {message("A is not a binary matrix");return(FALSE)}
  
  return(TRUE)
}

# checking index order vector in plot.nethist is valid or not
.is_valid_order <- function(ind_order, ind_nethist){
  return(setequal(ind_order, ind_nethist) & (length(ind_order)==length(ind_nethist)))
}

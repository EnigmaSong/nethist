#' @noRd

# Generalized inverse of rho_hat
.ginv <- function(x){
  non_zero_ind <- (x!=0)
  x[non_zero_ind] <- 1/x[non_zero_ind]
  return(x)
}

# checking index order vector in plot.nethist is valid or not
.is_valid_order <- function(ind_order, ind_nethist){
  return(setequal(ind_order, ind_nethist) & (length(ind_order)==length(ind_nethist)))
}

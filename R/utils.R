#' @noRd

# checking index order vector in plot.nethist is valid or not
.is_valid_order <- function(ind_order, ind_nethist){
  return(setequal(ind_order, ind_nethist) & (length(ind_order)==length(ind_nethist)))
}

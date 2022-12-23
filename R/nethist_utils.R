##' @noRd

check_input_error<-function(A, h, outfile, verbose){
  if(!.is_undirected_simple(A)) {
    stop("Network A must be an undirected simple network.")
  }
  if(is.integer(h)) {
    warning(paste0("User input h=", h, " is not an integer. Use the nearest integer."))
  }
  if((h <= 1)| (h > dim(A)[1])){
    stop(paste0("User input h=", h, 
                " is an invalid value. Use integers between 2 and ",
                dim(A)[1]))
  }
}
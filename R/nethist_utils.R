##' @noRd

check_input_error<-function(A, h, outfile, verbose){
  if(!.is_undirected_simple(A)) {
    stop("Network A must be an undirected simple network.")
  }
  if(is.integer(h)) {
    warning(paste0("User input h=", h, " is not an integer. Use the nearest integer."))
  }
  if((h <= 1)| (h > n)){
    stop(paste0("User input h=", h, 
                " is an invalid value. Use integers between 2 and ",
                dim(A)[1L]))
  }
}

get_bandwidth<- function(A, h, verbose){
  n <- dim(A)[1L]
  if(is.na(h)){
    h <- .oracbwplugin(A, min(4, sqrt(n)/8), 'degs', 1, rhoHat, verbose)$h
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
  return(h)
}
##' @noRd

check_input_error<-function(A, h, outfile, verbose){
  if(!.is_undirected_simple(A)) {
    stop("Network A must be an undirected simple network.")
  }
  if(!is.na(h)){
    if(is.integer(h)) {
      warning(paste0("User input h=", h, " is not an integer. Use the nearest integer."))
    }
    if((h <= 1)| (h > dim(A)[1L])){
      stop(paste0("User input h=", h, 
                  " is an invalid value. Use integers between 2 and ",
                  dim(A)[1L]))
    }
  }
}

get_bandwidth<- function(A, n, rhoHat, h, verbose){
  if(is.na(h)){
    h <- oracbwplugin(A, min(4, sqrt(n)/8), 'degs', 1, rhoHat, verbose)$h
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

oracbwplugin <- function(A,c,type, alpha,
                          rhoHat, verbose){
  #Assume A is symmetric, simple, and no self-loop
  if(missing(type)) type <- 'degs'
  if(missing(alpha)) alpha <- 1
  
  n <- dim(A)[1]
  midPt <- seq(round(n/2-c*sqrt(n),0), round(n/2+c*sqrt(n),0))
  rhoHat_inv <- .ginv(rhoHat)
  sampleSize <- n*(n-1)/2
  
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
  if(verbose){
    message(paste("M^2_hat =", round(estMSqrd,3), ", MISE bound_hat=", round(MISEfhatBnd,3)))
  }
  
  #Diagnostic plot (if the code is runned on interactive)
  if(interactive()){
    par(mfrow=c(1,2))
    plot(u, main = "Graphon projection for bandwidth estimation", type = 'l')
    plot(uMid, main = "Chosen patch of projection component (adjust using c)",type = 'l')
    par(mfrow=c(1,1))#Reset
  }
  
  return(list(h=h, estMSqrd=estMSqrd))
}

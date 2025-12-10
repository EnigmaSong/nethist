##' @noMd

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
  estMSqrd <- rep(0,L) #estmate of M^2
  
  #Rank-1 graphon estimate via fhat(x,y) = mult*u(x)*u(y)*pinv(rhoHat);
  for(l in 1:L){
    estMSqrd[l] <- estim_MSqrd(A[,,l], type, 
                           MoreArgs = list(n = n, midPt = midPt, 
                                           rhoHat_inv = rhoHat_inv[l]))
  }
  
  mean_estMSqrd <- mean(estMSqrd)
  mean_edgenum <- mean(sampleSize*rhoHat)  
  weigths <- rhoHat/sum(rhoHat)
  wmean_estMSqrd <- weighted.mean(estMSqrd, weigths)
  mean_inv_edgenum <- mean(1/(sampleSize*rhoHat))
  h <- ifelse(common_f,
              sqrt(n)*(2*mean_estMSqrd*sum(rhoHat))^(-0.25),
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
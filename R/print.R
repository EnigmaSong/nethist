##' Print Network histogram type objects
##'
##' Printing ``nethist", ``multinethist", and ``hnethist" 
##'
##' @param x one of ``nethist", ``multinethist", and ``hnethist" objects.
##' @param ... see print()
##' @details 
##' ... 
##' @returns 
##' Prints the objects
##' \itemize{
##' \item `thetahat' a probability array from multinetwork histogram ordered by cluster labels. 
##' }
##' @examples
##' \dontrun{
##'    set.seed(42)
##'    nethist_polblog <- nethist(polblog)
##'    nethist_polblog #print(nethist_polblog)
##'
##'    set.seed(42)
##'    mnethist_Ind_Vil <- multinethist(IndianVil)
##'    mnethist_Ind_Vil #print(mnethist_Ind_Vil)
##' }
##' @importFrom stats heatmap
##' @importFrom graphics text

##' @exportS3Method 
print.nethist <- function(x, ...){
  
  cat("\nthetahat:\n")
  print(x$thetahat, ...)
  
  cat(paste0("\nMethod: ", ifelse(x$method == "PLL", 
                                   "Profile Likelihood", 
                                   "Least Squared Error"),
             "\n"))
  
  cat(paste0("\n", ifelse(x$method == "PLL", 
                          paste("normalized likelihood:", x$normalized_LL, sep="\n"), 
                          paste("Mean squared error:", x$LSE, sep = "\n")),"\n"))
  
  cat("\nAvailable components:\n", sep = "\n")
  print(names(x))
  
  return(invisible(x))
}

##' @exportS3Method 
print.multinethist <- function(x, ...){
  cat("\nTheta_hat:\n")
  print(x$thetahat[,,], ...)
  
  cat(paste0("\nMethod: ", ifelse(x$method == "PLL", 
                                  "Profile Likelihood", 
                                  "Least Squared Error"),
             "\n"))
  
  cat(paste0("\n", ifelse(x$method == "PLL", 
                          paste("normalized likelihood:", x$normalized_LL, sep="\n"), 
                          paste("Mean squared error:", x$LSE, sep = "\n")),"\n"))
  
  cat("\nAvailable components:\n", sep = "\n")
  print(names(x))
  
  return(invisible(x))
}


##' @exportS3Method 
print.hnethist <- function(x, ...){
  cat("\nTheta_hat:\n")
  print(x$thetahat, ...)
  
  cat(paste0("\nMethod: ", ifelse(x$method == "PLL", 
                                  "Profile Likelihood", 
                                  "Least Squared Error"),
             "\n"))
  
  cat(paste0("\n", ifelse(x$method == "PLL", 
                          paste("normalized likelihood:", x$normalized_LL, sep="\n"), 
                          paste("Mean squared error:", x$LSE, sep = "\n")),"\n"))
  
  cat("\nAvailable components:\n", sep = "\n")
  print(names(x))
  
  return(invisible(x))
}

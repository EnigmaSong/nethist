#' Indian Village socioeconomic networks
#' 
#' A processed version of the Indian village dataset from Banerjee et al. (2013),
#' representing socio-economic networks at the individual level for village ID 40.
#' 
#' @name IndianVil
#' @docType data
#' @format \code{IndianVil} is an array with size of 231 x 231 x 12.
#' @details The array represents the socio-economic relationships of individuals from Village ID 40. 
#' Nodes with zero degrees across all layers were removed, reducing the total number of nodes from 241 to 231.
#'
#' @references Banerjee, A., Chandrasekhar, A. G., Duflo, E., & Jackson, M. O. (2013). The diffusion of microfinance. Science, 341(6144), 1236498.  
#' @references Song, Y. & Olhede, S.C. (2026+). Graph Limits for Sparse Multilayer Networks.
#' @source 
#'  https://doi.org/10.7910/DVN/U3BIHX
#' @keywords datasets
#' @examples
#' 
#'   data(IndianVil)
#'   
#'   IndianVil
#' 
NULL
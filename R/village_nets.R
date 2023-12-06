#' Indian village household social network data
#'
#' A subset of data from Indian villages household social networks. 
#' The array is constructed from the household network of Village ID 40 of the dataset. 
#' After removing nodes with zero degrees in all twelve layers, 
#' we obtain two multiplex networks with 231 households.
#'
#' @format ## `village_nets`
#' An array describes a multiplex network consists of 231 rows and columns and 12 layers:
#' \describe{
#'  Rows and columns of the array are households. Slices of the array are twelve types of layers.
#'  
#' }
#' @source Banerjee, Abhijit; Chandrasekhar, Arun G.; Duflo, Esther; Jackson, Matthew O., 2013, "The Diffusion of Microfinance", https://doi.org/10.7910/DVN/U3BIHX, Harvard Dataverse, V9
"village_nets"
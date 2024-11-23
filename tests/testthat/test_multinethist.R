set.seed(42)
data(kite, package="igraphdata")
kite <- array(igraph::as_adjacency_matrix(igraph::upgrade_graph(kite),sparse=FALSE), dim=c(10,10,1))
## Later remove this if we add this function in another public package
rnets_graphon<- function(n, num_vertice, graphon_fun, identical_latent_vars=TRUE){
  if(identical_latent_vars){
    latent_samples <- runif(num_vertice)
  }else{
    stop("identical_latent_vars= FALSE is not implemented yet.")
  }
  
  sampled.network <- rep(list(matrix(0, nrow = num_vertice, ncol = num_vertice)),n)
  prob_mat <- outer(latent_samples, latent_samples, graphon_fun)
  if(any((prob_mat>1)|(prob_mat<0))) {
    stop("Error: There is an entry that prob_mat > 1 or <0")
  }
  
  for(i in 1:n){
    for(row_ind in 1:(num_vertice-1)){
      for(col_ind in (row_ind+1):num_vertice){
        sampled.network[[i]][row_ind,col_ind] <- rbinom(1, size = 1, prob = prob_mat[row_ind,col_ind])
        sampled.network[[i]][col_ind,row_ind] <- sampled.network[[i]][row_ind,col_ind]
      }
    }
  }
  
  return(sampled.network)
}
####
sample_mnet<-rnets_graphon(5, 200, function(x,y) pmin(x,y))
array_mnet <- array(0, dim = c(200,200,5))
for(l in 1:length(sample_mnet)){
  array_mnet[,,l] = sample_mnet[[l]]
}

array_mnet_diffrho <- array(0, dim = c(200,200,5))
for(l in 1:length(sample_mnet)){
  array_mnet_diffrho[,,l] = rnets_graphon(5, 200, function(x,y) (0.25+0.15*l)*pmin(x,y))[[1]]
}

## Main methods ----
### Checking invalid input ----
# test_that("LSE is not used for multinethist (Temperally)",{
#     expect_error(multinethist(array_mnet, method = "LSE"))
# })
test_that("non-simple graph input for A (self-loop)", {
  A <- kite
  A[1,1,1] <- 1
  expect_error(multinethist(A))
})
test_that("non-simple graph input for A (asymmetric)", {
  A <- kite
  A[1,2,1] <- 1; A[2,1,1] <- 0
  expect_error(multinethist(A))
})
test_that("A is not supported object", {
  A <- structure(list(kite[,,1]), class = "test_object")
  expect_error(multinethist(A))
})

### Checking default value in argument ----
test_that("multinethist (one layer)", {
  expect_no_error(multinethist(kite,h=5, max_itr = 1000))
  expect_no_error(multinethist(kite,h=5, method = "LSE", max_itr = 1000))
})

test_that("multinethist (5-layers)", {
  expect_no_error(multinethist(array_mnet, max_itr = 1000))
  expect_no_error(multinethist(array_mnet,common_f=TRUE, max_itr = 1000))
})

test_that("multinethist (5-layers w/ different sparsity)", {
  expect_no_error(multinethist(array_mnet_diffrho, max_itr = 1000))
  expect_no_error(multinethist(array_mnet_diffrho,common_f=TRUE, max_itr = 1000))
})

test_that("nethist (checking default argument)", {
  expect_equal({set.seed(42); multinethist(kite,h=5, max_itr = 1000)},
               {set.seed(42); multinethist(kite,h=5, method = "PLL", max_itr = 1000)})
})


# test_that("No infinite loop", {
#   expect_no_error(multinethist(kite,h=10, max_itr = 1000))
# })

# Unit tests

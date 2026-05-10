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
test_that("multinethist (one layer, both PLL/LSE)", {
  expect_no_error(multinethist(kite,h=5, max_itr = 1000))
  expect_no_error(multinethist(kite,h=5, method = "LSE", max_itr = 1000))
})

test_that("multinethist (general/homogeneous, 5-layers)", {
  expect_no_error(multinethist(array_mnet, max_itr = 1000))
  expect_no_error(multinethist(array_mnet,common_f=TRUE, max_itr = 1000))
})

test_that("multinethist (general/homogeneous, 5-layers w/ different sparsity)", {
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

## Consistency: nethist(matrix) == multinethist(matrix) ----
test_that("nethist and multinethist give identical results for 2D matrix input", {
  A_kite <- kite[, , 1]
  result_nethist <- {set.seed(42); nethist(A_kite, h = 5, max_itr = 1000)}
  result_multi   <- {set.seed(42); multinethist(A_kite, h = 5, max_itr = 1000)}

  expect_equal(result_nethist$cluster,      result_multi$cluster)
  expect_equal(result_nethist$thetahat,     result_multi$thetahat)
  expect_equal(result_nethist$rho_hat,      result_multi$rho_hat)
  expect_equal(result_nethist$normalized_LL, result_multi$normalized_LL)
})

## network / list / combined_networks input ----
test_that("network input gives identical result to matrix (single layer)", {
  skip_if_not_installed("network")
  A <- kite[, , 1L]
  A_net <- network::as.network(A, directed = FALSE)
  expect_equal({set.seed(1L); multinethist(A_net, h = 5L, max_itr = 1000L)},
               {set.seed(1L); multinethist(A,     h = 5L, max_itr = 1000L)})
})

test_that("list of matrices gives identical result to array", {
  A1 <- array_mnet[, , 1L]
  A2 <- array_mnet[, , 2L]
  ref <- array(c(A1, A2), dim = c(nrow(A1), ncol(A1), 2L))
  expect_equal({set.seed(1L); multinethist(list(A1, A2), h = 10L, max_itr = 1000L)},
               {set.seed(1L); multinethist(ref,          h = 10L, max_itr = 1000L)})
})

test_that("list of igraph gives identical result to array", {
  skip_if_not_installed("igraph")
  A1 <- array_mnet[, , 1L]
  A2 <- array_mnet[, , 2L]
  g1 <- igraph::graph_from_adjacency_matrix(A1, mode = "undirected", diag = FALSE)
  g2 <- igraph::graph_from_adjacency_matrix(A2, mode = "undirected", diag = FALSE)
  ref <- array(c(A1, A2), dim = c(nrow(A1), ncol(A1), 2L))
  expect_equal({set.seed(1L); multinethist(list(g1, g2), h = 10L, max_itr = 1000L)},
               {set.seed(1L); multinethist(ref,          h = 10L, max_itr = 1000L)})
})

test_that("list of network gives identical result to array", {
  skip_if_not_installed("network")
  A1 <- array_mnet[, , 1L]
  A2 <- array_mnet[, , 2L]
  n1 <- network::as.network(A1, directed = FALSE)
  n2 <- network::as.network(A2, directed = FALSE)
  ref <- array(c(A1, A2), dim = c(nrow(A1), ncol(A1), 2L))
  expect_equal({set.seed(1L); multinethist(list(n1, n2), h = 10L, max_itr = 1000L)},
               {set.seed(1L); multinethist(ref,          h = 10L, max_itr = 1000L)})
})

test_that("combined_networks gives identical result to array", {
  skip_if_not_installed("network")
  skip_if_not_installed("ergm.multi")
  A1 <- array_mnet[, , 1L]
  A2 <- array_mnet[, , 2L]
  n1   <- network::as.network(A1, directed = FALSE)
  n2   <- network::as.network(A2, directed = FALSE)
  nets <- ergm.multi::Networks(list(n1, n2))
  ref  <- array(c(A1, A2), dim = c(nrow(A1), ncol(A1), 2L))
  expect_equal({set.seed(1L); multinethist(nets, h = 10L, max_itr = 1000L)},
               {set.seed(1L); multinethist(ref,  h = 10L, max_itr = 1000L)})
})

test_that("vertex count mismatch in list errors", {
  A1 <- array_mnet[, , 1L]
  A2 <- array_mnet[seq_len(nrow(A1) - 1L), seq_len(nrow(A1) - 1L), 1L]
  expect_error(multinethist(list(A1, A2), h = 10L))
})

test_that("combined_networks element inside list errors", {
  skip_if_not_installed("network")
  skip_if_not_installed("ergm.multi")
  A1   <- array_mnet[, , 1L]
  n1   <- network::as.network(A1, directed = FALSE)
  nets <- ergm.multi::Networks(list(n1, n1))
  expect_error(multinethist(list(nets, nets), h = 10L),
               regexp = "combined_networks")
})

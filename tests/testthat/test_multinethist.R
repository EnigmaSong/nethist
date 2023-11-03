set.seed(42)
data(kite, package="igraphdata")
kite <- array(igraph::as_adj(kite,sparse=FALSE), dim=c(10,10,1))

sample_mnet<-rnets_graphon(5, 200, function(x,y) pmin(x,y))
array_mnet <- array(0, dim = c(200,200,5))
for(l in 1:length(sample_mnet)){
  array_mnet[,,l] = sample_mnet[[l]]
}

array_mnet_diffrho <- array(0, dim = c(200,200,5))
for(l in 1:length(sample_mnet)){
  array_mnet_diffrho[,,l] = rnets_graphon(5, 200, function(x,y) (0.25+0.15*l)*pmin(x,y))[[1]]
}

## Main methods

test_that("multinethist (one layer)", {
  expect_no_error(multinethist(kite,h=5, max_itr = 1000))
})

test_that("multinethist (5-layers)", {
  expect_no_error(multinethist(array_mnet, max_itr = 1000))
  expect_no_error(multinethist(array_mnet,common_f=TRUE, max_itr = 1000))
})

test_that("multinethist (5-layers w/ different sparsity)", {
  expect_no_error(multinethist(array_mnet_diffrho, max_itr = 1000))
  expect_no_error(multinethist(array_mnet_diffrho,common_f=TRUE, max_itr = 1000))
})

test_that("No infinite loop", {
  expect_no_error(multinethist(kite,h=10, max_itr = 1000))
})

# Unit tests

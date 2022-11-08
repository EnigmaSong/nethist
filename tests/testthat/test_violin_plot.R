G<-polblog

#Network summary plot
test_that("violin plots with pre-specified subsample size",
          {expect_no_error({violin_netsummary(igraph::as_adj(G,sparse = FALSE), subsample_sizes = 100,
                  max_cycle_order = 7, save.plot = FALSE)})}
)

test_that("violin plots with auto-selected subsample size",
          {expect_no_error({violin_netsummary(igraph::as_adj(G,sparse = FALSE), Ns=5,
                             max_cycle_order = 7, save.plot = FALSE)})}
)

test_that("violin plots with igraph object",
          {expect_no_error({violin_netsummary(G, Ns=5,
                             max_cycle_order = 7, save.plot = FALSE)})}
)

test_that("violin plots with sparse matrix object",
          {expect_no_error({violin_netsummary(igraph::as_adj(G), Ns=5,
                             max_cycle_order = 7, save.plot = FALSE)})}
)
data(karate,package="igraphdata")
G<-igraph::upgrade_graph(karate)

#Network summary plot
test_that("violin plots with pre-specified subsample size",
          {suppressMessages(expect_no_error({violin_netsummary(igraph::as_adjacency_matrix(G,sparse = FALSE), subsample_sizes = 20,
                  max_cycle_order = 7, save.plot = FALSE)}))}
)

test_that("violin plots with auto-selected subsample size",
          {suppressMessages(expect_no_error({violin_netsummary(igraph::as_adjacency_matrix(G,sparse = FALSE), Ns=5,
                             max_cycle_order = 7, save.plot = FALSE)}))}
)

test_that("violin plots with igraph object",
          {suppressMessages(expect_no_error({violin_netsummary(G, Ns=5,
                             max_cycle_order = 7, save.plot = FALSE)}))}
)

test_that("violin plots with sparse matrix object",
          {suppressMessages(expect_no_error({violin_netsummary(igraph::as_adjacency_matrix(G), Ns=5,
                             max_cycle_order = 7, save.plot = FALSE)}))}
)

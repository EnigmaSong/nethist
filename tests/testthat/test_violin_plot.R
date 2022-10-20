G<-graph_from_edgelist(as.matrix(polblog),directed = FALSE)
G<-delete.vertices(G, degree(G)==0)

#Network summary plot
test_that("violin plots with pre-specified subsample size",
          {violin_netsummary(igraph::as_adj(G,sparse = FALSE), subsample_sizes = 100,
                  max_cycle_order = 7, save.plot = FALSE)}
)

test_that("violin plots with auto-selected subsample size",
          {violin_netsummary(igraph::as_adj(G,sparse = FALSE), Ns=5,
                             max_cycle_order = 7, save.plot = FALSE)}
)

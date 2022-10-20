G<-graph_from_edgelist(as.matrix(polblog),directed = FALSE)
G<-delete.vertices(G, degree(G)==0)
A<-igraph::as_adj(G,sparse=FALSE)
#Network histogram
test_that("network histogram estimation",
          {
            hist_G <- nethist(A, h = 72)
            
          }
)

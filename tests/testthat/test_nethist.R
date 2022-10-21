G<-igraph::graph_from_edgelist(as.matrix(polblog),directed = FALSE)
G<-igraph::delete.vertices(G, igraph::degree(G)==0)
A<-igraph::as_adj(G,sparse=FALSE)
h_used <- 72

hist_G <- nethist(A, h = h_used)
bin_size <- table(hist_G$cluster)
num_bins <- length(bin_size)

test_that("Check the estimated probablity matrix is symmetric",
          {
            expect_equal(isSymmetric(hist_G$p_mat), expected = TRUE)
          }
)

test_that("Check the size of bin is equal to h_used (except the last group)",
          {
            
            expect_equal(all(bin_size[num_bins-1] == h_used), expected = TRUE)
          }
)
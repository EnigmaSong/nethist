# Test code for .count_k_cycle and related hidden functions
test_that("Check cycle counts from a complete graph with n = 10.",
          {
            n <- 10
            G <- igraph::make_full_graph(n, directed = FALSE, loop = FALSE)
            true_num_cycles <- c(.ffct(n,3)/6, .ffct(n,4)/8, .ffct(n,5)/10, .ffct(n,6)/12, .ffct(n,7)/14)
            expect_equal(as.vector(.count_k_cycle(igraph::as_adj(G,sparse= FALSE), 7)), 
                         expected = true_num_cycles)
          }
)


test_that("Check cycle counts from a subset of political blog dataset",
          {
            G <- igraph::as_adj(polblog, sparse = FALSE)[1:50,1:50]
            true_num_cycles <- c(45,159,460,1223,2954)
            expect_equal(as.vector(.count_k_cycle(G, 7)), 
                         expected = true_num_cycles)
          }
)
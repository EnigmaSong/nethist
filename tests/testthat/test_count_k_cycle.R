# Test code for .count_k_cycle and related hidden functions
n <- 10

test_that("Check falling factorial is correct",
          {
            expect_equal(c(.ffct(10,3), .ffct(20, 5)), c(10*9*8, 20*19*18*17*16))
          }
)

test_that("Check cycle counts from a complete graph with n = 10.",
          {
            G <- igraph::make_full_graph(n, directed = FALSE, loop = FALSE)
            true_num_cycles <- c(.ffct(n,3)/6, .ffct(n,4)/8, .ffct(n,5)/10, .ffct(n,6)/12, .ffct(n,7)/14)
            expect_equal(as.vector(.count_k_cycle(igraph::as_adj(G,sparse= FALSE), 7)), 
                         expected = true_num_cycles)
          }
)
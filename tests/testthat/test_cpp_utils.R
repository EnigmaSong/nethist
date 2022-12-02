# Test code for utility functions written in Rcpp
#Binary adjacency mat
A1 <- matrix(c(0,1,1,0,0,
              1,0,0,1,1,
              1,0,0,0,0,
              0,1,0,0,0,
              0,1,0,0,0),5,5)
#adj with self-loop
A2 <- A1
A2[3,3] <- 1
#adj for non-simple graph
A3 <- A1
A3[3,4] <- 2
A3[4,3] <- 2

test_that("Check hamming distance of binary matrix",
          {
            manhattan_dist_A1 <- as.matrix(dist(A1,method = "manhattan"))
            dimnames(manhattan_dist_A1) <- NULL
            expect_equal(.hamming_dist_adj_mat(A1), 
                         manhattan_dist_A1)
          }
)

test_that("Check .is_undirected_simple() for binary matrix with self-loop",
          {
            expect_false(.is_undirected_simple(A2))
          }
)
test_that("Check .is_undirected_simple() for a multigraph",
          {
            expect_false(.is_undirected_simple(A3))
          }
)

test_that("Check falling factorial is correct",
          {
            expect_equal(c(.ffct(10,3), .ffct(20, 5)), c(10*9*8, 20*19*18*17*16))
          }
)
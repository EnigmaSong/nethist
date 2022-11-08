G<-igraph::sample_gnp(100, p = 0.1)
h_used <- 10

hist_G <- nethist(G, h = h_used)
bin_size <- table(hist_G$cluster)
num_bins <- length(bin_size)

test_that("Check any error when the input is igraph",
          {
            expect_no_error(nethist(G, h = h_used))
          }
)

test_that("Check any error when the input is sparse matrix",
          {
            expect_no_error(nethist(igraph::as_adj(G), h = h_used))
          }
)

test_that("Check any error when the input is dense matrix",
          {
            expect_no_error(nethist(igraph::as_adj(G, sparse = FALSE), h = h_used))
          }
)

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

test_that("Check rho_hat is between 0 and 1",
          {
            
            expect_true((hist_G$rho_hat <=1)&(hist_G$rho_hat>=0))
          }
)
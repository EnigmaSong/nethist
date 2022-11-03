G<-polblog
h_used <- 72

hist_G <- nethist(G, h = h_used)
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

test_that("Check rho_hat is between 0 and 1",
          {
            
            expect_true((hist_G$rho_hat <=1)&(hist_G$rho_hat>=0))
          }
)
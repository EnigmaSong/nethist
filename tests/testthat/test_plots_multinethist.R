set.seed(42)
# Delete this part later when this function is included in other public package
rnets_graphon <- function(n, num_vertice, graphon_fun, identical_latent_vars=TRUE){
  if(identical_latent_vars){
    latent_samples <- runif(num_vertice)
  }else{
    stop("identical_latent_vars= FALSE is not implemented yet.")
  }

  sampled.network <- rep(list(matrix(0, nrow = num_vertice, ncol = num_vertice)),n)
  prob_mat <- outer(latent_samples, latent_samples, graphon_fun)
  if(any((prob_mat>1)|(prob_mat<0))) {
    stop("Error: There is an entry that prob_mat > 1 or <0")
  }

  for(i in 1:n){
    for(row_ind in 1:(num_vertice-1)){
      for(col_ind in (row_ind+1):num_vertice){
        sampled.network[[i]][row_ind,col_ind] <- rbinom(1, size = 1, prob = prob_mat[row_ind,col_ind])
        sampled.network[[i]][col_ind,row_ind] <- sampled.network[[i]][row_ind,col_ind]
      }
    }
  }

  return(sampled.network)
}
#
sample_mnet<-rnets_graphon(5, 200, function(x,y) pmin(x,y))
array_mnet <- array(0, dim = c(200,200,5))
for(l in 1:length(sample_mnet)){
  array_mnet[,,l] = sample_mnet[[l]]
}

array_mnet_diffrho <- array(0, dim = c(200,200,5))
for(l in 1:length(sample_mnet)){
  array_mnet_diffrho[,,l] = rnets_graphon(5, 200, function(x,y) (0.25+0.15*l)*pmin(x,y))[[1]]
}

covar         <- factor(c(rep("a", 100), rep("b", 100)))
covar_numeric <- runif(200)

mnhist_id    <- multinethist(array_mnet_diffrho,
                              control = nethist_control(max_itr = 1e3))
mnhist_id_cf <- multinethist(array_mnet_diffrho,
                              control = nethist_control(max_itr = 1e3),
                              common_f = TRUE)
nhist_single <- nethist(sample_mnet[[1]],
                        control = nethist_control(max_itr = 1e3))


#Summary plots
test_that("check no errors in covariate_plots",
          {
            expect_no_error(covariate_plot(mnhist_id, covar))
            expect_no_error(covariate_plot(mnhist_id_cf, covar))
          }
          )


test_that("check no errors in covariate_plots with custom order",
          {
            expect_no_error(covariate_plot(mnhist_id, covar, idx_order = unique(mnhist_id$cluster)))
            expect_no_error(covariate_plot(mnhist_id_cf, covar, idx_order = unique(mnhist_id_cf$cluster)))
          }
)

#heatmap
test_that("check no errors in summary_plots",
          {
            expect_no_error(plot(mnhist_id, idx_order = seq_len(max(mnhist_id$cluster))))
            expect_no_error(plot(mnhist_id_cf, idx_order = seq_len(max(mnhist_id_cf$cluster))))
          }
)


test_that("check no errors in summary_plots with custom order",
          {
            expect_no_error(plot(mnhist_id, idx_order = unique(mnhist_id$cluster)))
            expect_no_error(plot(mnhist_id_cf, idx_order = unique(mnhist_id_cf$cluster)))
          }
)

# covariate_plot: numeric covariate (violin path)
test_that("covariate_plot works with numeric covariate",
          {
            expect_no_error(covariate_plot(mnhist_id, covar_numeric))
          }
)

# covariate_plot: covariate length mismatch
test_that("covariate_plot errors on covariate length mismatch",
          {
            expect_error(covariate_plot(mnhist_id, covar_numeric[1:50]))
          }
)

# covariate_plot: invalid idx_order warns and falls back
test_that("covariate_plot warns on invalid idx_order",
          {
            expect_warning(covariate_plot(mnhist_id, covar, idx_order = c(99L, 98L)))
          }
)

# covariate_plot: nethist class dispatches correctly
test_that("covariate_plot works with nethist object",
          {
            expect_no_error(covariate_plot(nhist_single, covar))
          }
)

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

covar <- factor(c(rep("a",100), rep("b",100)))

mnhist_id <- multinethist(array_mnet_diffrho, max_itr = 1e3)
mnhist_id_cf <- multinethist(array_mnet_diffrho, max_itr = 1e3, common_f = TRUE)


#Summary plots
test_that("check no errors in summary_plots",
          {
            expect_no_error(summary_plot(mnhist_id, covar))
            expect_no_error(summary_plot(mnhist_id_cf, covar))
          }
          )


test_that("check no errors in summary_plots with custom order",
          {
            expect_no_error(summary_plot(mnhist_id, covar, idx_order = unique(mnhist_id$cluster)))
            expect_no_error(summary_plot(mnhist_id_cf, covar, idx_order = unique(mnhist_id_cf$cluster)))
          }
)

#heatmap
test_that("check no errors in summary_plots",
          {
            expect_no_error(plot(mnhist_id, idx_order = 1:16))
            expect_no_error(plot(mnhist_id_cf, idx_order = 1:33))
          }
)


test_that("check no errors in summary_plots with custom order",
          {
            expect_no_error(plot(mnhist_id, idx_order = unique(mnhist_id$cluster)))
            expect_no_error(plot(mnhist_id_cf, idx_order = unique(mnhist_id_cf$cluster)))
          }
)

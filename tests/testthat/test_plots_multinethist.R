set.seed(42)

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
            expect_no_error(plot(mnhist_id))
            expect_no_error(plot(mnhist_id_cf))
          }
)


test_that("check no errors in summary_plots with custom order",
          {
            expect_no_error(plot(mnhist_id, idx_order = unique(mnhist_id$cluster)))
            expect_no_error(plot(mnhist_id_cf, idx_order = unique(mnhist_id_cf$cluster)))
          }
)

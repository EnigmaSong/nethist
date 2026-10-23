# Tests for covariate_plot().
data(karate, package="igraphdata")
set.seed(42)
nethist_karate <- multinethist(igraph::as_adjacency_matrix(igraph::upgrade_graph(karate), sparse = FALSE))

#Define covariates
#factor

#Add factor covariate
faction <- factor(igraph::vertex_attr(karate)$Faction)
numeric_covariate <- as.numeric(faction) + rnorm(igraph::vcount(karate), sd = 3)

test_that("factor covariate",
          {
            expect_no_error(covariate_plot(nethist_karate, faction))
          }
)

test_that("numeric covariate",
          {
            expect_no_error(covariate_plot(nethist_karate, numeric_covariate))
          }
)

test_that("Add title",
          {
            expect_no_error(covariate_plot(nethist_karate, faction, main = "Karate"))
          }
)

test_that("using user-defined cluster order",
          {
            expect_no_error(covariate_plot(nethist_karate, faction, idx_order = c(1,3,2)))
          }
)

test_that("Add y-axis labels",
          {
            expect_no_error(covariate_plot(nethist_karate, faction, ylab = "number karate"))
          }
)

test_that("add legend title",
          {
            expect_no_error(covariate_plot(nethist_karate, faction, legend_title = "faction"))
          }
)

test_that("summary_plot is deprecated in favour of covariate_plot",
          {
            expect_warning(summary_plot(nethist_karate, faction),
                           "deprecated")
          }
)

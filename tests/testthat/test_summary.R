data(karate, package="igraphdata")
set.seed(42)
nethist_karate <- nethist(karate)

#Define covariates
#factor

#Add factor covariate
faction <- factor(igraph::vertex_attr(karate)$Faction)
numeric_covariate <- as.numeric(faction) + rnorm(igraph::vcount(karate), sd = 3)

test_that("factor covariate",
          {
            expect_no_error(summary(nethist_karate, faction))
          }
)

test_that("numeric covariate",
          {
            expect_no_error(summary(nethist_karate, numeric_covariate))
          }
)

test_that("Add title",
          {
            expect_no_error(summary(nethist_karate, faction, main = "Karate"))
          }
)

test_that("using user-defined cluster order",
          {
            expect_no_error(summary(nethist_karate, faction, idx_order = c(1,4,2,3)))
          }
)

test_that("Add y-axis labels",
          {
            expect_no_error(summary(nethist_karate, faction, ylab = "number karate"))
          }
)

test_that("add legend title",
          {
            expect_no_error(summary(nethist_karate, faction, legend_title = "faction"))
          }
)

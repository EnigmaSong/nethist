data(karate, package = "igraphdata")
G <- igraph::upgrade_graph(karate)

#Network summary plot
test_that("netsummary_plot with pre-specified subsample size",
          {suppressMessages(expect_no_error({netsummary_plot(
            igraph::as_adjacency_matrix(G, sparse = FALSE),
            subsample_sizes = 20, max_cycle_order = 7, save_plot = FALSE)}))}
)

test_that("netsummary_plot with auto-selected subsample size",
          {suppressMessages(expect_no_error({netsummary_plot(
            igraph::as_adjacency_matrix(G, sparse = FALSE),
            n_subsample_sizes = 5, max_cycle_order = 7, save_plot = FALSE)}))}
)

test_that("netsummary_plot with igraph object",
          {suppressMessages(expect_no_error({netsummary_plot(
            G, n_subsample_sizes = 5, max_cycle_order = 7, save_plot = FALSE)}))}
)

test_that("netsummary_plot with sparse matrix object",
          {suppressMessages(expect_no_error({netsummary_plot(
            igraph::as_adjacency_matrix(G),
            n_subsample_sizes = 5, max_cycle_order = 7, save_plot = FALSE)}))}
)

test_that("netsummary_plot with network object",
          {suppressMessages(expect_no_error({netsummary_plot(
            network::as.network(igraph::as_adjacency_matrix(G, sparse = FALSE)),
            n_subsample_sizes = 5, max_cycle_order = 7, save_plot = FALSE)}))}
)

# Hybrid Network histogram estimation

Estimating hybrid network histogram for single-layer networks and
returning the indices of partitions.

## Usage

``` r
hnethist(A, h = NA, method = "LSE", control = nethist_control(), ...)
```

## Arguments

- A:

  An adjacency matrix or an igraph object. It must be an undirected and
  simple graph.

- h:

  A bandwidth parameter. If `NA`, the bandwidth is selected by Olhede
  and Wolfe (2014). If specified, the user-supplied value is used.

- method:

  Type of loss function for network histogram. Must be one of `LSE`
  (default) or `PLL` for single-layer hybrid network histogram. See
  Details.

- control:

  A control object from
  [`nethist_control`](https://enigmasong.github.io/nethist/reference/nethist_control.md).
  Governs `max_itr`, `greedy_swap_rule`, `greedy_stop_threshold`, and
  `verbose`.

- ...:

  Currently unused.

## Value

A list of class `c("hnethist", "nethist", "network_histogram")` (see
[network_histogram-class](https://enigmasong.github.io/nethist/reference/network_histogram-class.md))
with the following fields:

- `cluster` an integer vector of length n with vertex-level block
  assignments from the initial nethist fit.

- `thetahat` a k-by-k probability matrix of the selected hybrid model,
  where k is the number of nethist blocks.

- `rho_hat` estimated sparsity parameter from the initial nethist fit.

- `normalized_LL` normalized log-likelihood of the selected model,
  divided by the number of edges.

- `MSE` mean squared error of the selected model.

- `method` loss function used (`"LSE"` or `"PLL"`).

- `h` bandwidth used for the initial nethist fit.

- `blockcluster` a `kmeans` object describing how the k-by-k blocks were
  merged into `s` shapes.

- `BIC` BIC value of the selected model.

- `s` number of distinct shapes in the selected model.

- `details` list of all candidate models, one entry per number of shapes
  from 1 up to the maximum considered, each containing `s`,
  `blockcluster`, `thetahat`, `MSE`, `normalized_LL`, and `BIC`.

- `initial` the nethist object used as the starting point for block
  clustering.

## Details

Among the outputs, the best model is selected based on the BIC
criterion.The original reference provided theoretical guarantees for
LSE, but we also allow PLL for the initial nethist fit. The block
clustering and model selection steps are performed on the LSE loss
regardless of the initial method, as the theoretical results pertain to
LSE.Both stages are stochastic. The initial
[`nethist()`](https://enigmasong.github.io/nethist/reference/nethist.md)
fit uses a greedy search that swaps randomly chosen vertex pairs, and
the merging step calls
[`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html), which draws
its starting centres at random. Neither stage seeds itself, so both
consume the global random number stream and the selected number of
shapes `s` depends on its state. Call
[`set.seed()`](https://rdrr.io/r/base/Random.html) immediately before
the fit to make results reproducible, and expect different seeds to give
somewhat different values of `s`.

## References

Verdeyme, A. & Olhede, S. C. (2024). Hybrid of Node and Link Communities
for Graphon Estimation. arXiv:2401.05088

## Examples

``` r
# \donttest{
set.seed(42)
data(polblog)
hfit <- hnethist(polblog)
hfit
#> 
#> Hybrid network histogram fit (hnethist)
#> 
#> Nodes: 1224   Blocks: 16   Bandwidth: 73
#> Method: Least Squared Error
#> Selected shapes: 17 (BIC: 98687.27)
#> thetahat: 16 x 16 probability matrix
#> Mean squared error: 0.01704723
#> 
#> Use summary() for details.
plot(hfit)

# }
```

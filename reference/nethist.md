# Network histogram estimation for single-layer networks

Estimating a network histogram for a single-layer network and returning
the indices of partitions.

## Usage

``` r
nethist(A, h = NA, method = "PLL", control = nethist_control(), ...)
```

## Arguments

- A:

  An adjacency matrix or graph object. Accepted formats: a `matrix`,
  sparse `dgCMatrix`, `igraph` object, or `network` object. Must be an
  undirected simple graph.

- h:

  A bandwidth parameter. If `NA`, the bandwidth is selected by Olhede
  and Wolfe (2014). If specified, the user-supplied value is used.

- method:

  Type of loss function. One of `"PLL"` (default, profile
  log-likelihood) or `"LSE"` (least squares). See Details.

- control:

  A control object from
  [`nethist_control`](https://enigmasong.github.io/nethist/reference/nethist_control.md).
  Governs `max_itr`, `greedy_swap_rule`, `greedy_stop_threshold`, and
  `verbose`.

- ...:

  **\[deprecated\]** Pass `max_itr`, `greedy_swap_rule`,
  `greedy_stop_threshold`, or `verbose` via
  `control = nethist_control(...)` instead.

## Value

An object of class `c("nethist", "network_histogram")` (see
[network_histogram-class](https://enigmasong.github.io/nethist/reference/network_histogram-class.md))
with the following fields:

- `cluster` an integer vector of length \\n\\ with block assignments.

- `thetahat` a \\k \times k\\ probability matrix ordered by group
  labels.

- `rho_hat` estimated sparsity parameter.

- `normalized_LL` normalized log-likelihood.

- `MSE` mean squared error.

- `method` loss function used (`"PLL"` or `"LSE"`).

- `h` bandwidth used for estimation.

## Details

`method = "PLL"` is for Olhede and Wolfe (2014). `method = "LSE"` is for
Gao et al. (2015).

Note that `cluster` labels are not ordered: vertices in cluster 1 are
not necessarily more similar to cluster 2 than to cluster 10. Users may
specify a custom display order in
[`plot.nethist`](https://enigmasong.github.io/nethist/reference/plot.nethist.md).

The greedy search is stochastic. At each iteration it draws a pair of
vertices uniformly at random and swaps their block labels if the
objective improves, so it converges to a local optimum that depends on
the state of the random number stream. Call
[`set.seed`](https://rdrr.io/r/base/Random.html) immediately before the
fit to make results reproducible, and expect different seeds to give
somewhat different partitions.

## References

Olhede, S. C. & Wolfe, P. J. (2014). Network Histograms and Universality
of Blockmodel Approximation. PNAS, 111(41), 14722-14727.
doi:10.1073/pnas.1400374111

Gao, C., Lu, Y., & Zhou, H. H. (2015). Rate-Optimal Graphon Estimation.
The Annals of Statistics, 43(6), 2624-2652. doi:10.1214/15-AOS1354

## See also

[`multinethist`](https://enigmasong.github.io/nethist/reference/multinethist.md),
[`plot.nethist`](https://enigmasong.github.io/nethist/reference/plot.nethist.md),
[`nethist_control`](https://enigmasong.github.io/nethist/reference/nethist_control.md)

## Examples

``` r
# \donttest{
set.seed(42)
data(polblog)
fit <- nethist(polblog)
fit
#> 
#> Network histogram fit (nethist)
#> 
#> Nodes: 1224   Blocks: 16   Bandwidth: 73
#> Method: Profile Likelihood
#> thetahat: 16 x 16 probability matrix
#> Normalized log-likelihood: -2.969964
#> 
#> Use summary() for details.
plot(fit)


fit_h <- nethist(polblog, h = 72)
# }
```

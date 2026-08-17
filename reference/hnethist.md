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

A list of class `c("hnethist", "nethist")` with the following fields:

- `cluster` an integer vector of length n with vertex-level block
  assignments from the initial nethist fit.

- `thetahat` a k-by-k probability matrix of the selected hybrid model,
  where k is the number of nethist blocks.

- `rho_hat` estimated sparsity parameter from the initial nethist fit.

- `normalized_LL` normalized log-likelihood from the initial nethist
  fit.

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
LSE.

## References

Verdeyme, A. & Olhede, S. C. (2024). Hybrid of Node and Link Communities
for Graphon Estimation. arXiv:2401.05088

# Changelog

## nethist 1.0.0

update: 2026/08/04

#### New features

- [`multinethist()`](https://enigmasong.github.io/nethist/reference/multinethist.md)
  - New implementation of Song and Olhede (2026+) for multi-layer
    network histogram estimation.
  - `common_f` option assumes a common graphon across layers.
  - Accepts a wider range of inputs: 3D array, list of
    `igraph`/`network` objects, and
    [`ergm.multi::combined_networks`](https://rdrr.io/pkg/ergm.multi/man/combine_networks.html)
    object.
- [`nethist()`](https://enigmasong.github.io/nethist/reference/nethist.md)
  - New `method = "LSE"` option, implementing Gao et al. (2015) for
    single-layer network histograms.
  - Now implemented as a single-layer wrapper around
    [`multinethist()`](https://enigmasong.github.io/nethist/reference/multinethist.md).
- [`hnethist()`](https://enigmasong.github.io/nethist/reference/hnethist.md):
  new hybrid network histogram method for shape-based model selection.
- [`fitted()`](https://rdrr.io/r/stats/fitted.values.html): new S3
  method for `nethist` and `multinethist` objects.
- New dataset: `IndianVil`, a multilayer network.

#### Renamed functions

- [`summary()`](https://rdrr.io/r/base/summary.html) -\>
  [`covariate_plot()`](https://enigmasong.github.io/nethist/reference/covariate_plot.md)
- [`violin_netsummary()`](https://enigmasong.github.io/nethist/reference/netsummary_plot.md)
  -\>
  [`netsummary_plot()`](https://enigmasong.github.io/nethist/reference/netsummary_plot.md)

#### Bug fixes

- [`netsummary_plot()`](https://enigmasong.github.io/nethist/reference/netsummary_plot.md)
  (formerly
  [`violin_netsummary()`](https://enigmasong.github.io/nethist/reference/netsummary_plot.md)):
  return 0 instead of NaN when a subsampled network is empty.

#### Other changes

- [`netsummary_plot()`](https://enigmasong.github.io/nethist/reference/netsummary_plot.md):
  x-axis values are now shown as igraph plots of v-shapes and cycles
  instead of labels.
- Renamed output field `LSE` to `MSE` in `nethist`/`multinethist`
  objects.
- Documentation terminology unified from “node” to “vertex”.
- [`multinethist()`](https://enigmasong.github.io/nethist/reference/multinethist.md):
  reject plain `matrix`/`Matrix` elements in list input, since vertex
  ordering across layers cannot be verified.

## nethist 0.2.4

update: 2023/04/04

-bug fix in plot()

## nethist 0.2.3

update: 2023/03/15

- Documentation updates
- Internal functions (not visible to users)

## nethist 0.2.2

update: 2022/12/08

- [`plot.nethist()`](https://enigmasong.github.io/nethist/reference/plot.nethist.md)
  - You can display the estimated probabilities in `p_mat` of `nethist`
    objects by setting `prob=TRUE`. Use `prob.col` and `prob.cex` to set
    the color and size of text, respectively.
- [`nethist()`](https://enigmasong.github.io/nethist/reference/nethist.md)
  - Check the convergence criterion every 5 iterations from the fifth
    iteration (5,10,15,…) instead of every 5 iterations from the first
    iteration (1,6,11,…).
  - minor speed improvement.

## nethist 0.2.1

update: 2022/11/26

- Initialization step of
  [`nethist()`](https://enigmasong.github.io/nethist/reference/nethist.md)
  is faster than earlier by computing hamming distance from the (simple)
  adjacency matrix instead of Manhattan distance.

## nethist 0.2.0

update: 2022/11/19

- Computation speed of
  [`nethist()`](https://enigmasong.github.io/nethist/reference/nethist.md)
  and
  [`violin_netsummary()`](https://enigmasong.github.io/nethist/reference/netsummary_plot.md)
  are improved.

#### Bug fixes

- [`nethist()`](https://enigmasong.github.io/nethist/reference/nethist.md):
  computation of the diagonal part o in greedy search *Rcpp* routine is
  fixed.
- [`nethist()`](https://enigmasong.github.io/nethist/reference/nethist.md):
  computation of log-likelihood when there is a smaller size of group is
  fixed, in greedy search *Rcpp* routine.

## nethist 0.1.1

update: 2022/11/08

- Documentations update: provide extra about default values of
  arguments, rewrite examples.

#### Bug fixes

- Resolving default value issue of *R* in
  [`violin_netsummary()`](https://enigmasong.github.io/nethist/reference/netsummary_plot.md).

## nethist 0.1.0

update: 2022/11/04

- Initial release of `nethist` R package (development)

# Changelog

## nethist 1.1.0

#### New features

- [`logLik()`](https://rdrr.io/r/stats/logLik.html): new S3 method for
  `nethist`, `multinethist`, and `hnethist` objects. Returns the
  unnormalized Bernoulli log-likelihood as a `logLik` object with `df`
  and `nobs` attributes, so [`AIC()`](https://rdrr.io/r/stats/AIC.html)
  and [`BIC()`](https://rdrr.io/r/stats/AIC.html) work directly on a
  fit. Use `logLik(fit, normalized = TRUE)` for the per-edge normalized
  value stored in `fit$normalized_LL`.
- [`hnethist()`](https://enigmasong.github.io/nethist/reference/hnethist.md):
  added examples to the help page.
- [`summary()`](https://rdrr.io/r/base/summary.html): new S3 method for
  `nethist`, `multinethist`, and `hnethist` objects. Reports the fit
  dimensions, the edge density, the log-likelihood with its degrees of
  freedom, `AIC`, `BIC`, and, for `hnethist`, the shape selection.

#### Bug fixes

- [`hnethist()`](https://enigmasong.github.io/nethist/reference/hnethist.md):
  `normalized_LL` is now divided by the number of edges `sum(A)/2`,
  matching
  [`nethist()`](https://enigmasong.github.io/nethist/reference/nethist.md)
  and
  [`multinethist()`](https://enigmasong.github.io/nethist/reference/multinethist.md).
  Previously it was divided by `sum(A)`, so reported values were half
  the size of the corresponding
  [`nethist()`](https://enigmasong.github.io/nethist/reference/nethist.md)
  values. Model selection and `BIC` are unaffected.
- [`print()`](https://rdrr.io/r/base/print.html) methods for `nethist`,
  `multinethist`, and `hnethist` now respect the `digits` argument when
  printing the normalized likelihood and the mean squared error.
- [`netsummary_plot()`](https://enigmasong.github.io/nethist/reference/netsummary_plot.md):
  the example using the suggested package **network** is now guarded by
  [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html).
- [`logLik()`](https://rdrr.io/r/stats/logLik.html) on a `multinethist`
  fit that collapses to a single block now reports the correct degrees
  of freedom.
- [`multinethist()`](https://enigmasong.github.io/nethist/reference/multinethist.md):
  `method = "LSE"` with multilayer input now raises an error. Previously
  it was accepted silently, the fit was computed by profile likelihood,
  and the returned object was labelled `"LSE"`.

#### Class hierarchy

- All fitted objects now inherit from a common base class,
  `network_histogram`: `nethist` fits have class
  `c("nethist", "network_histogram")`, `multinethist` fits
  `c("multinethist", "network_histogram")`, and `hnethist` fits
  `c("hnethist", "nethist", "network_histogram")`. Previously `hnethist`
  inherited from `nethist` while `multinethist` was unrelated to either.
- [`logLik()`](https://rdrr.io/r/stats/logLik.html) and the
  [`print()`](https://rdrr.io/r/base/print.html) method of the summary
  objects are defined once on `network_histogram`.
  [`summary.multinethist()`](https://enigmasong.github.io/nethist/reference/summary.network_histogram.md)
  and
  [`summary.hnethist()`](https://enigmasong.github.io/nethist/reference/summary.network_histogram.md)
  extend the base method through
  [`NextMethod()`](https://rdrr.io/r/base/UseMethod.html). Methods whose
  behaviour depends on the layer dimension, namely
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  [`plot3d()`](https://enigmasong.github.io/nethist/reference/plot3d.md),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html), and
  [`print()`](https://rdrr.io/r/base/print.html), remain specific to the
  subclasses.
- See
  [`?"network_histogram-class"`](https://enigmasong.github.io/nethist/reference/network_histogram-class.md)
  for the documented hierarchy.
- Code that tests for an exact class vector, such as
  `identical(class(fit), c("hnethist", "nethist"))`, needs updating. Use
  [`inherits()`](https://rdrr.io/r/base/class.html) instead.

#### User-visible changes

- The [`print()`](https://rdrr.io/r/base/print.html) methods for
  `nethist`, `multinethist`, and `hnethist` are now compact. They
  identify the fit and report `thetahat` by its dimensions only, instead
  of printing the whole probability matrix or array, and the
  `Available components` listing has been dropped. The detailed output
  moved to the new [`summary()`](https://rdrr.io/r/base/summary.html)
  methods. Use `object$thetahat` or
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) for the
  estimated probabilities.

#### Documentation

- [`nethist()`](https://enigmasong.github.io/nethist/reference/nethist.md),
  [`multinethist()`](https://enigmasong.github.io/nethist/reference/multinethist.md),
  and
  [`hnethist()`](https://enigmasong.github.io/nethist/reference/hnethist.md)
  now state in their help pages that the greedy search is stochastic and
  converges to a local optimum, and that
  [`set.seed()`](https://rdrr.io/r/base/Random.html) should be called
  immediately before the fit for reproducible results.
  [`hnethist()`](https://enigmasong.github.io/nethist/reference/hnethist.md)
  additionally notes that its merging step draws
  [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html) starting
  centres from the same stream, so the selected number of shapes depends
  on it.
- The `@returns` sections of
  [`nethist()`](https://enigmasong.github.io/nethist/reference/nethist.md),
  [`multinethist()`](https://enigmasong.github.io/nethist/reference/multinethist.md),
  and
  [`hnethist()`](https://enigmasong.github.io/nethist/reference/hnethist.md)
  now state the full class vector including `network_histogram`, and
  [`multinethist()`](https://enigmasong.github.io/nethist/reference/multinethist.md)
  lists the `method` field.
- Removed the leftover prompt boilerplate from the `polblog` and
  `polblog_edgelist` help pages, as requested by the CRAN team.
- The `polblog` and `polblog_edgelist` examples now call
  [`igraph::delete_vertices()`](https://r.igraph.org/reference/delete_vertices.html)
  instead of the deprecated
  [`igraph::delete.vertices()`](https://r.igraph.org/reference/delete.vertices.html).

## nethist 1.0.0

CRAN release: 2026-09-02

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

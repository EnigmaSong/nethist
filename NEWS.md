# nethist 1.1.0

### New features

-   `logLik()`: new S3 method for `nethist`, `multinethist`, and `hnethist`
    objects. Returns the unnormalized Bernoulli log-likelihood as a `logLik`
    object with `df` and `nobs` attributes, so `AIC()` and `BIC()` work
    directly on a fit. Use `logLik(fit, normalized = TRUE)` for the per-edge
    normalized value stored in `fit$normalized_LL`.
-   `hnethist()`: added examples to the help page.
-   `summary()`: new S3 method for `nethist`, `multinethist`, and `hnethist`
    objects. Reports the fit dimensions, the edge density, the log-likelihood
    with its degrees of freedom, `AIC`, `BIC`, and, for `hnethist`, the shape
    selection.

### Bug fixes

-   `hnethist()`: `normalized_LL` is now divided by the number of edges
    `sum(A)/2`, matching `nethist()` and `multinethist()`. Previously it was
    divided by `sum(A)`, so reported values were half the size of the
    corresponding `nethist()` values. Model selection and `BIC` are unaffected.
-   `print()` methods for `nethist`, `multinethist`, and `hnethist` now respect
    the `digits` argument when printing the normalized likelihood and the mean
    squared error.
-   `netsummary_plot()`: the example using the suggested package **network** is
    now guarded by `requireNamespace()`.
-   `logLik()` on a `multinethist` fit that collapses to a single block now
    reports the correct degrees of freedom.
-   `multinethist()`: `method = "LSE"` with multilayer input now raises an
    error. Previously it was accepted silently, the fit was computed by
    profile likelihood, and the returned object was labelled `"LSE"`.

### Class hierarchy

-   All fitted objects now inherit from a common base class,
    `network_histogram`: `nethist` fits have class
    `c("nethist", "network_histogram")`, `multinethist` fits
    `c("multinethist", "network_histogram")`, and `hnethist` fits
    `c("hnethist", "nethist", "network_histogram")`. Previously `hnethist`
    inherited from `nethist` while `multinethist` was unrelated to either.
-   `logLik()` and the `print()` method of the summary objects are defined once
    on `network_histogram`. `summary.multinethist()` and `summary.hnethist()`
    extend the base method through `NextMethod()`. Methods whose behaviour
    depends on the layer dimension, namely `plot()`, `plot3d()`, `fitted()`,
    and `print()`, remain specific to the subclasses.
-   See `?"network_histogram-class"` for the documented hierarchy.
-   Code that tests for an exact class vector, such as
    `identical(class(fit), c("hnethist", "nethist"))`, needs updating. Use
    `inherits()` instead.

### User-visible changes

-   The `print()` methods for `nethist`, `multinethist`, and `hnethist` are now
    compact. They identify the fit and report `thetahat` by its dimensions only,
    instead of printing the whole probability matrix or array, and the
    `Available components` listing has been dropped. The detailed output moved
    to the new `summary()` methods. Use `object$thetahat` or `fitted()` for the
    estimated probabilities.

### Documentation

-   `nethist()`, `multinethist()`, and `hnethist()` now state in their help
    pages that the greedy search is stochastic and converges to a local
    optimum, and that `set.seed()` should be called immediately before the fit
    for reproducible results. `hnethist()` additionally notes that its merging
    step draws `stats::kmeans()` starting centres from the same stream, so the
    selected number of shapes depends on it.
-   The `@returns` sections of `nethist()`, `multinethist()`, and `hnethist()`
    now state the full class vector including `network_histogram`, and
    `multinethist()` lists the `method` field.
-   Removed the leftover prompt boilerplate from the `polblog` and
    `polblog_edgelist` help pages, as requested by the CRAN team.
-   The `polblog` and `polblog_edgelist` examples now call
    `igraph::delete_vertices()` instead of the deprecated
    `igraph::delete.vertices()`.

# nethist 1.0.0

update: 2026/08/04

### New features

-   `multinethist()`
    -   New implementation of Song and Olhede (2026+) for multi-layer
        network histogram estimation.
    -   `common_f` option assumes a common graphon across layers.
    -   Accepts a wider range of inputs: 3D array, list of `igraph`/`network`
        objects, and `ergm.multi::combined_networks` object.
-   `nethist()`
    -   New `method = "LSE"` option, implementing Gao et al. (2015) for
        single-layer network histograms.
    -   Now implemented as a single-layer wrapper around `multinethist()`.
-   `hnethist()`: new hybrid network histogram method for shape-based model
    selection.
-   `fitted()`: new S3 method for `nethist` and `multinethist` objects.
-   New dataset: `IndianVil`, a multilayer network.

### Renamed functions

-   `summary()` -> `covariate_plot()`
-   `violin_netsummary()` -> `netsummary_plot()`

### Bug fixes

-   `netsummary_plot()` (formerly `violin_netsummary()`): return 0 instead of
    NaN when a subsampled network is empty.

### Other changes

-   `netsummary_plot()`: x-axis values are now shown as igraph plots of
    v-shapes and cycles instead of labels.
-   Renamed output field `LSE` to `MSE` in `nethist`/`multinethist` objects.
-   Documentation terminology unified from "node" to "vertex".
-   `multinethist()`: reject plain `matrix`/`Matrix` elements in list input,
    since vertex ordering across layers cannot be verified.

# nethist 0.2.4

update: 2023/04/04

-bug fix in plot()


# nethist 0.2.3

update: 2023/03/15
 
-  Documentation updates
-  Internal functions (not visible to users)

# nethist 0.2.2

update: 2022/12/08

-   `plot.nethist()`
    -   You can display the estimated probabilities in `p_mat` of `nethist` objects by setting `prob=TRUE`. Use `prob.col` and `prob.cex` to set the color and size of text, respectively.
-   `nethist()`
    -   Check the convergence criterion every 5 iterations from the fifth iteration (5,10,15,...) instead of every 5 iterations from the first iteration (1,6,11,...).
    -   minor speed improvement.

# nethist 0.2.1

update: 2022/11/26

-   Initialization step of `nethist()` is faster than earlier by computing hamming distance from the (simple) adjacency matrix instead of Manhattan distance.

# nethist 0.2.0

update: 2022/11/19

-   Computation speed of `nethist()` and `violin_netsummary()` are improved.

### Bug fixes

-   `nethist()`: computation of the diagonal part o in greedy search *Rcpp* routine is fixed.
-   `nethist()`: computation of log-likelihood when there is a smaller size of group is fixed, in greedy search *Rcpp* routine.

# nethist 0.1.1

update: 2022/11/08

-   Documentations update: provide extra about default values of arguments, rewrite examples.

### Bug fixes

-   Resolving default value issue of *R* in `violin_netsummary()`.

# nethist 0.1.0

update: 2022/11/04

-   Initial release of `nethist` R package (development)

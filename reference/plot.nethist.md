# Network histogram plot

Drawing
[`lattice::levelplot()`](https://rdrr.io/pkg/lattice/man/levelplot.html)
using a `nethist` object.

## Usage

``` r
# S3 method for class 'nethist'
plot(
  x,
  type = "nethist",
  idx_order = 1:max(x$cluster),
  power = 0.25,
  col.regions = function(n) grDevices::hcl.colors(n, palette = "Reds 3", rev = TRUE),
  colorkey = FALSE,
  prob = FALSE,
  digits = 2,
  prob.cex = 0.1 + 0.5/log10(max(x$cluster)),
  prob.col = "white",
  at = NULL,
  y = NA,
  ...
)
```

## Arguments

- x:

  a nethist object from
  [`nethist()`](https://enigmasong.github.io/nethist/reference/nethist.md).

- type:

  One of `nethist` (default) or `prob`. `"pmat"` is a deprecated alias
  for `"prob"`.

- idx_order:

  A numeric vector for index label order, which must be a permutation of
  `x$cluster`. If `NA`, it uses `1:max(x$cluster)`.

- power:

  A positive number for the power transform applied to the graphon
  estimate. Only used when `type = "nethist"`. Default is `0.25`.

- col.regions:

  A function taking an integer `n` and returning `n` colors. Default is
  viridis via
  [`grDevices::hcl.colors()`](https://rdrr.io/r/grDevices/palettes.html).

- colorkey:

  Logical. Whether to draw a color legend. Default `FALSE`.

- prob:

  Logical. Whether to print block probabilities on the plot. Default
  `FALSE`.

- digits:

  Integer. Number of decimal places for probabilities.

- prob.cex:

  Numeric. `cex` for probability labels.

- prob.col:

  Color for probability labels. Default `"white"`.

- at:

  A numeric vector of breakpoints for the color scale. If `NULL`
  (default), breakpoints are computed automatically from the data range.
  Specify a fixed vector to compare plots from different fits on a
  common scale.

- y:

  A dummy variable for S3 dispatch. Never used.

- ...:

  Additional arguments passed to
  [`lattice::levelplot()`](https://rdrr.io/pkg/lattice/man/levelplot.html).

## Value

Called for its side effects (plotting). Returns `NULL` invisibly.

## Examples

``` r
# \donttest{
set.seed(2022)
A <- igraph::sample_gnp(200, 0.05)
hist_A <- nethist(A)
plot(hist_A)

plot(hist_A, power = 0.5)

plot(hist_A, type = "prob", prob = TRUE)

# }
```

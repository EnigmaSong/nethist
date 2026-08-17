# Network histogram plot

Drawing
[`lattice::levelplot()`](https://rdrr.io/pkg/lattice/man/levelplot.html)
using a `multinethist` object.

## Usage

``` r
# S3 method for class 'multinethist'
plot(
  x,
  y = NA,
  type = "nethist",
  idx_order = 1:max(x$cluster),
  power = 0.25,
  col.regions = function(n) grDevices::hcl.colors(n, palette = "Reds 3", rev = TRUE),
  colorkey = FALSE,
  prob = FALSE,
  digits = 2,
  prob.cex = 0.1 + 0.5/log10(max(x$cluster)),
  prob.col = "white",
  layout = NULL,
  layer_titles = NULL,
  at = NULL,
  ...
)
```

## Arguments

- x:

  a multinethist object from
  [`multinethist()`](https://enigmasong.github.io/nethist/reference/multinethist.md).

- y:

  A dummy variable for S3 dispatch. Never used.

- type:

  One of `nethist` (default) or `prob`. `"MNhist"` is a deprecated alias
  for `"nethist"`.

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

- layout:

  An integer vector `c(nrows, ncols)` specifying the panel grid for
  multi-layer plots, following the `mfrow` convention. If `NULL`
  (default), each layer is plotted as a separate figure.

- layer_titles:

  A character vector of length equal to the number of layers plotted,
  giving each panel's title. If `NULL` (default), titles default to
  `"Layer 1"`, `"Layer 2"`, etc.

- at:

  A numeric vector of breakpoints for the color scale. If `NULL`
  (default), breakpoints are computed automatically from the data range.
  Specify a fixed vector to compare plots from different fits on a
  common scale.

- ...:

  Additional arguments passed to
  [`lattice::levelplot()`](https://rdrr.io/pkg/lattice/man/levelplot.html).

## Value

Called for its side effects (plotting). Returns `NULL` invisibly.

## Examples

``` r
# \donttest{
set.seed(42)
data(IndianVil)
mnhist_Ind_vil <- multinethist(IndianVil)
plot(mnhist_Ind_vil)












plot(mnhist_Ind_vil, power = 0.5)












plot(mnhist_Ind_vil, layout = c(3,4))

plot(mnhist_Ind_vil, layer_titles = paste0("Network ", seq_along(mnhist_Ind_vil$rho_hat)))












# }
```

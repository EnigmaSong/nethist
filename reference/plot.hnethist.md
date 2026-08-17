# Plot an hnethist object

Plots a heatmap or a BIC curve for an `hnethist` object.

## Usage

``` r
# S3 method for class 'hnethist'
plot(x, type = "nethist", at = NULL, ...)
```

## Arguments

- x:

  an `hnethist` object from
  [`hnethist()`](https://enigmasong.github.io/nethist/reference/hnethist.md).

- type:

  One of `nethist` (default), `prob`, or `BIC`. When `type = "BIC"`,
  plots BIC values against the number of shapes `s`. The selected model
  is marked with a dashed vertical line, and the rightmost point
  (largest `s`, corresponding to the initial nethist fit) is labelled.
  `type = "bic"` is also accepted.

- at:

  A numeric vector of breakpoints for the color scale. Passed to
  [`plot.nethist()`](https://enigmasong.github.io/nethist/reference/plot.nethist.md).
  See
  [`plot.nethist()`](https://enigmasong.github.io/nethist/reference/plot.nethist.md)
  for details.

- ...:

  Additional arguments passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) when
  `type = "BIC"`, or to
  [`plot.nethist()`](https://enigmasong.github.io/nethist/reference/plot.nethist.md)
  otherwise.

## Value

Called for its side effects (plotting). Returns `NULL` invisibly.

## See also

[`plot.nethist()`](https://enigmasong.github.io/nethist/reference/plot.nethist.md),
[`hnethist()`](https://enigmasong.github.io/nethist/reference/hnethist.md)

## Examples

``` r
# \donttest{
set.seed(2022)
A <- igraph::as_adjacency_matrix(
  igraph::sample_gnp(100, 0.3), sparse = FALSE)
fit <- suppressMessages(hnethist(A))
plot(fit)

plot(fit, type = "BIC")

# }
```

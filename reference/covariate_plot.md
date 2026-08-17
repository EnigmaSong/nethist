# Bin summary by covariate

Drawing a bin summary plot of covariates given a `multinethist` object
with a user-specified order.

## Usage

``` r
covariate_plot(
  object,
  covariate,
  idx_order = 1:max(object$cluster),
  main = NA,
  xlab = NULL,
  ylab = NA,
  legend_title = NA,
  stat = "count",
  position = "stack"
)

summary_plot(object, covariate, ...)
```

## Arguments

- object:

  a `multinethist` object from
  [`multinethist()`](https://enigmasong.github.io/nethist/reference/multinethist.md).

- covariate:

  a vector for univariate covariate. If it is a factor, a stacked bar
  chart is drawn. If it is numeric, a violin plot is drawn.

- idx_order:

  A numeric vector for index label order, which must be a permutation of
  `object$cluster`. If `NA`, it uses `1:max(object$cluster)`.

- main:

  title of summary plot. If NA, the plot has no title.

- xlab:

  label of x-axis. If `NULL`, no label is shown.

- ylab:

  label of y-axis. If NA, y-axis label is "covariate"

- legend_title:

  title of legend. If NA, the legend title is "covariate"

- stat:

  variables pass to
  [`ggplot2::geom_bar()`](https://ggplot2.tidyverse.org/reference/geom_bar.html).
  Only used for a factor covariate.

- position:

  variables pass to
  [`ggplot2::geom_bar()`](https://ggplot2.tidyverse.org/reference/geom_bar.html).
  Only used for a factor covariate.

- ...:

  currently unused.

## Value

A `ggplot` object. Printed as a side effect. Returns the plot invisibly.

## Details

When `covariate` is a factor, a stacked bar chart is drawn with bins
ordered by `idx_order`. When `covariate` is numeric, a violin plot is
drawn.

## Examples

``` r
# \donttest{
set.seed(42)
data(polblog)
nethist_polblog <- multinethist(polblog)
x <- factor(c(rep("Liberal", 586), rep("Conservative", 638)))
covariate_plot(nethist_polblog, x)

# }
```

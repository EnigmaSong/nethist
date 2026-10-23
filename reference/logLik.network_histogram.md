# Log-likelihood of network histogram fits

Extracts the Bernoulli log-likelihood of a fitted network histogram. The
method is defined on the common base class, so it applies to `nethist`,
`multinethist`, and `hnethist` objects alike; see
[network_histogram-class](https://enigmasong.github.io/nethist/reference/network_histogram-class.md).

## Usage

``` r
# S3 method for class 'network_histogram'
logLik(object, normalized = FALSE, ...)
```

## Arguments

- object:

  a `network_histogram` object, that is, a fit returned by
  [`nethist`](https://enigmasong.github.io/nethist/reference/nethist.md),
  [`multinethist`](https://enigmasong.github.io/nethist/reference/multinethist.md),
  or
  [`hnethist`](https://enigmasong.github.io/nethist/reference/hnethist.md).

- normalized:

  logical. If `FALSE` (default), the unnormalized log-likelihood is
  returned as a `logLik` object. If `TRUE`, the per-edge normalized
  value stored in `object$normalized_LL` is returned as a plain numeric.

- ...:

  currently unused.

## Value

With `normalized = FALSE`, an object of class `logLik` carrying `df` and
`nobs` attributes, so that
[`stats::AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`stats::BIC()`](https://rdrr.io/r/stats/AIC.html) work directly on the
fit. With `normalized = TRUE`, an unclassed numeric equal to
`object$normalized_LL`.

## Details

The normalized value is not given class `logLik` on purpose. Passing a
per-edge normalized value to [`AIC()`](https://rdrr.io/r/stats/AIC.html)
or [`BIC()`](https://rdrr.io/r/stats/AIC.html) would return a
meaningless number without any error.

The degrees of freedom count the free block probabilities:
`k * (k + 1) / 2` for `nethist`, where `k` is the number of blocks; the
same for `multinethist` fitted with `common_f = TRUE`, and `L` times
that for `L` layers otherwise; and `s`, the number of selected shapes,
for `hnethist`.

The number of observations is the number of Bernoulli edge variables,
i.e. `n * (n - 1) / 2` for a single layer and `L` times that for `L`
layers.

The log-likelihood is evaluated at the fitted probability matrix
regardless of `method`, so it is also available for `method = "LSE"`
fits.

## See also

[network_histogram-class](https://enigmasong.github.io/nethist/reference/network_histogram-class.md),
[`summary.network_histogram`](https://enigmasong.github.io/nethist/reference/summary.network_histogram.md)

## Examples

``` r
# \donttest{
set.seed(42)
data(polblog)
fit <- suppressMessages(nethist(polblog))
logLik(fit)
#> 'log Lik.' -49642.94 (df=136)
logLik(fit, normalized = TRUE)
#> [1] -2.969964
AIC(fit)
#> [1] 99557.89
BIC(fit)
#> [1] 101125.4
# }
```

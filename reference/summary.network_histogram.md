# Summarize network histogram fits

Detailed summary of a fitted network histogram. The `print` methods of
the fitted objects give a short identification only; the model fit
statistics are reported here.

## Usage

``` r
# S3 method for class 'network_histogram'
summary(object, ...)

# S3 method for class 'multinethist'
summary(object, ...)

# S3 method for class 'hnethist'
summary(object, ...)

# S3 method for class 'summary.network_histogram'
print(x, digits = getOption("digits"), ...)
```

## Arguments

- object:

  a `network_histogram` object, that is, a fit returned by
  [`nethist`](https://enigmasong.github.io/nethist/reference/nethist.md),
  [`multinethist`](https://enigmasong.github.io/nethist/reference/multinethist.md),
  or
  [`hnethist`](https://enigmasong.github.io/nethist/reference/hnethist.md).

- ...:

  additional arguments. Currently unused.

- x:

  a `summary.network_histogram` object.

- digits:

  number of significant digits used when printing numeric quantities.

## Value

[`summary()`](https://rdrr.io/r/base/summary.html) returns an object
inheriting from `summary.network_histogram`: a list holding the fitted
object together with the number of nodes, the number of blocks, the
log-likelihood, the degrees of freedom, `AIC`, and `BIC`. Multilayer
fits add `n_layers`, and `hnethist` fits add the selected number of
shapes `s` and the number of candidate models `n_candidate`. The `print`
method returns its input invisibly.

## Details

The `summary` method is defined on the common base class, so it applies
to `nethist`, `multinethist`, and `hnethist` objects alike. The subclass
methods add the fields specific to them and are otherwise deferred to
the base method via
[`NextMethod()`](https://rdrr.io/r/base/UseMethod.html); see
[network_histogram-class](https://enigmasong.github.io/nethist/reference/network_histogram-class.md).

The estimated probabilities themselves are not printed, since the matrix
grows with the number of blocks. Use `object$thetahat` or
[`fitted()`](https://enigmasong.github.io/nethist/reference/fitted.nethist.md)
to obtain them.

## See also

[network_histogram-class](https://enigmasong.github.io/nethist/reference/network_histogram-class.md),
[`print.nethist`](https://enigmasong.github.io/nethist/reference/print.nethist.md),
[`logLik.network_histogram`](https://enigmasong.github.io/nethist/reference/logLik.network_histogram.md)

## Examples

``` r
# \donttest{
set.seed(42)
data(polblog)
fit <- suppressMessages(nethist(polblog))
fit
#> 
#> Network histogram fit (nethist)
#> 
#> Nodes: 1224   Blocks: 16   Bandwidth: 73
#> Method: Profile Likelihood
#> thetahat: 16 x 16 probability matrix
#> Normalized log-likelihood: -3.052758
#> 
#> Use summary() for details.
summary(fit)
#> 
#> Network histogram fit (nethist)
#> 
#> Nodes: 1224   Blocks: 16   Bandwidth: 73
#> Method: Profile Likelihood
#> Edge density (rho_hat): 0.02233205
#> 
#> Model fit:
#>   Normalized log-likelihood: -3.052758
#>   Log-likelihood: -51026.85  (df = 136)
#>   AIC: 102325.7   BIC: 103893.2
#>   Mean squared error: 0.01814148
#> 
#> thetahat: 16 x 16 probability matrix (use fitted() or $thetahat)
# }
```

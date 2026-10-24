# Print nethist objects

Compact one-screen description of `nethist`, `multinethist`, and
`hnethist` objects. The estimated probabilities are reported by their
dimensions only; use
[`summary.network_histogram`](https://enigmasong.github.io/nethist/reference/summary.network_histogram.md)
for the model fit statistics and `object$thetahat` for the values
themselves.

## Usage

``` r
# S3 method for class 'nethist'
print(x, digits = getOption("digits"), ...)

# S3 method for class 'multinethist'
print(x, digits = getOption("digits"), ...)

# S3 method for class 'hnethist'
print(x, digits = getOption("digits"), ...)
```

## Arguments

- x:

  a `nethist`, `multinethist`, or `hnethist` object.

- digits:

  number of significant digits used when printing the loss.

- ...:

  additional arguments. Currently unused.

## Value

The input object, invisibly.

## See also

[`summary.network_histogram`](https://enigmasong.github.io/nethist/reference/summary.network_histogram.md)

## Examples

``` r
# \donttest{
set.seed(42)
data(polblog)
fit <- suppressMessages(nethist(polblog))
print(fit)
#> 
#> Network histogram fit (nethist)
#> 
#> Nodes: 1224   Blocks: 16   Bandwidth: 73
#> Method: Profile Likelihood
#> thetahat: 16 x 16 probability matrix
#> Normalized log-likelihood: -2.969964
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
#>   Normalized log-likelihood: -2.969964
#>   Log-likelihood: -49642.94  (df = 136)
#>   AIC: 99557.89   BIC: 101125.4
#>   Mean squared error: 0.01769046
#> 
#> thetahat: 16 x 16 probability matrix (use fitted() or $thetahat)
# }
```

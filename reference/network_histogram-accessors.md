# Extract components of a network histogram fit

Accessor functions for fitted network histograms. They return the
quantities stored in a fit, so that code does not rely on the internal
structure of the object, which may change in future versions.

## Usage

``` r
block_id(object, ...)

# S3 method for class 'network_histogram'
block_id(object, ...)

bandwidth(object, ...)

# S3 method for class 'network_histogram'
bandwidth(object, ...)

shape_id(object, ...)

# Default S3 method
shape_id(object, ...)

# S3 method for class 'hnethist'
shape_id(object, ...)
```

## Arguments

- object:

  a fit returned by
  [`nethist`](https://enigmasong.github.io/nethist/reference/nethist.md),
  [`multinethist`](https://enigmasong.github.io/nethist/reference/multinethist.md),
  or
  [`hnethist`](https://enigmasong.github.io/nethist/reference/hnethist.md).
  For `shape_id()`, a fit returned by
  [`hnethist`](https://enigmasong.github.io/nethist/reference/hnethist.md).

- ...:

  currently unused.

## Value

`block_id()` returns an integer vector of length \\n\\ with the block
label of each vertex. `bandwidth()` returns a single number.
`shape_id()` returns an integer \\k \times k\\ matrix whose entry \\(a,
b)\\ is the shape of the block formed by groups \\a\\ and \\b\\.

## Details

`block_id()` returns the block assignment of each vertex. `bandwidth()`
returns the bandwidth, that is, the group size used to construct the
network histogram. `shape_id()` returns, for a hybrid network histogram,
the shape assigned to each block as a symmetric \\k \times k\\ matrix,
where \\k\\ is the number of blocks. Shapes are numbered in increasing
order of their estimated connection probability, so shape 1 has the
lowest probability.

The estimated block probabilities are returned by
[`coef()`](https://enigmasong.github.io/nethist/reference/coef.nethist.md).

## See also

[`coef.nethist`](https://enigmasong.github.io/nethist/reference/coef.nethist.md),
[network_histogram-class](https://enigmasong.github.io/nethist/reference/network_histogram-class.md)

## Examples

``` r
# \donttest{
set.seed(42)
data(polblog)
fit <- hnethist(polblog)
bandwidth(fit)
#> [1] 73
table(block_id(fit))
#> 
#>   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16 
#>  73  73  73  73  73  73  73  73  73  73  73  73  73  73  73 129 
shape_id(fit)
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
#>  [1,]   14    9   12    7   11   11   10    7    5     9     8     7     3
#>  [2,]    9   15    7   13    6    6    4   11   10     4     1     1     6
#>  [3,]   12    7   10    4   10    7    8    4    2     7     4     1     1
#>  [4,]    7   13    4   11    4    1    2   10    8     3     1     1     7
#>  [5,]   11    6   10    4    9    7    8    4    3     6     5     1     1
#>  [6,]   11    6    7    1    7    4    5    2    1     4     2     1     1
#>  [7,]   10    4    8    2    8    5   10    3    2     5     6     2     1
#>  [8,]    7   11    4   10    4    2    3    9    7     3     1     1     7
#>  [9,]    5   10    2    8    3    1    2    7    4     1     1     1     4
#> [10,]    9    4    7    3    6    4    5    3    1     3     2     1     1
#> [11,]    8    1    4    1    5    2    6    1    1     2     1     1     1
#> [12,]    7    1    1    1    1    1    2    1    1     1     1     1     1
#> [13,]    3    6    1    7    1    1    1    7    4     1     1     1     5
#> [14,]    4    8    2    5    2    1    1    5    3     1     1     1     3
#> [15,]    3    6    2    2    2    1    4    1    2     1     1     1     1
#> [16,]    2    1    3    3    4    1    2    2    2     1     1     1     2
#>       [,14] [,15] [,16]
#>  [1,]     4     3     2
#>  [2,]     8     6     1
#>  [3,]     2     2     3
#>  [4,]     5     2     3
#>  [5,]     2     2     4
#>  [6,]     1     1     1
#>  [7,]     1     4     2
#>  [8,]     5     1     2
#>  [9,]     3     2     2
#> [10,]     1     1     1
#> [11,]     1     1     1
#> [12,]     1     1     1
#> [13,]     3     1     2
#> [14,]     1     2     1
#> [15,]     2     1     1
#> [16,]     1     1     2
# }
```

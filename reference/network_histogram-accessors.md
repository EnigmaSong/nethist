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
#>  [1,]   19   13   18   11   17   16   15    9   14     6    11     9     3
#>  [2,]   13   19   10   18    9    7    5   17    7    15     4     1    10
#>  [3,]   18   10   15    6   14   10   12    6    9     3     6     1     2
#>  [4,]   11   18    6   16    6    2    3   14    5    11     4     1     5
#>  [5,]   17    9   14    6   14   11   12    6    8     4     7     1     3
#>  [6,]   16    7   10    2   11    6    7    3    6     1     3     1     1
#>  [7,]   15    5   12    3   12    7   14    4    7     1     7     2     1
#>  [8,]    9   17    6   14    6    3    4   13    5     9     2     1     3
#>  [9,]   14    7    9    5    8    6    7    5    5     1     2     1     1
#> [10,]    6   15    3   11    4    1    1    9    1     6     1     1     3
#> [11,]   11    4    6    4    7    3    7    2    2     1     4     1     1
#> [12,]    9    1    1    1    1    1    2    1    1     1     1     1     1
#> [13,]    3   10    2    5    3    1    1    3    1     3     1     1     3
#> [14,]    6   12    2   10    3    1    1   10    1     4     2     2     2
#> [15,]    6    2    6    1    8    2    5    1    3     1     2     1     1
#> [16,]    2    5    1    6    1    1    2    7    2     3     1     1     2
#>       [,14] [,15] [,16]
#>  [1,]     6     6     2
#>  [2,]    12     2     5
#>  [3,]     2     6     1
#>  [4,]    10     1     6
#>  [5,]     3     8     1
#>  [6,]     1     2     1
#>  [7,]     1     5     2
#>  [8,]    10     1     7
#>  [9,]     1     3     2
#> [10,]     4     1     3
#> [11,]     2     2     1
#> [12,]     2     1     1
#> [13,]     2     1     2
#> [14,]     7     2     5
#> [15,]     2     1     1
#> [16,]     5     1     4
# }
```

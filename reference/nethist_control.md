# Control parameters for network histogram algorithms

Constructs a control object for
[`nethist`](https://enigmasong.github.io/nethist/reference/nethist.md),
[`multinethist`](https://enigmasong.github.io/nethist/reference/multinethist.md),
and
[`hnethist`](https://enigmasong.github.io/nethist/reference/hnethist.md).

## Usage

``` r
nethist_control(
  algorithm = "greedy",
  max_itr = 5e+06,
  greedy_swap_rule = "single_random",
  greedy_stop_threshold = 20000,
  verbose = FALSE,
  ...
)

# S3 method for class 'nethist_control'
print(x, ...)
```

## Arguments

- algorithm:

  character. Optimization algorithm. Currently only `"greedy"` is
  implemented.

- max_itr:

  integer. Maximum number of iterations. Default is \\5 \times 10^6\\.

- greedy_swap_rule:

  character. Vertex-pair selection rule for the greedy search. At each
  iteration, two vertices are drawn and their group labels are swapped
  if the move improves the objective. Currently only `"single_random"`
  (one pair drawn uniformly at random) is implemented.

- greedy_stop_threshold:

  integer. Early stopping criterion for the greedy search; the algorithm
  terminates if the objective has not improved for this many consecutive
  iterations. Default is 20,000.

- verbose:

  logical. Print progress messages during fitting.

- ...:

  Accepts deprecated argument names `swap_rule` and
  `consecutive_iter_threshold` with a warning.

- x:

  a `nethist_control` object.

## Value

An object of class `"nethist_control"`.

## See also

[`nethist`](https://enigmasong.github.io/nethist/reference/nethist.md),
[`multinethist`](https://enigmasong.github.io/nethist/reference/multinethist.md),
[`hnethist`](https://enigmasong.github.io/nethist/reference/hnethist.md)

## Examples

``` r
# default control object
ctrl <- nethist_control()
print(ctrl)
#> nethist control parameters:
#>   algorithm            : greedy 
#>   max_itr              : 5e+06 
#>   greedy_swap_rule     : single_random 
#>   greedy_stop_threshold: 20000 
#>   verbose              : FALSE 

# reduce iteration limit for quick testing
ctrl <- nethist_control(max_itr = 1e4, greedy_stop_threshold = 100)

# \donttest{
data(polblog)
fit <- nethist(polblog, control = nethist_control(max_itr = 1e4))
# }
```

# Political blog

A pre-processed version of the Political Blog dataset from Olhede and
Wolfe (2014), based on the original data by Adamic and Glance (2005).
The provided igraph object represents the graph after removing
zero-degree nodes and simplifying the graph structure. The graph is
derived from the edge list `polblog_edgelist`.

## Usage

``` r
data(polblog)
```

## Format

`polblog` is an igraph object, which is pre-processed data symmetrize
and simplified without self-loops.

## Details

The adjancency matrix from `polblog_edgelist` is sparse, assymetric, so
it needs to be symmetrized if a method requires a symmetric matrix.
Nodes with zero degree need to be removed. All edges are considered as
undirected edges. Then, we get the undirected version of dataset as
`polblog`.

First 586 blogs are liberal; remaining 638 are conservative.

## Source

http://www-personal.umich.edu/~mejn/netdata/

https://github.com/p-wolfe/network-histogram-code

## References

Adamic, L. A., & Glance, N. (2005). The political blogosphere and the
2004 US election: divided they blog. In Proceedings of the 3rd
international workshop on Link discovery (pp. 36-43)

## Examples

``` r
data(polblog_edgelist)
data(polblog)

#From polblog_edgelist to polblog
G<- igraph::graph_from_edgelist(as.matrix(polblog_edgelist), directed = FALSE)
G<- igraph::delete.vertices(G, igraph::degree(G) == 0)
#> Warning: `delete.vertices()` was deprecated in igraph 2.0.0.
#> ℹ Please use `delete_vertices()` instead.
G<- igraph::simplify(G)

polblog
#> IGRAPH c8f356f U--- 1224 16715 -- 
#> + edges from c8f356f:
#>  [1] 1--   2 1--  19 1--  21 1--  48 1--  57 1--  72 1--  97 1-- 127 1-- 156
#> [10] 1-- 193 1-- 204 1-- 205 1-- 254 1-- 287 1-- 341 1-- 377 1-- 390 1-- 446
#> [19] 1-- 451 1-- 498 1-- 499 1-- 501 1-- 518 1-- 909 1--1008 1--1178 2--  16
#> [28] 2--  22 2--  26 2--  46 2--  48 2--  84 2-- 119 2-- 126 2-- 127 2-- 131
#> [37] 2-- 141 2-- 147 2-- 166 2-- 191 2-- 198 2-- 254 2-- 261 2-- 284 2-- 304
#> [46] 2-- 306 2-- 339 2-- 352 2-- 373 2-- 377 2-- 387 2-- 428 2-- 429 2-- 440
#> [55] 2-- 442 2-- 443 2-- 444 2-- 448 2-- 487 2-- 498 2-- 499 2-- 501 2-- 507
#> [64] 2-- 514 2-- 518 2-- 522 2-- 537 2-- 584 2-- 586 2-- 909 3-- 405 3-- 562
#> [73] 3-- 563 3--1180 4-- 574 5-- 359 6--  37 6--  48 6-- 147 6-- 148 6-- 164
#> + ... omitted several edges
```

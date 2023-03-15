# nethist 0.2.3

update: 2023/03/15
 
-  Documentation updates
-  Internal functions (not visible to users)

# nethist 0.2.2

update: 2022/12/08

-   `plot.nethist()`
    -   You can display the estimated probabilities in `p_mat` of `nethist` objects by setting `prob=TRUE`. Use `prob.col` and `prob.cex` to set the color and size of text, respectively.
-   `nethist()`
    -   Check the convergence criterion every 5 iterations from the fifth iteration (5,10,15,...) instead of every 5 iterations from the first iteration (1,6,11,...).
    -   minor speed improvement.

# nethist 0.2.1

update: 2022/11/26

-   Initialization step of `nethist()` is faster than earlier by computing hamming distance from the (simple) adjacency matrix instead of Manhattan distance.

# nethist 0.2.0

update: 2022/11/19

-   Computation speed of `nethist()` and `violin_netsummary()` are improved.

### Bug fixes

-   `nethist()`: computation of the diagonal part o in greedy search *Rcpp* routine is fixed.
-   `nethist()`: computation of log-likelihood when there is a smaller size of group is fixed, in greedy search *Rcpp* routine.

# nethist 0.1.1

update: 2022/11/08

-   Documentations update: provide extra about default values of arguments, rewrite examples.

### Bug fixes

-   Resolving default value issue of *R* in `violin_netsummary()`.

# nethist 0.1.0

update: 2022/11/04

-   Initial release of `nethist` R package (development)

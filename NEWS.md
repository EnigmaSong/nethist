# nethist 0.2.2

update: 2022/12/31


* minor changes not affected to users (e.g test codes).

# nethist 0.2.1

update: 2022/11/26

* Initialization step of `nethist()` is faster than earlier by computing hamming distance from the (simple) adjacency matrix instead of Manhattan distance.


# nethist 0.2.0 

update: 2022/11/19

*   Computation speed of `nethist()` and `violin_netsummary()` are improved. 

### Bug fixes

*   `nethist()`: computation of the diagonal part o in greedy search  *Rcpp* routine is fixed.
*   `nethist()`: computation of log-likelihood when there is a smaller size of group is fixed, in greedy search  *Rcpp* routine.

# nethist 0.1.1 

update: 2022/11/08

*   Documentations update: provide extra about default values of arguments, rewrite examples.

### Bug fixes

*   Resolving default value issue of *R* in `violin_netsummary()`.


# nethist 0.1.0 

update: 2022/11/04

*   Initial release of nethist R package
   

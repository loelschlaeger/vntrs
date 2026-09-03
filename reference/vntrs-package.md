# vntrs: Variable Neighborhood Trust Region Search

Implements a variable neighborhood trust region search (VNTRS) algorithm
for nonlinear global optimization, based on Bierlaire et al. (2009) "A
Heuristic for Nonlinear Global Optimization"
[doi:10.1287/ijoc.1090.0343](https://doi.org/10.1287/ijoc.1090.0343) .
The method combines neighborhood exploration with a trust-region
framework to search the solution space efficiently. It can terminate a
local search early when the iterates converge toward a previously
visited local optimum or when further improvement within the current
region is unlikely. The algorithm can also be used to identify multiple
local optima.

## See also

Useful links:

- <https://loelschlaeger.de/vntrs/>

- Report bugs at <https://github.com/loelschlaeger/vntrs/issues>

## Author

**Maintainer**: Lennart Oelschläger <oelschlaeger.lennart@gmail.com>
([ORCID](https://orcid.org/0000-0001-5421-9313))

Authors:

- Lennart Oelschläger <oelschlaeger.lennart@gmail.com>
  ([ORCID](https://orcid.org/0000-0001-5421-9313))

## Examples

``` r
### Example 1: Rosenbrock function; one global minimum; no local minima
set.seed(1)
rosenbrock <- function(x) 100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
vntrs(f = rosenbrock, npar = 2)
#>   p1 p2        value global
#> 1  1  1 1.681733e-16   TRUE

### Example 2: Six-hump camel function; two global minima; four local minima
set.seed(1)
camel <- function(x) {
  (4 - 2.1*x[1]^2 + x[1]^4/3) * x[1]^2 + x[1]*x[2] + (-4 + 4*x[2]^2) * x[2]^2
}
vntrs(
  f = camel, npar = 2,
  lower = c(-3, -2), upper = c(3, 2), # search bounds
  collect_all = TRUE,                 # also collect local optima
  neighborhoods = 10                  # number of neighborhoods
)
#>            p1         p2      value global
#> 1 -0.08984203  0.7126564 -1.0316285   TRUE
#> 2  0.08984201 -0.7126564 -1.0316285   TRUE
#> 3  1.70360672 -0.7960836 -0.2154638  FALSE
#> 4 -1.70360674  0.7960836 -0.2154638  FALSE
#> 5 -1.60710473 -0.5686515  2.1042503  FALSE
#> 6  1.60710475  0.5686517  2.1042503  FALSE
```

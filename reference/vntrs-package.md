# vntrs: Variable Neighborhood Trust Region Search

Implements the variable neighborhood trust region search (VNTRS)
algorithm for nonlinear global optimization, following Bierlaire et al.
(2009) "A Heuristic for Nonlinear Global Optimization"
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

- <https://github.com/loelschlaeger/vntrs/>

- Report bugs at <https://github.com/loelschlaeger/vntrs/issues>

## Author

**Maintainer**: Lennart Oelschläger <oelschlaeger.lennart@gmail.com>
([ORCID](https://orcid.org/0000-0001-5421-9313))

## Examples

``` r
rosenbrock <- function(x) 100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
vntrs(f = rosenbrock, npar = 2)
#>   p1 p2        value global
#> 1  1  1 5.916025e-20   TRUE
```


<!-- README.md is generated from README.Rmd. Please edit that file -->

# vntrs <a href="https://loelschlaeger.de/vntrs/"><img src="man/figures/logo.png" align="right" height="139" alt="vntrs website" /></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/loelschlaeger/vntrs/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/loelschlaeger/vntrs/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/loelschlaeger/vntrs/graph/badge.svg)](https://app.codecov.io/gh/loelschlaeger/vntrs)
[![CRAN
status](https://www.r-pkg.org/badges/version/vntrs)](https://CRAN.R-project.org/package=vntrs)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/grand-total/vntrs)](https://CRAN.R-project.org/package=vntrs)
<!-- badges: end -->

This R package implements the variable neighborhood trust region search
(VNTRS) algorithm for nonlinear global optimization, following Bierlaire
et al. (2009) “A Heuristic for Nonlinear Global Optimization”. The
method combines neighborhood exploration with a trust-region framework
to search the solution space efficiently. It can terminate a local
search early when the iterates converge toward a previously visited
local optimum or when further improvement within the current region is
unlikely. The algorithm can also be used to identify multiple local
optima.

## Installation

You can install the released package version from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("vntrs")
```

## How to get started

1.  Specify a function `f` that computes the objective value. It may
    also return the gradient and Hessian. Return either a numeric value
    or a named list with `value`, `gradient`, and `hessian` components.
    Omitted derivatives are approximated by finite differences.

2.  Call `vntrs::vntrs(f = f, npar = npar, minimize = minimize)`, where

- `npar` is the number of parameters of `f` and

- `minimize` determines whether `f` should be minimized
  (`minimize = TRUE`, the default) or maximized (`minimize = FALSE`).

Optionally, the algorithm can be tuned by setting dedicated control
arguments (for example `init_runs` or `neighbors`, see `help("vntrs")`
for details).

## Example

The example below minimizes the well-known [six-hump camel
function](https://www.sfu.ca/~ssurjano/camel6.html), which has two
global minima.

``` r
set.seed(1)

camel <- function(x) {
  x1 <- x[1]
  x2 <- x[2]
  value <- (4 - 2.1 * x1^2 + x1^4 / 3) * x1^2 + x1 * x2 + (-4 + 4 * x2^2) * x2^2
  gradient <- c(
    8 * x1 - 8.4 * x1^3 + 2 * x1^5 + x2,
    x1 - 8 * x2 + 16 * x2^3
  )
  hessian <- matrix(
    c(
      8 - 25.2 * x1^2 + 10 * x1^4, 1,
      1, -8 + 48 * x2^2
    ),
    nrow = 2,
    byrow = TRUE
  )
  list(value = value, gradient = gradient, hessian = hessian)
}

vntrs::vntrs(
  f = camel,           # objective supplies value, gradient, Hessian
  npar = 2,            # two variables (x1 and x2)
  init_runs = 5,       # start from 5 random points
  neighborhoods = 5,   # try 5 neighborhood radii per trust region
  neighbors = 5,       # evaluate 5 trial points per neighborhood
  lower = c(-3, -2),   # lower search bounds for x1 and x2
  upper = c(3, 2),     # upper search bounds for x1 and x2
  collect_all = TRUE,  # also look for local optima
  quiet = FALSE        # show status messages
)
#> Initialize VNTRS.
#> * Apply local search at 5 random starting points.
#> ** Run 1 [0 s] [found optimum] [optimum is unknown]
#> ** Run 2 [0 s] [found optimum]
#> ** Run 3 [0 s] [found optimum]
#> ** Run 4 [0 s] [found optimum]
#> ** Run 5 [0 s] [found optimum]
#> Start VNTRS.
#> * Select neighborhood 1.
#> ** Neighbor 1 [0 s]
#> ** Neighbor 2 [0 s] [found optimum] [optimum is unknown]
#> ** Neighbor 3 [0 s] [found optimum]
#> ** Neighbor 4 [0 s]
#> ** Neighbor 5 [0 s]
#> * Reset neighborhood, because better optimum was found.
#> * Select neighborhood 1.
#> ** Neighbor 1 [0 s]
#> ** Neighbor 2 [0 s] [found optimum]
#> ** Neighbor 3 [0 s] [found optimum]
#> ** Neighbor 4 [0 s] [found optimum]
#> ** Neighbor 5 [0 s] [found optimum]
#> * Select neighborhood 2.
#> ** Neighbor 1 [0 s] [found optimum]
#> ** Neighbor 2 [0 s]
#> ** Neighbor 3 [0 s]
#> ** Neighbor 4 [0 s] [found optimum]
#> ** Neighbor 5 [0 s] [found optimum]
#> * Select neighborhood 3.
#> ** Neighbor 1 [0 s] [found optimum] [optimum is unknown]
#> ** Neighbor 2 [0 s] [found optimum]
#> ** Neighbor 3 [0 s] [found optimum] [optimum is unknown]
#> ** Neighbor 4 [0 s] [found optimum]
#> ** Neighbor 5 [0 s] [found optimum]
#> * Select neighborhood 4.
#> ** Neighbor 1 [0 s] [found optimum]
#> ** Neighbor 2 [0 s] [found optimum]
#> ** Neighbor 3 [0 s] [found optimum]
#> ** Neighbor 4 [0 s] [found optimum]
#> ** Neighbor 5 [0 s] [found optimum]
#> * Select neighborhood 5.
#> ** Neighbor 1 [0 s] [found optimum] [optimum is unknown]
#> ** Neighbor 2 [0 s] [found optimum]
#> ** Neighbor 3 [0 s] [found optimum]
#> ** Neighbor 4 [0 s] [found optimum]
#> ** Neighbor 5 [0 s] [found optimum]
#> Done.
#>            p1         p2      value global
#> 1  0.08984201 -0.7126564 -1.0316285   TRUE
#> 2 -0.08984201  0.7126564 -1.0316285   TRUE
#> 3 -1.70360672  0.7960836 -0.2154638  FALSE
#> 4  1.60710475  0.5686515  2.1042503  FALSE
#> 5  1.70360671 -0.7960836 -0.2154638  FALSE
```

# Variable neighborhood trust region search

Run the variable neighborhood trust region search algorithm. Features:

- expanding-neighborhood exploration,

- local trust-region optimization,

- optional parameter bounds,

- analytical derivatives are not required,

- early stopping of inferior or previously explored search paths,

- collection of multiple local optima.

## Usage

``` r
vntrs(
  f,
  npar,
  lower = NULL,
  upper = NULL,
  minimize = TRUE,
  collect_all = FALSE,
  init_runs = 5L,
  init_min = -1,
  init_max = 1,
  init_iterlim = 20L,
  neighborhoods = 5L,
  neighbors = 5L,
  beta = 0.05,
  iterlim = 100L,
  identical_tolerance = 0.001,
  inferior_tolerance = 1e-06,
  interruption_gradient_tolerance = 0.001,
  gradient_tolerance = 1e-06,
  time_limit = NULL,
  quiet = TRUE
)
```

## Arguments

- f:

  \[`function`\]  
  A function that accepts a `numeric` parameter vector and returns
  either

  - a `numeric` objective value, or

  - a `list` with `value` and optional `gradient` and `hessian`
    components.

  Missing derivatives are approximated by finite differences.

  **Recommendation:** Provide analytical derivatives when available.

- npar:

  \[`integer(1)`\]  
  The length of the parameter vector passed to `f`.

- lower, upper:

  \[`numeric(npar)` \| `NULL`\]  
  Lower and upper bounds for the parameters. `NULL` represents no bound;
  use `-Inf` or `Inf` for individual unbounded parameters. Starting
  points and neighborhood points are restricted to these bounds.

- minimize:

  \[`logical(1)`\]  
  If `TRUE`, minimize `f`; if `FALSE`, maximize it.

- collect_all:

  \[`logical(1)`\]  
  If `TRUE`, do not interrupt a local search solely because it appears
  to approach an optimum inferior to the best one found so far. This
  increases the chance of returning multiple local optima, but usually
  requires more function evaluations. Duplicate optima are still merged.

  **Recommendation:** Use `FALSE` when only the best solution matters
  and `TRUE` when local optima are themselves of interest.

- init_runs:

  \[`integer(1)`\]  
  Number of random starting points used to initialize the search. A
  short local search is run from each point, and the best result seeds
  the neighborhood search. More runs improve exploration but increase
  cost approximately linearly.

  **Recommendation:** Start with `5`; increase it for strongly
  multimodal problems.

- init_min, init_max:

  \[`numeric(1)`\]  
  For a parameter without two finite bounds `lower` and `upper`, an
  initial value is uniformly sampled between `init_min` and `init_max`
  and then restricted by any one-sided bound.

  If a parameter has finite `lower` and `upper` bounds, it is instead
  sampled uniformly across that interval.

  **Recommendation:** Use `-1` and `1` for reasonably scaled unbounded
  parameters; otherwise choose a range containing plausible solutions.

- init_iterlim:

  \[`integer(1)`\]  
  Maximum number of trust-region iterations for each local search
  started from one of the `init_runs` random points.

  **Recommendation:** Start with `20`; increase it when initialization
  searches stop before reaching an optimum.

- neighborhoods:

  \[`integer(1)`\]  
  Number of neighborhoods tried around the best solution found so far.
  Neighborhood \\k\\ has scale \\1.5^{k-1}\\. The algorithm tries at
  most `neighborhoods` sizes without improvement. Whenever a better
  solution is found, it starts again with the first neighborhood.

  **Recommendation:** Start with `5`; use more for broader global
  exploration at additional computational cost.

- neighbors:

  \[`integer(1)`\]  
  Number of new starting points generated for each neighborhood. Each of
  these trial points is a parameter vector obtained by moving away from
  the best solution found so far in a sampled direction. A separate
  local trust-region search starts from every trial point. Trying more
  points makes it more likely to reach different regions of the
  parameter space, but also increases the number of local searches.

  **Recommendation:** Start with `5`; increase it when broader
  exploration is worth the additional computation time.

- beta:

  \[`numeric(1)`\]  
  Controls which directions are more likely to be selected when
  generating trial points. The local model describes the curvature
  around the best solution through eigenvectors (directions) and
  eigenvalues (strength of curvature). The sampling weight of a
  direction is proportional to `exp(beta * lambda / d)`, where `lambda`
  is its eigenvalue and `d` is the current neighborhood scale. With
  `beta = 0`, all positive and negative eigenvector directions are
  equally likely. Larger values concentrate the search on directions in
  which the objective bends more strongly.

  This is the weighting rule for \\\beta\\ in Bierlaire et al. (2009),
  Section 2.3, Equation (17).

  **Recommendation:** Start with `0.05`; use `0` to disable the
  curvature preference, for example when the objective is noisy.

- iterlim:

  \[`integer(1)`\]  
  Maximum number of trust-region iterations for each local search during
  neighborhood exploration.

  **Recommendation:** Start with `100`; increase it if otherwise
  promising local searches frequently stop before convergence.

- identical_tolerance:

  \[`numeric(1)`\]  
  Relative parameter tolerance for deciding whether two solutions
  represent the same optimum. For parameter vectors \\x\\ and \\y\\, the
  algorithm computes the difference in each component relative to
  \\\max(1, \|x_i\|, \|y_i\|)\\. The points are treated as the same
  optimum when the largest of these scaled differences does not exceed
  `identical_tolerance`. The later point is then not stored as a
  separate optimum, and a local search approaching such a point may be
  stopped early. This rule also applies when `collect_all = TRUE`.

  **Recommendation:** Start with `1e-3`; increase it when repeated
  searches return slightly different parameter vectors for what is
  clearly the same optimum; decrease it when genuinely different optima
  may lie close together.

- inferior_tolerance:

  \[`numeric(1)`\]  
  Controls when the algorithm gives up early on a local search that is
  unlikely to improve the best known solution. It compares objective
  values in absolute units and is used only after at least one optimum
  has been found. The search is stopped in either of the following
  situations:

  - Its value is worse than the best known value by more than
    `inferior_tolerance`.

  - It is within Euclidean distance `1` of a known optimum but does not
    improve the best known value by more than `inferior_tolerance`.

  Setting `collect_all = TRUE` disables both objective-value rules
  above, allowing searches toward inferior local optima to continue.
  Searches that approach an already identified optimum can still be
  stopped using `identical_tolerance`.

  **Recommendation:** Start with `1e-6`; decrease it when improvements
  below `1e-6` are meaningful.

- interruption_gradient_tolerance:

  \[`numeric(1)`\]  
  Non-negative threshold for the ordinary Euclidean gradient norm used
  by the premature-interruption rule. Once an optimum is known,
  objective-value checks controlled by `inferior_tolerance` are
  considered only when the current gradient norm is no greater than this
  value. This indicates that the local search is approaching a
  stationary point and is therefore unlikely to leave the current
  region.

  Setting `collect_all = TRUE` disables the interruption rules, so
  `interruption_gradient_tolerance` then has no effect.

  **Recommendation:** Use `1e-3`, the threshold used for premature
  interruption in Bierlaire et al. (2009). Use a smaller value to make
  premature interruption more conservative.

- gradient_tolerance:

  \[`numeric(1)`\]  
  First-order convergence tolerance for each local trust-region search.
  A local search satisfies the first-order convergence condition when
  the (projected; to deal with parameter bounds) gradient is no greater
  than `gradient_tolerance`. Additional curvature checks are used to
  avoid accepting saddle points.

  **Recommendation:** Start with `1e-6`; increase it for noisy
  objectives or when faster, less precise local searches are sufficient;
  decrease it for smooth objectives with reliable derivatives when
  greater local accuracy is required.

- time_limit:

  \[`numeric(1)` \| `NULL`\]  
  Optional approximate time limit in seconds. It is checked between
  local searches, so a running objective evaluation or local search is
  not interrupted and the elapsed time may exceed the limit.

- quiet:

  \[`logical(1)`\]  
  If `TRUE`, suppress progress messages. Warnings are still emitted.

## Value

A `data.frame` summarizing the identified optima, ordered from best to
worst by objective value, or `NULL` if none could be determined.

## References

Bierlaire et al. (2009) "A Heuristic for Nonlinear Global Optimization"
[doi:10.1287/ijoc.1090.0343](https://doi.org/10.1287/ijoc.1090.0343) .

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

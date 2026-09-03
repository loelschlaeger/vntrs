#' Variable neighborhood trust region search
#'
#' @description
#' Run the variable neighborhood trust region search algorithm. Features:
#' \itemize{
#'   \item expanding-neighborhood exploration,
#'   \item local trust-region optimization,
#'   \item optional parameter bounds,
#'   \item analytical derivatives are not required,
#'   \item early stopping of inferior or previously explored search paths,
#'   \item collection of multiple local optima.
#' }
#'
#' @references
#' Bierlaire et al. (2009) "A Heuristic for Nonlinear Global Optimization"
#' \doi{10.1287/ijoc.1090.0343}.
#'
#' @param f \[`function`\]\cr
#' A function that accepts a \code{numeric} parameter vector and returns either
#'
#' - a \code{numeric} objective value, or
#' - a \code{list} with \code{value} and optional \code{gradient} and
#'   \code{hessian} components.
#'
#' Missing derivatives are approximated by finite differences.
#'
#' \strong{Recommendation:} Provide analytical derivatives when available.
#'
#' @param npar \[`integer(1)`\]\cr
#' The length of the parameter vector passed to \code{f}.
#'
#' @param lower,upper \[`numeric(npar)` | `NULL`\]\cr
#' Lower and upper bounds for the parameters. \code{NULL} represents no bound;
#' use \code{-Inf} or \code{Inf} for individual unbounded parameters. Starting
#' points and neighborhood points are restricted to these bounds.
#'
#' @param minimize \[`logical(1)`\]\cr
#' If \code{TRUE}, minimize \code{f}; if \code{FALSE}, maximize it.
#'
#' @param collect_all \[`logical(1)`\]\cr
#' If \code{TRUE}, do not interrupt a local search solely because it appears to
#' approach an optimum inferior to the best one found so far. This increases
#' the chance of returning multiple local optima, but usually requires more
#' function evaluations. Duplicate optima are still merged.
#'
#' \strong{Recommendation:} Use \code{FALSE} when only the best solution matters
#' and \code{TRUE} when local optima are themselves of interest.
#'
#' @param init_runs \[`integer(1)`\]\cr
#' Number of random starting points used to initialize the search. A short local
#' search is run from each point, and the best result seeds the neighborhood
#' search. More runs improve exploration but increase cost approximately
#' linearly.
#'
#' \strong{Recommendation:} Start with \code{5}; increase it for strongly
#' multimodal problems.
#'
#' @param init_min,init_max \[`numeric(1)`\]\cr
#' For a parameter without two finite bounds \code{lower} and \code{upper}, an
#' initial value is uniformly sampled between \code{init_min} and
#' \code{init_max} and then restricted by any one-sided bound.
#'
#' If a parameter has finite \code{lower} and \code{upper} bounds, it is
#' instead sampled uniformly across that interval.
#'
#' \strong{Recommendation:} Use \code{-1} and \code{1} for reasonably scaled
#' unbounded parameters; otherwise choose a range containing plausible
#' solutions.
#'
#' @param init_iterlim \[`integer(1)`\]\cr
#' Maximum number of trust-region iterations for each local search started from
#' one of the \code{init_runs} random points.
#'
#' \strong{Recommendation:} Start with \code{20}; increase it when
#' initialization searches stop before reaching an optimum.
#'
#' @param neighborhoods \[`integer(1)`\]\cr
#' Number of neighborhoods tried around the best solution found so far.
#' Neighborhood \eqn{k} has scale \eqn{1.5^{k-1}}. The algorithm tries at most
#' \code{neighborhoods} sizes without improvement. Whenever a better solution is
#' found, it starts again with the first neighborhood.
#'
#' \strong{Recommendation:} Start with \code{5}; use more for broader global
#' exploration at additional computational cost.
#'
#' @param neighbors \[`integer(1)`\]\cr
#' Number of new starting points generated for each neighborhood. Each of
#' these trial points is a parameter vector obtained by moving away from the
#' best solution found so far in a sampled direction. A separate local
#' trust-region search starts from every trial point. Trying more points makes
#' it more likely to reach different regions of the parameter space, but also
#' increases the number of local searches.
#'
#' \strong{Recommendation:} Start with \code{5}; increase it when broader
#' exploration is worth the additional computation time.
#'
#' @param beta \[`numeric(1)`\]\cr
#' Controls which directions are more likely to be selected when generating
#' trial points. The local model describes the curvature around the best solution
#' through eigenvectors (directions) and eigenvalues (strength of curvature).
#' The sampling weight of a direction is proportional to
#' \code{exp(beta * lambda / d)}, where \code{lambda} is its eigenvalue and
#' \code{d} is the current neighborhood scale. With \code{beta = 0}, all
#' positive and negative eigenvector directions are equally likely. Larger
#' values concentrate the search on directions in which the objective bends
#' more strongly.
#'
#' This is the weighting rule for \eqn{\beta} in Bierlaire et al. (2009), Section
#' 2.3, Equation (17).
#'
#' \strong{Recommendation:} Start with \code{0.05}; use \code{0} to disable the
#' curvature preference, for example when the objective is noisy.
#'
#' @param iterlim \[`integer(1)`\]\cr
#' Maximum number of trust-region iterations for each local search during
#' neighborhood exploration.
#'
#' \strong{Recommendation:} Start with \code{100}; increase it if otherwise
#' promising local searches frequently stop before convergence.
#'
#' @param identical_tolerance \[`numeric(1)`\]\cr
#' Relative parameter tolerance for deciding whether two solutions represent
#' the same optimum. For parameter vectors \eqn{x} and \eqn{y}, the algorithm
#' computes the difference in each component relative to
#' \eqn{\max(1, |x_i|, |y_i|)}. The points are treated as the same optimum when
#' the largest of these scaled differences does not exceed
#' \code{identical_tolerance}. The later point is then not stored as a separate
#' optimum, and a local search approaching such a point may be stopped early.
#' This rule also applies when \code{collect_all = TRUE}.
#'
#' \strong{Recommendation:} Start with \code{1e-3}; increase it when repeated
#' searches return slightly different parameter vectors for what is clearly the
#' same optimum; decrease it when genuinely different optima may lie close
#' together.
#'
#' @param inferior_tolerance \[`numeric(1)`\]\cr
#' Controls when the algorithm gives up early on a local search that is unlikely
#' to improve the best known solution. It compares objective values in absolute
#' units and is used only after at least one optimum has been found. For
#' minimization, the search may then be stopped in either of the following
#' situations:
#' \itemize{
#'   \item Its value is worse than the best known value by more than
#'   \code{inferior_tolerance}.
#'   \item It is within Euclidean distance \code{1} of a known optimum but does
#'   not improve the best known value by more than \code{inferior_tolerance}.
#' }
#' For maximization, "worse" and "improve" have the opposite direction.
#'
#' Setting \code{collect_all = TRUE} disables both objective-value rules above,
#' allowing searches toward inferior local optima to continue. Searches that
#' approach an already identified optimum can still be stopped using
#' \code{identical_tolerance}.
#'
#' \strong{Recommendation:} Start with \code{1e-6}; decrease it when
#' improvements below \code{1e-6} are meaningful.
#'
#' @param interruption_gradient_tolerance \[`numeric(1)`\]\cr
#' Non-negative threshold for the ordinary Euclidean gradient norm used by the
#' premature-interruption rule. Once an optimum is known, objective-value checks
#' controlled by \code{inferior_tolerance} are considered only when the current
#' gradient norm is no greater than this value. This indicates that the local
#' search is approaching a stationary point and is therefore unlikely to leave
#' the current region.
#'
#' Setting \code{collect_all = TRUE} disables the interruption rules, so
#' \code{interruption_gradient_tolerance} then has no effect.
#'
#' \strong{Recommendation:} Use \code{1e-3}, the threshold used for
#' premature interruption in Bierlaire et al. (2009). Use a smaller value to
#' make premature interruption more conservative.
#'
#' @param gradient_tolerance \[`numeric(1)`\]\cr
#' First-order convergence tolerance for each local trust-region search.
#' A local search satisfies the first-order convergence condition when the
#' (projected; to deal with parameter bounds) gradient is no greater than
#' \code{gradient_tolerance}. Additional curvature checks are used to avoid
#' accepting saddle points.
#'
#' \strong{Recommendation:} Start with \code{1e-6}; increase it for noisy
#' objectives or when faster, less precise local searches are sufficient;
#' decrease it for smooth objectives with reliable derivatives
#' when greater local accuracy is required.
#'
#' @param time_limit \[`numeric(1)` | `NULL`\]\cr
#' Optional approximate time limit in seconds. It is checked between local
#' searches, so a running objective evaluation or local search is not interrupted
#' and the elapsed time may exceed the limit.
#'
#' @param quiet \[`logical(1)`\]\cr
#' If \code{TRUE}, suppress progress messages. Warnings are still emitted.
#'
#' @return
#' A \code{data.frame} summarizing the identified optima, ordered from best to
#' worst by objective value, or \code{NULL} if none could be determined.
#'
#' @export
#'
#' @examples
#' ### Example 1: Rosenbrock function; one global minimum; no local minima
#' set.seed(1)
#' rosenbrock <- function(x) 100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
#' vntrs(f = rosenbrock, npar = 2)
#'
#' ### Example 2: Six-hump camel function; two global minima; four local minima
#' set.seed(1)
#' camel <- function(x) {
#'   (4 - 2.1*x[1]^2 + x[1]^4/3) * x[1]^2 + x[1]*x[2] + (-4 + 4*x[2]^2) * x[2]^2
#' }
#' vntrs(
#'   f = camel, npar = 2,
#'   lower = c(-3, -2), upper = c(3, 2), # search bounds
#'   collect_all = TRUE,                 # also collect local optima
#'   neighborhoods = 10                  # number of neighborhoods
#' )

vntrs <- function(
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
    identical_tolerance = 1e-3,
    inferior_tolerance = 1e-6,
    interruption_gradient_tolerance = 1e-3,
    gradient_tolerance = 1e-6,
    time_limit = NULL,
    quiet = TRUE
  ) {
  oeli::input_check_response(
    check = checkmate::check_function(f),
    var_name = "f"
  )
  oeli::input_check_response(
    check = checkmate::check_count(npar, positive = TRUE),
    var_name = "npar"
  )
  npar <- as.integer(npar)
  if (is.null(lower)) {
    lower <- rep.int(-Inf, npar)
  } else {
    oeli::input_check_response(
      check = oeli::check_numeric_vector(
        lower, any.missing = FALSE, len = npar
      ),
      var_name = "lower"
    )
  }
  lower <- as.numeric(lower)
  if (is.null(upper)) {
    upper <- rep.int(Inf, npar)
  } else {
    oeli::input_check_response(
      check = oeli::check_numeric_vector(
        upper, any.missing = FALSE, len = npar
      ),
      var_name = "upper"
    )
  }
  upper <- as.numeric(upper)
  invalid_bounds <- is.finite(lower) & is.finite(upper) & lower > upper
  if (any(invalid_bounds)) {
    stop("Please ensure 'lower' <= 'upper'.", call. = FALSE)
  }
  oeli::input_check_response(
    check = checkmate::check_flag(minimize),
    var_name = "minimize"
  )
  minimize <- isTRUE(minimize)
  oeli::input_check_response(
    check = checkmate::check_flag(collect_all),
    var_name = "collect_all"
  )
  collect_all <- isTRUE(collect_all)
  oeli::input_check_response(
    check = checkmate::check_count(init_runs, positive = TRUE),
    var_name = "init_runs"
  )
  init_runs <- as.integer(init_runs)
  oeli::input_check_response(
    check = checkmate::check_number(init_min, finite = TRUE),
    var_name = "init_min"
  )
  oeli::input_check_response(
    check = checkmate::check_number(init_max, lower = init_min, finite = TRUE),
    var_name = "init_max"
  )
  oeli::input_check_response(
    check = checkmate::check_count(init_iterlim, positive = TRUE),
    var_name = "init_iterlim"
  )
  init_iterlim <- as.integer(init_iterlim)
  oeli::input_check_response(
    check = checkmate::check_count(neighborhoods, positive = TRUE),
    var_name = "neighborhoods"
  )
  neighborhoods <- as.integer(neighborhoods)
  oeli::input_check_response(
    check = checkmate::check_count(neighbors, positive = TRUE),
    var_name = "neighbors"
  )
  neighbors <- as.integer(neighbors)
  oeli::input_check_response(
    check = checkmate::check_number(beta, finite = TRUE, lower = 0),
    var_name = "beta"
  )
  oeli::input_check_response(
    check = checkmate::check_count(iterlim, positive = TRUE),
    var_name = "iterlim"
  )
  iterlim <- as.integer(iterlim)
  oeli::input_check_response(
    check = checkmate::check_number(
      identical_tolerance, finite = TRUE, lower = 0
    ),
    var_name = "identical_tolerance"
  )
  oeli::input_check_response(
    check = checkmate::check_number(
      inferior_tolerance, finite = TRUE, lower = 0
    ),
    var_name = "inferior_tolerance"
  )
  oeli::input_check_response(
    check = checkmate::check_number(
      interruption_gradient_tolerance, finite = TRUE, lower = 0
    ),
    var_name = "interruption_gradient_tolerance"
  )
  oeli::input_check_response(
    check = checkmate::check_number(
      gradient_tolerance, finite = TRUE, lower = .Machine$double.xmin
    ),
    var_name = "gradient_tolerance"
  )
  has_time_limit <- !is.null(time_limit)
  if (has_time_limit) {
    oeli::input_check_response(
      check = checkmate::check_number(time_limit, finite = TRUE, lower = 0),
      var_name = "time_limit"
    )
    if (time_limit <= 0) {
      stop("Please ensure 'time_limit' is positive.", call. = FALSE)
    }
    time_limit <- as.numeric(time_limit)
  } else {
    time_limit <- 0
  }
  oeli::input_check_response(
    check = checkmate::check_flag(quiet),
    var_name = "quiet"
  )
  quiet <- isTRUE(quiet)
  .Call(
    `_vntrs_vntrs_cpp`, f, npar, minimize, init_runs, init_min, init_max,
    init_iterlim, neighborhoods, neighbors, beta, iterlim, identical_tolerance,
    inferior_tolerance, interruption_gradient_tolerance, has_time_limit,
    time_limit, lower, upper, quiet, collect_all, gradient_tolerance
  )
}

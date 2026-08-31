#' Variable neighborhood trust region search
#'
#' @description
#' Run the variable neighborhood trust region search algorithm.
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
#' @param npar \[`integer(1)`\]\cr
#' The number of parameters of \code{f}.
#'
#' @param minimize \[`logical(1)`\]\cr
#' If \code{TRUE}, minimize \code{f}; otherwise, maximize it.
#'
#' @param init_runs \[`integer(1)`\]\cr
#' Number of random starting points for the initialization stage.
#'
#' @param init_min,init_max \[`numeric(1)`\]\cr
#' Lower and upper bounds for the uniform initialization range.
#'
#' @param init_iterlim \[`integer(1)`\]\cr
#' Maximum trust-region iterations during initialization.
#'
#' @param neighborhoods \[`integer(1)`\]\cr
#' Number of neighborhood expansions to try.
#'
#' @param neighbors \[`integer(1)`\]\cr
#' Number of trial points sampled in each neighborhood.
#'
#' @param beta \[`numeric(1)`\]\cr
#' Non-negative scaling factor for neighborhood expansion.
#'
#' @param iterlim \[`integer(1)`\]\cr
#' Maximum trust-region iterations during the main search.
#'
#' @param tolerance \[`numeric(1)`\]\cr
#' Minimum Euclidean distance for two optima to be treated as distinct.
#'
#' @param inferior_tolerance \[`numeric(1)`\]\cr
#' Maximum objective-value gap from the best known solution before a local
#' optimum is discarded early.
#'
#' @param gradient_tolerance \[`numeric(1)`\]\cr
#' Positive convergence tolerance for the scaled gradient.
#'
#' @param time_limit \[`numeric(1)` | `NULL`\]\cr
#' Optional time limit in seconds. If reached, the search stops early with a
#' warning.
#'
#' @param lower,upper \[`numeric(npar)` | `NULL`\]\cr
#' Optional lower and upper parameter bounds. Use \code{NULL} for unbounded
#' dimensions.
#'
#' @param collect_all \[`logical(1)`\]\cr
#' If \code{TRUE}, keep all converged local optima and disable early stopping
#' for optima that are inferior to the best known solution.
#'
#' @param quiet \[`logical(1)`\]\cr
#' If \code{TRUE}, suppress progress messages.
#'
#' @return
#' A \code{data.frame} summarizing the identified optima or \code{NULL} if none
#' could be determined.
#'
#' @export
#'
#' @examples
#' rosenbrock <- function(x) 100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
#' vntrs(f = rosenbrock, npar = 2)

vntrs <- function(
    f,
    npar,
    minimize = TRUE,
    init_runs = 5L,
    init_min = -1,
    init_max = 1,
    init_iterlim = 20L,
    neighborhoods = 5L,
    neighbors = 5L,
    beta = 0.05,
    iterlim = 100L,
    tolerance = 1e-6,
    inferior_tolerance = 1e-6,
    time_limit = NULL,
    lower = NULL,
    upper = NULL,
    collect_all = FALSE,
    quiet = TRUE,
    gradient_tolerance = 1e-6
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
  oeli::input_check_response(
    check = checkmate::check_flag(minimize),
    var_name = "minimize"
  )
  minimize <- isTRUE(minimize)
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
    check = checkmate::check_number(tolerance, finite = TRUE, lower = 0),
    var_name = "tolerance"
  )
  oeli::input_check_response(
    check = checkmate::check_number(
      inferior_tolerance, finite = TRUE, lower = 0
    ),
    var_name = "inferior_tolerance"
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
    check = checkmate::check_flag(collect_all),
    var_name = "collect_all"
  )
  collect_all <- isTRUE(collect_all)
  oeli::input_check_response(
    check = checkmate::check_flag(quiet),
    var_name = "quiet"
  )
  quiet <- isTRUE(quiet)
  oeli::input_check_response(
    check = checkmate::check_number(
      gradient_tolerance, finite = TRUE, lower = .Machine$double.xmin
    ),
    var_name = "gradient_tolerance"
  )
  .Call(
    `_vntrs_vntrs_cpp`, f, npar, minimize, init_runs, init_min, init_max,
    init_iterlim, neighborhoods, neighbors, beta, iterlim, tolerance,
    inferior_tolerance, has_time_limit, time_limit, lower, upper, quiet,
    collect_all, gradient_tolerance
  )
}

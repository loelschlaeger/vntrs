#include <RcppArmadillo.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <cfloat>
#include <string>
#include <vector>
#include "objective_validation.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

struct Controls {
  int npar;
  int init_runs;
  double init_min;
  double init_max;
  int init_iterlim;
  int neighborhoods;
  int neighbors;
  double beta;
  int iterlim;
  double tolerance;
  double inferior_tolerance;
  bool has_time_limit;
  double time_limit;
  int cores;
  arma::vec par_lower;
  arma::vec par_upper;
  bool collect_all_optima;
};

static Controls create_controls(
    int npar,
    int init_runs,
    double init_min,
    double init_max,
    int init_iterlim,
    int neighborhoods,
    int neighbors,
    double beta,
    int iterlim,
    double tolerance,
    double inferior_tolerance,
    bool has_time_limit,
    double time_limit,
    int cores,
    NumericVector par_lower_vec,
    NumericVector par_upper_vec,
    bool collect_all_optima
) {

  Controls cfg;

  if (init_runs <= 0) {
    stop("'init_runs' must be positive.");
  }
  if (init_iterlim <= 0) {
    stop("'init_iterlim' must be positive.");
  }
  if (neighborhoods <= 0) {
    stop("'neighborhoods' must be positive.");
  }
  if (neighbors <= 0) {
    stop("'neighbors' must be positive.");
  }
  if (!R_finite(beta) || beta < 0.0) {
    stop("'beta' must be finite and greater than or equal to zero.");
  }
  if (iterlim <= 0) {
    stop("'iterlim' must be positive.");
  }
  if (!R_finite(tolerance) || tolerance < 0.0) {
    stop("'tolerance' must be finite and greater than or equal to zero.");
  }
  if (!R_finite(inferior_tolerance) || inferior_tolerance < 0.0) {
    stop("'inferior_tolerance' must be finite and greater than or equal to zero.");
  }
  if (init_max < init_min) {
    stop("'init_max' must be greater than or equal to 'init_min'.");
  }
  if (has_time_limit) {
    if (!R_finite(time_limit) || time_limit <= 0.0) {
      stop("'time_limit' must be finite and positive.");
    }
  } else {
    time_limit = 0.0;
  }
  if (cores <= 0) {
    stop("'cores' must be positive.");
  }

  if (par_lower_vec.size() != npar) {
    stop("'lower' must have length matching 'npar'.");
  }
  if (par_upper_vec.size() != npar) {
    stop("'upper' must have length matching 'npar'.");
  }

  for (int i = 0; i < npar; ++i) {
    double low = par_lower_vec[i];
    double upp = par_upper_vec[i];
    if (!R_finite(low)) {
      low = R_NegInf;
    }
    if (!R_finite(upp)) {
      upp = R_PosInf;
    }
    if (std::isfinite(low) && std::isfinite(upp) && low > upp) {
      stop("Each element of 'lower' must be less than or equal to the corresponding element of 'upper'.");
    }
    par_lower_vec[i] = low;
    par_upper_vec[i] = upp;
  }

  cfg.npar = npar;
  cfg.init_runs = init_runs;
  cfg.init_min = init_min;
  cfg.init_max = init_max;
  cfg.init_iterlim = init_iterlim;
  cfg.neighborhoods = neighborhoods;
  cfg.neighbors = neighbors;
  cfg.beta = beta;
  cfg.iterlim = iterlim;
  cfg.tolerance = tolerance;
  cfg.inferior_tolerance = inferior_tolerance;
  cfg.has_time_limit = has_time_limit;
  cfg.time_limit = time_limit;
  cfg.cores = cores;
  cfg.par_lower = as<arma::vec>(par_lower_vec);
  cfg.par_upper = as<arma::vec>(par_upper_vec);
  cfg.collect_all_optima = collect_all_optima;

  try {
    Environment parallel_env = Environment::namespace_env("parallel");
    Function detectCores = parallel_env["detectCores"];
    SEXP detected = detectCores();
    if (!Rf_isNull(detected)) {
      int max_cores = as<int>(detected);
      if (max_cores > 0 && cfg.cores > max_cores) {
        warning("'cores' reduced to available cores.");
        cfg.cores = max_cores;
      }
    }
  } catch (const std::exception& ex) {
    warning("Failed to query available cores: %s", ex.what());
  }

  return cfg;
}

static void check_function(Function f, int npar, const Controls& controls) {
  int test_runs = 10;
  NumericVector samples = runif(test_runs * npar, controls.init_min, controls.init_max);
  NumericMatrix points(test_runs, npar);
  for (int i = 0, idx = 0; i < test_runs; ++i) {
    for (int j = 0; j < npar; ++j, ++idx) {
      double value = std::round(samples[idx] * 10.0) / 10.0;
      points(i, j) = value;
    }
  }

  for (int run = 0; run < test_runs; ++run) {
    NumericVector point = points(run, _);
    arma::vec point_vec = Rcpp::as<arma::vec>(point);
    vntrs::parse_objective(
      f,
      point_vec,
      /*require_finite_gradient=*/true,
      /*require_finite_hessian=*/true);
  }
}

struct OptimaStorage {
  int npar;
  std::vector<arma::vec> arguments;
  std::vector<double> values;

  explicit OptimaStorage(int npar_) : npar(npar_) {}

  bool empty() const { return values.empty(); }
  std::size_t size() const { return values.size(); }

  void append(const arma::vec& argument, double value) {
    arguments.push_back(argument);
    values.push_back(value);
  }

  bool unique(const arma::vec& argument, double tolerance) const {
    if (arguments.empty()) {
      return true;
    }
    if (!argument.is_finite()) {
      return false;
    }
    double tol_sq = tolerance * tolerance;
    for (std::size_t i = 0; i < arguments.size(); ++i) {
      const arma::vec& current = arguments[i];
      arma::vec diff = current - argument;
      double dist_sq = arma::dot(diff, diff);
      if (!std::isfinite(dist_sq)) {
        dist_sq = std::numeric_limits<double>::infinity();
      }
      if (dist_sq < tol_sq) {
        return false;
      }
    }
    return true;
  }

  int best_index(bool minimize) const {
    if (values.empty()) {
      return -1;
    }
    int idx = -1;
    double best = minimize ? std::numeric_limits<double>::infinity() : -std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < values.size(); ++i) {
      double val = values[i];
      if (!std::isfinite(val)) {
        continue;
      }
      if (minimize) {
        if (val < best) {
          best = val;
          idx = static_cast<int>(i);
        }
      } else {
        if (val > best) {
          best = val;
          idx = static_cast<int>(i);
        }
      }
    }
    return idx;
  }

  double best_value(bool minimize) const {
    int idx = best_index(minimize);
    if (idx < 0) {
      return NA_REAL;
    }
    return values[idx];
  }

  arma::vec best_argument(bool minimize) const {
    int idx = best_index(minimize);
    if (idx < 0) {
      return arma::vec(npar, arma::fill::value(NA_REAL));
    }
    return arguments[idx];
  }
};

static bool should_interrupt(
   Function f,
   const arma::vec& point,
   const OptimaStorage& storage,
   bool minimize,
   double inferior_tolerance,
   bool quiet,
   bool collect_all_optima
) {
  if (collect_all_optima) {
    return false;
  }
  if (storage.empty()) {
    return false;
  }
  if (!point.is_finite()) {
    return false;
  }

  arma::vec best_argument = storage.best_argument(minimize);
  double best_value = storage.best_value(minimize);
  if (!std::isfinite(best_value)) {
    return false;
  }

  for (std::size_t i = 0; i < storage.arguments.size(); ++i) {
    arma::vec diff = storage.arguments[i] - point;
    double dist_sq = arma::dot(diff, diff);
    if (!std::isfinite(dist_sq)) {
      dist_sq = std::numeric_limits<double>::infinity();
    }
    if (dist_sq <= 1.0) {
      if (!quiet) {
        Rcpp::Rcout << " [optimum already visited]";
      }
      return true;
    }
  }

  vntrs::ObjectiveComponents eval = vntrs::parse_objective(
    f,
    point,
    /*require_finite_gradient=*/false,
    /*require_finite_hessian=*/false);
  arma::vec gradient = eval.gradient;
  double value = eval.value;
  if (!gradient.is_finite() || !std::isfinite(value)) {
    return false;
  }

  double grad_norm_sq = arma::dot(gradient, gradient);
  if (grad_norm_sq <= std::pow(1e-3, 2)) {
    if (minimize) {
      if (value > best_value + inferior_tolerance) {
        if (!quiet) {
          Rcpp::Rcout << " [optimum inferior to best known]";
        }
        return true;
      }
    } else {
      if (value < best_value - inferior_tolerance) {
        if (!quiet) {
          Rcpp::Rcout << " [optimum inferior to best known]";
        }
        return true;
      }
    }
  }

  return false;
}

// Declaration from trust_region.cpp
List trust_region_cpp(
    Function objfun,
    NumericVector parinit,
    double rinit,
    double rmax,
    int iterlim,
    bool minimize,
    double tol,
    double eta,
    NumericVector lower,
    NumericVector upper
);

static List run_local(
    Function f,
    const arma::vec& parinit,
    bool minimize,
    const Controls& controls,
    OptimaStorage& storage,
    bool quiet
) {
  int batches = storage.empty() ? 1 : controls.iterlim;
  arma::vec current = parinit;
  List last;

  for (int b = 0; b < batches; ++b) {
    int iterlim = std::max(1, controls.iterlim / (storage.empty() ? 1 : batches));
    last = trust_region_cpp(
      f,
      wrap(current),
      1.0,
      10.0,
      iterlim,
      minimize,
      1e-6,
      0.1,
      wrap(controls.par_lower),
      wrap(controls.par_upper)
    );

    bool converged = as<bool>(last["converged"]);
    arma::vec argument = as<arma::vec>(last["argument"]);

    if (b < batches - 1) {
      if (converged) {
        break;
      }
      if (should_interrupt(
            f, argument, storage, minimize, controls.inferior_tolerance,
            quiet, controls.collect_all_optima)
          )
        {
        NumericVector arg_out(storage.npar, NA_REAL);
        return List::create(
          Named("success") = false,
          Named("value") = NA_REAL,
          Named("argument") = arg_out
        );
      }
      current = argument;
    } else {
      current = argument;
    }
  }

  arma::vec argument = as<arma::vec>(last["argument"]);
  NumericVector argument_out(storage.npar, NA_REAL);
  if (argument.n_elem == static_cast<arma::uword>(storage.npar) && argument.is_finite()) {
    argument_out = wrap(argument);
  }
  double value = as<double>(last["value"]);
  if (!std::isfinite(value)) {
    value = NA_REAL;
  }

  return List::create(
    Named("success") = as<bool>(last["converged"]),
    Named("value") = value,
    Named("argument") = argument_out
  );
}

static std::vector<arma::vec> select_neighbors(
    Function f, const arma::vec& x, double expansion, const Controls& controls
) {
  vntrs::ObjectiveComponents eval = vntrs::parse_objective(
    f,
    x,
    /*require_finite_gradient=*/false,
    /*require_finite_hessian=*/false);
  arma::mat hessian = eval.hessian;
  if (!hessian.is_finite()) {
    hessian.eye(controls.npar, controls.npar);
  }
  arma::mat sym = 0.5 * (hessian + hessian.t());
  if (!sym.is_finite()) {
    sym.eye(controls.npar, controls.npar);
  }
  arma::vec eigval;
  arma::mat eigvec;
  bool eigen_ok = arma::eig_sym(eigval, eigvec, sym);
  if (!eigen_ok || !eigval.is_finite() || !eigvec.is_finite()) {
    eigval.set_size(controls.npar);
    eigval.fill(1.0);
    eigvec.eye(controls.npar, controls.npar);
  }

  arma::vec scaled = controls.beta * eigval / expansion;
  if (!scaled.is_finite()) {
    scaled.zeros();
  } else {
    double max_val = scaled.max();
    scaled -= max_val;
  }
  arma::vec weights = arma::exp(scaled);
  double sum_weights = arma::accu(weights);
  if (!std::isfinite(sum_weights) || sum_weights <= 0.0) {
    weights.fill(1.0 / eigval.n_elem);
  } else {
    weights /= sum_weights;
  }

  std::vector<double> cumulative(weights.n_elem, 0.0);
  double cumulative_sum = 0.0;
  for (arma::uword i = 0; i < weights.n_elem; ++i) {
    cumulative_sum += weights(i);
    cumulative[i] = cumulative_sum;
  }
  cumulative.back() = 1.0;

  std::vector<arma::vec> neighbors;
  neighbors.reserve(controls.neighbors);
  for (int n = 0; n < controls.neighbors; ++n) {
    double alpha = R::runif(0.75, 1.0);
    double direction = (R::runif(0.0, 1.0) < 0.5) ? -1.0 : 1.0;
    double draw = R::runif(0.0, 1.0);
    arma::uword index = 0;
    while (index < cumulative.size() && draw > cumulative[index]) {
      ++index;
    }
    if (index >= eigvec.n_cols) {
      index = eigvec.n_cols - 1;
    }
    arma::vec neighbor = x + expansion * alpha * direction * eigvec.col(index);
    for (arma::uword j = 0; j < neighbor.n_elem; ++j) {
      double low = controls.par_lower(j);
      double upp = controls.par_upper(j);
      if (std::isfinite(low) && neighbor(j) < low) {
        neighbor(j) = low;
      }
      if (std::isfinite(upp) && neighbor(j) > upp) {
        neighbor(j) = upp;
      }
    }
    neighbors.push_back(neighbor);
  }
  return neighbors;
}

static arma::vec generate_start(int npar, const Controls& controls) {
  arma::vec point(npar);
  for (int i = 0; i < npar; ++i) {
    double value = R::runif(controls.init_min, controls.init_max);
    if (std::isfinite(controls.par_lower(i)) && value < controls.par_lower(i)) {
      value = controls.par_lower(i);
    }
    if (std::isfinite(controls.par_upper(i)) && value > controls.par_upper(i)) {
      value = controls.par_upper(i);
    }
    point(i) = value;
  }
  return point;
}

static List initialize_search(
    Function f, int npar, bool minimize, const Controls& controls,
    OptimaStorage& storage, bool quiet
) {
  std::vector<List> results;
  results.reserve(controls.init_runs);
  if (!quiet) {
    Rcpp::Rcout << "* Apply local search at " << controls.init_runs <<
      " random starting points.\n";
  }

  for (int run = 0; run < controls.init_runs; ++run) {
    arma::vec start = generate_start(npar, controls);
    auto run_start = std::chrono::steady_clock::now();
    List local = trust_region_cpp(
      f,
      wrap(start),
      1.0,
      10.0,
      std::max(1, controls.init_iterlim),
      minimize,
      1e-6,
      0.1,
      wrap(controls.par_lower),
      wrap(controls.par_upper)
    );
    auto run_end = std::chrono::steady_clock::now();
    double duration = std::chrono::duration<double>(run_end - run_start).count();

    bool success = as<bool>(local["converged"]);
    arma::vec argument = as<arma::vec>(local["argument"]);
    double value = as<double>(local["value"]);

    if (!quiet) {
      Rcpp::Rcout << "** Run " << (run + 1);
      Rcpp::Rcout << " [" << std::round(duration) << " s]";
    }

    if (success && argument.is_finite() && std::isfinite(value)) {
      if (!quiet) {
        Rcpp::Rcout << " [found optimum]";
      }
      if (storage.unique(argument, controls.tolerance)) {
        storage.append(argument, value);
        if (!quiet) {
          Rcpp::Rcout << " [optimum is unknown]";
        }
      }
    }

    if (!quiet) {
      Rcpp::Rcout << "\n";
    }

    results.push_back(local);
  }

  if (results.empty()) {
    NumericVector missing_x(npar, NA_REAL);
    return List::create(
      Named("success") = false,
      Named("x_best") = missing_x
    );
  }

  std::vector<double> candidate_values(results.size());
  for (std::size_t i = 0; i < results.size(); ++i) {
    double value = as<double>(results[i]["value"]);
    if (!std::isfinite(value)) {
      candidate_values[i] = minimize ? std::numeric_limits<double>::infinity() : -std::numeric_limits<double>::infinity();
    } else {
      candidate_values[i] = value;
    }
  }

  int best_index = 0;
  if (minimize) {
    best_index = std::distance(candidate_values.begin(), std::min_element(candidate_values.begin(), candidate_values.end()));
  } else {
    best_index = std::distance(candidate_values.begin(), std::max_element(candidate_values.begin(), candidate_values.end()));
  }

  List best_run = results[best_index];
  bool best_success = as<bool>(best_run["converged"]);
  arma::vec best_argument = as<arma::vec>(best_run["argument"]);
  double best_value = as<double>(best_run["value"]);
  if (best_success && best_argument.is_finite() && std::isfinite(best_value)) {
    if (storage.unique(best_argument, controls.tolerance)) {
      storage.append(best_argument, best_value);
    }
    return List::create(
      Named("success") = true,
      Named("x_best") = wrap(best_argument)
    );
  }

  if (!quiet) {
    Rcpp::Rcout << "* Continue the best run " << (best_index + 1) << ".\n";
  }

  arma::vec restart = best_argument;
  if (!restart.is_finite()) {
    restart = generate_start(npar, controls);
  }
  List extended = trust_region_cpp(
    f,
    wrap(restart),
    1.0,
    10.0,
    controls.iterlim,
    minimize,
    1e-6,
    0.1,
    wrap(controls.par_lower),
    wrap(controls.par_upper)
  );
  bool extended_success = as<bool>(extended["converged"]);
  arma::vec ext_argument = as<arma::vec>(extended["argument"]);
  double ext_value = as<double>(extended["value"]);
  if (extended_success && ext_argument.is_finite() && std::isfinite(ext_value)) {
    if (!quiet) {
      Rcpp::Rcout << " [found optimum]\n";
    }
    if (storage.unique(ext_argument, controls.tolerance)) {
      storage.append(ext_argument, ext_value);
    }
    return List::create(
      Named("success") = true,
      Named("x_best") = wrap(ext_argument)
    );
  }

  if (!quiet) {
    Rcpp::Rcout << " [failed]\n";
  }

  NumericVector missing_x(npar, NA_REAL);
  return List::create(
    Named("success") = false,
    Named("x_best") = missing_x
  );
}

// [[Rcpp::export]]
SEXP vntrs_cpp(
    Function f,
    int npar,
    bool minimize,
    int init_runs,
    double init_min,
    double init_max,
    int init_iterlim,
    int neighborhoods,
    int neighbors,
    double beta,
    int iterlim,
    double tolerance,
    double inferior_tolerance,
    bool has_time_limit,
    double time_limit,
    int cores,
    NumericVector lower,
    NumericVector upper,
    bool quiet,
    bool collect_all_optima
) {

  if (npar <= 0) {
    stop("'npar' must be positive.");
  }

  Controls controls = create_controls(
    npar,
    init_runs,
    init_min,
    init_max,
    init_iterlim,
    neighborhoods,
    neighbors,
    beta,
    iterlim,
    tolerance,
    inferior_tolerance,
    has_time_limit,
    time_limit,
    cores,
    lower,
    upper,
    collect_all_optima
  );
  check_function(f, npar, controls);

  OptimaStorage storage(npar);

  auto start_time = std::chrono::steady_clock::now();

  if (!quiet) {
    Rcpp::Rcout << "Initialize VNTRS.\n";
  }
  List init = initialize_search(f, npar, minimize, controls, storage, quiet);
  bool init_success = as<bool>(init["success"]);
  arma::vec x_best = as<arma::vec>(init["x_best"]);
  if (!init_success || !x_best.is_finite()) {
    warning("No optima found.");
    return R_NilValue;
    if (storage.empty()) {
      warning("No optima found.");
      return R_NilValue;
    }
    arma::vec fallback = storage.best_argument(minimize);
    if (!fallback.is_finite()) {
      warning("No optima found.");
      return R_NilValue;
    }
    x_best = fallback;
  }

  if (!quiet) {
    Rcpp::Rcout << "Start VNTRS.\n";
  }

  int k = 1;
  bool stop_loop = false;
  while (k <= controls.neighborhoods) {
    if (stop_loop) {
      break;
    }
    if (!quiet) {
      Rcpp::Rcout << "* Select neighborhood " << k << ".\n";
    }
    double expansion = std::pow(1.5, k - 1);
    std::vector<arma::vec> neighbors = select_neighbors(f, x_best, expansion, controls);

    for (std::size_t j = 0; j < neighbors.size(); ++j) {
      if (controls.has_time_limit) {
        auto now = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(now - start_time).count();
        if (elapsed > controls.time_limit) {
          warning("Stopped early because 'time_limit' was reached.");
          stop_loop = true;
          break;
        }
      }

      if (!quiet) {
        Rcpp::Rcout << "** Neighbor " << (j + 1);
      }

      auto neighbor_start = std::chrono::steady_clock::now();
      List local = run_local(f, neighbors[j], minimize, controls, storage, quiet);
      auto neighbor_end = std::chrono::steady_clock::now();
      double duration = std::chrono::duration<double>(neighbor_end - neighbor_start).count();
      if (!quiet) {
        Rcpp::Rcout << " [" << std::round(duration) << " s]";
      }

      bool success = as<bool>(local["success"]);
      NumericVector argument = local["argument"];
      double value = as<double>(local["value"]);
      bool valid_argument = argument.size() == npar;
      bool finite_argument = true;
      for (int idx = 0; idx < argument.size(); ++idx) {
        if (!R_finite(argument[idx])) {
          finite_argument = false;
          break;
        }
      }
      if (success && valid_argument && finite_argument && std::isfinite(value)) {
        if (!quiet) {
          Rcpp::Rcout << " [found optimum]";
        }
        arma::vec arg_vec = as<arma::vec>(argument);
        if (storage.unique(arg_vec, controls.tolerance)) {
          if (!quiet) {
            Rcpp::Rcout << " [optimum is unknown]";
          }
          storage.append(arg_vec, value);
        }
      }
      if (!quiet) {
        Rcpp::Rcout << "\n";
      }
    }

    if (storage.empty()) {
      break;
    }
    arma::vec x_new = storage.best_argument(minimize);
    arma::vec diff = x_new - x_best;
    double dist_sq = arma::dot(diff, diff);
    if (!std::isfinite(dist_sq) || dist_sq > controls.tolerance * controls.tolerance) {
      if (!quiet) {
        Rcpp::Rcout << "* Reset neighborhood, because better optimum was found.\n";
      }
      x_best = x_new;
      k = 1;
    } else {
      ++k;
    }
  }

  if (storage.empty()) {
    warning("No optima found.");
    return R_NilValue;
  }

  double best_value = storage.best_value(minimize);
  if (!std::isfinite(best_value)) {
    warning("No finite optima found.");
    return R_NilValue;
  }

  std::size_t n_opt = storage.size();
  NumericMatrix args(n_opt, npar);
  NumericVector values(n_opt);
  for (std::size_t i = 0; i < n_opt; ++i) {
    arma::vec arg = storage.arguments[i];
    for (int j = 0; j < npar; ++j) {
      args(i, j) = arg(j);
    }
    values[i] = storage.values[i];
  }

  double tol_value = std::sqrt(DBL_EPSILON) * std::max(1.0, std::fabs(best_value));
  LogicalVector global(n_opt);
  for (std::size_t i = 0; i < n_opt; ++i) {
    global[i] = std::fabs(values[i] - best_value) <= tol_value;
  }

  List df(npar + 2);
  CharacterVector names(npar + 2);
  for (int j = 0; j < npar; ++j) {
    df[j] = args(_, j);
    names[j] = std::string("p") + std::to_string(j + 1);
  }
  df[npar] = values;
  names[npar] = "value";
  df[npar + 1] = global;
  names[npar + 1] = "global";

  df.attr("names") = names;
  df.attr("class") = "data.frame";
  df.attr("row.names") = IntegerVector::create(NA_INTEGER, static_cast<int>(n_opt));

  if (!quiet) {
    Rcpp::Rcout << "Done.\n";
  }

  return df;
}

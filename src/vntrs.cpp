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
  double interruption_gradient_tolerance;
  double gradient_tolerance;
  bool has_time_limit;
  double time_limit;
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
    double interruption_gradient_tolerance,
    bool has_time_limit,
    double time_limit,
    NumericVector par_lower_vec,
    NumericVector par_upper_vec,
    bool collect_all_optima,
    double gradient_tolerance
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
  if (!R_finite(interruption_gradient_tolerance) ||
      interruption_gradient_tolerance < 0.0) {
    stop("'interruption_gradient_tolerance' must be finite and greater than or equal to zero.");
  }
  if (!R_finite(gradient_tolerance) || gradient_tolerance <= 0.0) {
    stop("'gradient_tolerance' must be finite and positive.");
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
  cfg.interruption_gradient_tolerance = interruption_gradient_tolerance;
  cfg.gradient_tolerance = gradient_tolerance;
  cfg.has_time_limit = has_time_limit;
  cfg.time_limit = time_limit;
  cfg.par_lower = as<arma::vec>(par_lower_vec);
  cfg.par_upper = as<arma::vec>(par_upper_vec);
  cfg.collect_all_optima = collect_all_optima;

  return cfg;
}

static arma::vec initialization_center(int npar, const Controls& controls) {
  arma::vec point(npar);
  for (int i = 0; i < npar; ++i) {
    double low = controls.par_lower(i);
    double upp = controls.par_upper(i);
    if (std::isfinite(low) && std::isfinite(upp)) {
      point(i) = 0.5 * low + 0.5 * upp;
    } else {
      point(i) = 0.5 * controls.init_min + 0.5 * controls.init_max;
      if (std::isfinite(low)) {
        point(i) = std::max(point(i), low);
      }
      if (std::isfinite(upp)) {
        point(i) = std::min(point(i), upp);
      }
    }
  }
  return point;
}

static void check_function(Function f, int npar, const Controls& controls) {
  arma::vec point = initialization_center(npar, controls);
  vntrs::parse_objective(
    f, point, /*need_gradient=*/false, controls.par_lower, controls.par_upper
  );
}

static bool same_point(const arma::vec& x, const arma::vec& y,
                       double tolerance) {
  arma::vec scale = arma::max(
    arma::max(arma::abs(x), arma::abs(y)), arma::ones<arma::vec>(x.n_elem)
  );
  return arma::max(arma::abs(x - y) / scale) <= tolerance;
}

struct OptimaStorage {
  int npar;
  std::vector<arma::vec> arguments;
  std::vector<double> values;
  std::vector<arma::mat> models;

  explicit OptimaStorage(int npar_) : npar(npar_) {}

  bool empty() const { return values.empty(); }
  std::size_t size() const { return values.size(); }

  void append(const arma::vec& argument,
              double value,
              const arma::mat& model) {
    arguments.push_back(argument);
    values.push_back(value);
    if (model.n_rows == static_cast<arma::uword>(npar) &&
        model.n_cols == static_cast<arma::uword>(npar) && model.is_finite()) {
      models.push_back(model);
    } else {
      models.push_back(arma::eye<arma::mat>(npar, npar));
    }
  }

  bool unique(const arma::vec& argument, double tolerance) const {
    if (arguments.empty()) {
      return true;
    }
    if (!argument.is_finite()) {
      return false;
    }
    for (std::size_t i = 0; i < arguments.size(); ++i) {
      if (same_point(arguments[i], argument, tolerance)) {
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

  arma::mat best_model(bool minimize) const {
    int idx = best_index(minimize);
    if (idx < 0 || static_cast<std::size_t>(idx) >= models.size()) {
      return arma::eye<arma::mat>(npar, npar);
    }
    return models[idx];
  }
};

static bool should_interrupt(
   const arma::vec& point,
   double value,
   const arma::vec& gradient,
   const OptimaStorage& storage,
   bool minimize,
   double inferior_tolerance,
   double interruption_gradient_tolerance,
   double optimum_tolerance,
   bool quiet,
   bool collect_all_optima
) {
  if (storage.empty()) {
    return false;
  }
  if (!point.is_finite()) {
    return false;
  }

  double best_value = storage.best_value(minimize);
  if (!std::isfinite(best_value)) {
    return false;
  }

  bool near_known_optimum = false;
  for (std::size_t i = 0; i < storage.arguments.size(); ++i) {
    arma::vec diff = storage.arguments[i] - point;
    double dist_sq = arma::dot(diff, diff);
    if (!std::isfinite(dist_sq)) {
      dist_sq = std::numeric_limits<double>::infinity();
    }
    near_known_optimum = near_known_optimum || dist_sq <= 1.0;
    if (same_point(storage.arguments[i], point, optimum_tolerance)) {
      if (!quiet) {
        Rcpp::Rcout << " [optimum already visited]";
      }
      return true;
    }
  }
  if (collect_all_optima) return false;

  if (!gradient.is_finite() || !std::isfinite(value)) {
    return false;
  }

  double grad_norm_sq = arma::dot(gradient, gradient);
  if (grad_norm_sq <= std::pow(interruption_gradient_tolerance, 2)) {
    bool no_meaningful_improvement = minimize
      ? value >= best_value - inferior_tolerance
      : value <= best_value + inferior_tolerance;
    if (near_known_optimum && no_meaningful_improvement) {
      if (!quiet) {
        Rcpp::Rcout << " [optimum already visited]";
      }
      return true;
    }
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
    NumericVector upper,
    Nullable<NumericMatrix> initial_model = R_NilValue,
    double initial_delta = -1.0
);

static List run_local(
    Function f,
    const arma::vec& parinit,
    bool minimize,
    const Controls& controls,
    OptimaStorage& storage,
    bool quiet,
    int iteration_limit
) {
  vntrs::ObjectiveComponents start = vntrs::parse_objective(
    f, parinit, /*need_gradient=*/false, controls.par_lower,
    controls.par_upper, /*allow_nonfinite=*/true
  );
  if (!std::isfinite(start.value)) {
    return List::create(
      Named("success") = false,
      Named("value") = NA_REAL,
      Named("argument") = NumericVector(storage.npar, NA_REAL)
    );
  }
  arma::vec current = parinit;
  int remaining = std::max(1, iteration_limit);
  List last;
  bool has_state = false;
  bool warm_start = !storage.empty();
  while (remaining > 0) {
    int chunk = storage.empty() ? remaining : std::min(20, remaining);
    if (!has_state) {
      if (warm_start) {
        NumericMatrix warm_model = wrap(storage.best_model(minimize));
        last = trust_region_cpp(
          f,
          wrap(current),
          1.0,
          10.0,
          chunk,
          minimize,
          controls.gradient_tolerance,
          0.1,
          wrap(controls.par_lower),
          wrap(controls.par_upper),
          warm_model,
          1.0
        );
      } else {
        last = trust_region_cpp(
          f,
          wrap(current),
          1.0,
          10.0,
          chunk,
          minimize,
          controls.gradient_tolerance,
          0.1,
          wrap(controls.par_lower),
          wrap(controls.par_upper)
        );
      }
      has_state = true;
    } else {
      NumericMatrix previous_model = last["model"];
      double previous_radius = as<double>(last["radius"]);
      last = trust_region_cpp(
        f,
        wrap(current),
        1.0,
        10.0,
        chunk,
        minimize,
        controls.gradient_tolerance,
        0.1,
        wrap(controls.par_lower),
        wrap(controls.par_upper),
        previous_model,
        previous_radius
      );
    }
    remaining -= chunk;
    bool converged = as<bool>(last["converged"]);
    current = as<arma::vec>(last["argument"]);
    if (converged || remaining == 0) {
      break;
    }
    if (should_interrupt(
          current, as<double>(last["value"]),
          as<arma::vec>(last["gradient"]), storage, minimize,
          controls.inferior_tolerance,
          controls.interruption_gradient_tolerance, controls.tolerance,
          quiet, controls.collect_all_optima)) {
      NumericVector arg_out(storage.npar, NA_REAL);
      return List::create(
        Named("success") = false,
        Named("value") = NA_REAL,
        Named("argument") = arg_out
      );
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
    Named("argument") = argument_out,
    Named("model") = last["model"],
    Named("radius") = last["radius"]
  );
}

static std::vector<arma::vec> select_neighbors(
    const arma::vec& x,
    double expansion,
    const Controls& controls,
    const arma::mat& curvature
) {
  arma::mat hessian = curvature;
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
    double low = controls.par_lower(i);
    double upp = controls.par_upper(i);
    double value;
    if (std::isfinite(low) && std::isfinite(upp)) {
      value = R::runif(low, upp);
    } else {
      value = R::runif(controls.init_min, controls.init_max);
      if (std::isfinite(low) && value < low) {
        value = low;
      }
      if (std::isfinite(upp) && value > upp) {
        value = upp;
      }
    }
    point(i) = value;
  }
  return point;
}

static arma::vec repair_start(Function f, arma::vec point,
                              const Controls& controls) {
  arma::vec center = initialization_center(controls.npar, controls);
  for (int attempt = 0; attempt < 32; ++attempt) {
    vntrs::ObjectiveComponents value = vntrs::parse_objective(
      f, point, /*need_gradient=*/false, controls.par_lower,
      controls.par_upper, /*allow_nonfinite=*/true
    );
    if (std::isfinite(value.value)) return point;
    point = 0.5 * (point + center);
  }
  return center;
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
    arma::vec start = repair_start(
      f, generate_start(npar, controls), controls
    );
    if (!quiet) {
      Rcpp::Rcout << "** Run " << (run + 1);
    }
    auto run_start = std::chrono::steady_clock::now();
    List local = run_local(
      f,
      start,
      minimize,
      controls,
      storage,
      quiet,
      controls.init_iterlim
    );
    auto run_end = std::chrono::steady_clock::now();
    double duration = std::chrono::duration<double>(run_end - run_start).count();

    bool success = as<bool>(local["success"]);
    arma::vec argument = as<arma::vec>(local["argument"]);
    double value = as<double>(local["value"]);

    if (!quiet) {
      Rcpp::Rcout << " [" << std::round(duration) << " s]";
    }

    if (success && argument.is_finite() && std::isfinite(value)) {
      if (!quiet) {
        Rcpp::Rcout << " [found optimum]";
      }
      if (storage.unique(argument, controls.tolerance)) {
        arma::mat local_model = as<arma::mat>(local["model"]);
        storage.append(argument, value, local_model);
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
  bool best_success = as<bool>(best_run["success"]);
  arma::vec best_argument = as<arma::vec>(best_run["argument"]);
  double best_value = as<double>(best_run["value"]);
  if (best_success && best_argument.is_finite() && std::isfinite(best_value)) {
    if (storage.unique(best_argument, controls.tolerance)) {
      arma::mat best_model = as<arma::mat>(best_run["model"]);
      storage.append(best_argument, best_value, best_model);
    }
    return List::create(
      Named("success") = true,
      Named("x_best") = wrap(best_argument)
    );
  }

  if (!quiet) {
    Rcpp::Rcout << "* Continue the best run " << (best_index + 1) << ".";
  }

  arma::vec restart = best_argument;
  if (!restart.is_finite()) {
    restart = repair_start(f, generate_start(npar, controls), controls);
  }
  List extended = run_local(
    f, restart, minimize, controls, storage, quiet, controls.iterlim
  );
  bool extended_success = as<bool>(extended["success"]);
  arma::vec ext_argument = as<arma::vec>(extended["argument"]);
  double ext_value = as<double>(extended["value"]);
  if (extended_success && ext_argument.is_finite() && std::isfinite(ext_value)) {
    if (!quiet) {
      Rcpp::Rcout << " [found optimum]\n";
    }
    if (storage.unique(ext_argument, controls.tolerance)) {
      arma::mat ext_model = as<arma::mat>(extended["model"]);
      storage.append(ext_argument, ext_value, ext_model);
    }
    return List::create(
      Named("success") = true,
      Named("x_best") = wrap(ext_argument)
    );
  }

  if (!quiet) {
    Rcpp::Rcout << " [failed]\n";
  }

  if (!storage.empty()) {
    return List::create(
      Named("success") = true,
      Named("x_best") = wrap(storage.best_argument(minimize))
    );
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
    double interruption_gradient_tolerance,
    bool has_time_limit,
    double time_limit,
    NumericVector lower,
    NumericVector upper,
    bool quiet,
    bool collect_all_optima,
    double gradient_tolerance
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
    interruption_gradient_tolerance,
    has_time_limit,
    time_limit,
    lower,
    upper,
    collect_all_optima,
    gradient_tolerance
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
  }
  double x_best_value = storage.best_value(minimize);

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
    std::vector<arma::vec> neighbors = select_neighbors(
      x_best, expansion, controls, storage.best_model(minimize)
    );

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
      List local = run_local(
        f, neighbors[j], minimize, controls, storage, quiet, controls.iterlim
      );
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
          arma::mat local_model = as<arma::mat>(local["model"]);
          storage.append(arg_vec, value, local_model);
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
    double x_new_value = storage.best_value(minimize);
    arma::vec diff = x_new - x_best;
    double dist_sq = arma::dot(diff, diff);
    double value_tolerance = std::sqrt(DBL_EPSILON) * std::max(
      1.0, std::max(std::fabs(x_best_value), std::fabs(x_new_value))
    );
    double improvement = minimize
      ? x_best_value - x_new_value : x_new_value - x_best_value;
    bool meaningfully_better = std::isfinite(improvement) &&
      improvement > value_tolerance;
    bool equivalent_new_optimum = std::isfinite(dist_sq) &&
      std::fabs(x_new_value - x_best_value) <= value_tolerance &&
      dist_sq > 1.0;
    if (meaningfully_better || equivalent_new_optimum) {
      if (!quiet) {
        Rcpp::Rcout << "* Reset neighborhood, because better optimum was found.\n";
      }
      x_best = x_new;
      x_best_value = x_new_value;
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
  std::vector<std::size_t> order(n_opt);
  for (std::size_t i = 0; i < n_opt; ++i) {
    order[i] = i;
  }
  std::stable_sort(
    order.begin(),
    order.end(),
    [&storage, minimize](std::size_t left, std::size_t right) {
      return minimize
        ? storage.values[left] < storage.values[right]
        : storage.values[left] > storage.values[right];
    }
  );

  NumericMatrix args(n_opt, npar);
  NumericVector values(n_opt);
  for (std::size_t i = 0; i < n_opt; ++i) {
    std::size_t source = order[i];
    arma::vec arg = storage.arguments[source];
    for (int j = 0; j < npar; ++j) {
      args(i, j) = arg(j);
    }
    values[i] = storage.values[source];
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

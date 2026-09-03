#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include "objective_validation.h"
#include "optimization_metrics.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

static void apply_bounds(arma::vec& x, const arma::vec& lower,
                         const arma::vec& upper) {
  for (arma::uword i = 0; i < x.n_elem; ++i) {
    if (std::isfinite(lower(i))) x(i) = std::max(x(i), lower(i));
    if (std::isfinite(upper(i))) x(i) = std::min(x(i), upper(i));
  }
}

static arma::mat positive_definite_model(const arma::mat& candidate) {
  arma::uword n = candidate.n_rows;
  arma::mat model = 0.5 * (candidate + candidate.t());
  if (!model.is_finite()) {
    return arma::eye<arma::mat>(n, n);
  }

  arma::mat chol_factor;
  if (arma::chol(chol_factor, model)) {
    return model;
  }

  arma::vec eigenvalues;
  bool eigen_ok = arma::eig_sym(eigenvalues, model);
  if (!eigen_ok || !eigenvalues.is_finite()) {
    return arma::eye<arma::mat>(n, n);
  }
  double scale = std::max(1.0, arma::norm(model, "fro"));
  double floor = std::sqrt(std::numeric_limits<double>::epsilon()) * scale;
  double shift = std::max(0.0, floor - eigenvalues.min());
  model.diag() += shift;
  if (!arma::chol(chol_factor, model)) {
    return scale * arma::eye<arma::mat>(n, n);
  }
  return model;
}

static void project_onto_critical_cone(arma::vec& direction,
                                       const arma::ivec& cone_sign) {
  for (arma::uword i = 0; i < direction.n_elem; ++i) {
    if (cone_sign(i) > 0 && direction(i) < 0.0) {
      direction(i) = 0.0;
    } else if (cone_sign(i) < 0 && direction(i) > 0.0) {
      direction(i) = 0.0;
    }
  }
}

static bool negative_curvature_from_seed(const arma::mat& curvature,
                                         const arma::ivec& cone_sign,
                                         arma::vec direction,
                                         double curvature_tolerance,
                                         double spectral_scale) {
  project_onto_critical_cone(direction, cone_sign);
  double norm = arma::norm(direction, 2);
  if (!std::isfinite(norm) || norm <= std::numeric_limits<double>::epsilon()) {
    return false;
  }
  direction /= norm;
  double step_size = 0.25 / std::max(1.0, spectral_scale);
  for (int iteration = 0; iteration < 40; ++iteration) {
    arma::vec product = curvature * direction;
    double rayleigh = arma::dot(direction, product);
    if (rayleigh < -curvature_tolerance) {
      return true;
    }
    direction -= step_size * (product - rayleigh * direction);
    project_onto_critical_cone(direction, cone_sign);
    norm = arma::norm(direction, 2);
    if (!std::isfinite(norm) || norm <= std::numeric_limits<double>::epsilon()) {
      return false;
    }
    direction /= norm;
  }
  return false;
}

static bool second_order_satisfied(const arma::mat& curvature,
                                   const arma::vec& gradient,
                                   const arma::vec& x,
                                   const arma::vec& lower,
                                   const arma::vec& upper) {
  if (!curvature.is_finite() || curvature.n_rows != x.n_elem ||
      curvature.n_cols != x.n_elem) {
    return false;
  }
  arma::uvec cone_indices(x.n_elem);
  arma::ivec cone_sign(x.n_elem);
  arma::uword cone_dimension = 0;
  double gradient_tolerance = std::sqrt(
    std::numeric_limits<double>::epsilon()
  ) * std::max(1.0, arma::norm(gradient, "inf"));
  for (arma::uword i = 0; i < x.n_elem; ++i) {
    double lower_tol = std::isfinite(lower(i))
      ? 1e-8 * (1.0 + std::fabs(lower(i))) : 0.0;
    double upper_tol = std::isfinite(upper(i))
      ? 1e-8 * (1.0 + std::fabs(upper(i))) : 0.0;
    bool at_lower = std::isfinite(lower(i)) &&
      x(i) <= lower(i) + lower_tol;
    bool at_upper = std::isfinite(upper(i)) &&
      x(i) >= upper(i) - upper_tol;
    bool fixed = at_lower && at_upper;
    bool strictly_blocked = (at_lower && gradient(i) > gradient_tolerance) ||
      (at_upper && gradient(i) < -gradient_tolerance);
    if (!fixed && !strictly_blocked) {
      cone_indices(cone_dimension) = i;
      cone_sign(cone_dimension) = at_lower ? 1 : (at_upper ? -1 : 0);
      ++cone_dimension;
    }
  }
  if (cone_dimension == 0) {
    return true;
  }
  cone_indices.resize(cone_dimension);
  cone_sign.resize(cone_dimension);
  arma::mat cone_curvature = curvature.submat(cone_indices, cone_indices);
  cone_curvature = 0.5 * (cone_curvature + cone_curvature.t());
  arma::vec eigenvalues;
  arma::mat eigenvectors;
  if (!arma::eig_sym(eigenvalues, eigenvectors, cone_curvature) ||
      !eigenvalues.is_finite() || !eigenvectors.is_finite()) {
    return false;
  }
  double tolerance = std::sqrt(std::numeric_limits<double>::epsilon()) *
    std::max(1.0, arma::norm(cone_curvature, "fro"));
  if (eigenvalues.min() >= -tolerance) {
    return true;
  }

  double spectral_scale = arma::abs(eigenvalues).max();
  for (arma::uword i = 0; i < eigenvalues.n_elem; ++i) {
    if (eigenvalues(i) >= -tolerance) {
      continue;
    }
    if (negative_curvature_from_seed(
          cone_curvature, cone_sign, eigenvectors.col(i), tolerance,
          spectral_scale) ||
        negative_curvature_from_seed(
          cone_curvature, cone_sign, -eigenvectors.col(i), tolerance,
          spectral_scale)) {
      return false;
    }
  }

  for (arma::uword i = 0; i < cone_dimension; ++i) {
    arma::vec seed(cone_dimension, arma::fill::zeros);
    seed(i) = cone_sign(i) < 0 ? -1.0 : 1.0;
    if (negative_curvature_from_seed(
          cone_curvature, cone_sign, seed, tolerance, spectral_scale)) {
      return false;
    }
    if (cone_sign(i) == 0) {
      if (negative_curvature_from_seed(
            cone_curvature, cone_sign, -seed, tolerance, spectral_scale)) {
        return false;
      }
    }
  }
  return true;
}

static arma::mat damped_bfgs_update(const arma::mat& current,
                                    const arma::vec& step,
                                    const arma::vec& gradient_change,
                                    bool first_update) {
  arma::mat model = current;
  double sy = arma::dot(step, gradient_change);
  double yy = arma::dot(gradient_change, gradient_change);
  if (first_update && std::isfinite(sy) && std::isfinite(yy) &&
      sy > std::numeric_limits<double>::epsilon() && yy > 0.0) {
    double scale = yy / sy;
    scale = std::max(1e-8, std::min(scale, 1e12));
    model.eye(step.n_elem, step.n_elem);
    model *= scale;
  }

  arma::vec model_step = model * step;
  double sBs = arma::dot(step, model_step);
  if (!std::isfinite(sBs) || sBs <= std::numeric_limits<double>::epsilon()) {
    return positive_definite_model(model);
  }

  double theta = 1.0;
  if (!std::isfinite(sy) || sy < 0.2 * sBs) {
    double denominator = sBs - sy;
    if (!std::isfinite(denominator) || denominator <= 0.0) {
      return positive_definite_model(model);
    }
    theta = 0.8 * sBs / denominator;
  }
  arma::vec corrected = theta * gradient_change +
    (1.0 - theta) * model_step;
  double sr = arma::dot(step, corrected);
  if (!std::isfinite(sr) || sr <= std::numeric_limits<double>::epsilon()) {
    return positive_definite_model(model);
  }

  model -= (model_step * model_step.t()) / sBs;
  model += (corrected * corrected.t()) / sr;
  return positive_definite_model(model);
}

// [[Rcpp::export]]
List trust_region_cpp(
    Function objfun,
    NumericVector parinit,
    double rinit,
    double rmax,
    int iterlim,
    bool minimize,
    double tol,
    bool relative_scale,
    double eta,
    NumericVector lower,
    NumericVector upper,
    Nullable<NumericMatrix> initial_model = R_NilValue,
    double initial_delta = -1.0
) {
  arma::vec x = as<arma::vec>(parinit);
  if (x.n_elem == 0) {
    stop("'parinit' must be a numeric vector with positive length.");
  }

  arma::vec lower_vec = as<arma::vec>(lower);
  arma::vec upper_vec = as<arma::vec>(upper);
  if (lower_vec.n_elem != x.n_elem || upper_vec.n_elem != x.n_elem) {
    stop("Bounds must have the same length as 'parinit'.");
  }
  for (arma::uword i = 0; i < x.n_elem; ++i) {
    if (std::isfinite(lower_vec(i)) && std::isfinite(upper_vec(i)) &&
        lower_vec(i) > upper_vec(i)) {
      stop("Each lower bound must be less than or equal to its upper bound.");
    }
  }
  apply_bounds(x, lower_vec, upper_vec);

  double direction = minimize ? 1.0 : -1.0;
  bool continuing = initial_model.isNotNull();
  double delta = (continuing && std::isfinite(initial_delta) &&
                  initial_delta > 0.0) ? initial_delta : rinit;
  delta = std::min(std::max(delta, 1e-12), rmax);
  int npar = static_cast<int>(x.n_elem);

  vntrs::ObjectiveComponents components = vntrs::parse_objective(
    objfun, x, /*need_gradient=*/true, lower_vec, upper_vec
  );
  double value = components.value;
  arma::vec gradient_original = components.gradient;
  bool exact_hessian = components.hessian_supplied;
  arma::mat supplied_curvature;
  if (exact_hessian) {
    supplied_curvature = direction * components.hessian;
  }
  arma::vec gradient = direction * gradient_original;
  arma::mat model;
  if (exact_hessian) {
    model = positive_definite_model(direction * components.hessian);
  } else if (continuing) {
    NumericMatrix model_matrix(initial_model);
    if (model_matrix.nrow() != npar || model_matrix.ncol() != npar) {
      stop("'initial_model' must have dimension matching 'parinit'.");
    }
    model = as<arma::mat>(model_matrix);
    if (!model.is_finite()) {
      stop("'initial_model' must contain finite values.");
    }
    model = positive_definite_model(model);
  } else if (components.gradient_supplied) {
    model = positive_definite_model(
      direction * vntrs::approximate_hessian_from_gradient(
        objfun, x, components.gradient, lower_vec, upper_vec
      )
    );
  } else {
    double initial_scale = std::max(
      1.0, std::min(1e8, arma::norm(gradient, 2) / std::max(delta, 1.0))
    );
    model = initial_scale * arma::eye<arma::mat>(npar, npar);
  }

  int iter = 0;
  int accepted_updates = 0;
  bool converged = false;

  while (iter < iterlim) {
    double gradient_measure = vntrs::gradient_measure(
      gradient, x, value, lower_vec, upper_vec, relative_scale
    );
    if (!std::isfinite(gradient_measure)) {
      break;
    }
    if (gradient_measure <= tol) {
      converged = true;
      break;
    }
    ++iter;

    model = positive_definite_model(model);
    arma::mat chol_factor;
    bool factorized = arma::chol(chol_factor, model);
    arma::vec newton_step(npar, arma::fill::value(NA_REAL));
    if (factorized) {
      arma::vec y = arma::solve(
        arma::trimatl(chol_factor.t()), -gradient,
        arma::solve_opts::fast
      );
      newton_step = arma::solve(
        arma::trimatu(chol_factor), y,
        arma::solve_opts::fast
      );
    }

    double gradient_norm = arma::norm(gradient, 2);
    if (!std::isfinite(gradient_norm) || gradient_norm == 0.0) {
      break;
    }
    double gBg = arma::as_scalar(gradient.t() * model * gradient);
    double alpha = delta / gradient_norm;
    if (std::isfinite(gBg) && gBg > 0.0) {
      alpha = std::min(
        arma::dot(gradient, gradient) / gBg, delta / gradient_norm
      );
    }
    arma::vec cauchy_step = -alpha * gradient;

    arma::vec step = cauchy_step;
    double newton_norm = arma::norm(newton_step, 2);
    if (newton_step.is_finite() && std::isfinite(newton_norm)) {
      if (newton_norm <= delta) {
        step = newton_step;
      } else {
        arma::vec difference = newton_step - cauchy_step;
        double a = arma::dot(difference, difference);
        double b = 2.0 * arma::dot(cauchy_step, difference);
        double c = arma::dot(cauchy_step, cauchy_step) - delta * delta;
        double tau = 0.0;
        if (a > std::numeric_limits<double>::epsilon()) {
          double discriminant = std::max(0.0, b * b - 4.0 * a * c);
          tau = (-b + std::sqrt(discriminant)) / (2.0 * a);
          tau = std::max(0.0, std::min(1.0, tau));
        }
        step = cauchy_step + tau * difference;
      }
    }

    arma::vec candidate = x + step;
    apply_bounds(candidate, lower_vec, upper_vec);
    step = candidate - x;
    double step_norm = arma::norm(step, 2);
    if (!std::isfinite(step_norm) ||
        step_norm <= std::numeric_limits<double>::epsilon() *
          std::max(1.0, arma::norm(x, 2))) {
      break;
    }

    double predicted = -(
      arma::dot(gradient, step) +
      0.5 * arma::as_scalar(step.t() * model * step)
    );
    if (!std::isfinite(predicted) || predicted <= 0.0) {
      predicted = -arma::dot(gradient, step);
    }

    vntrs::ObjectiveComponents next = vntrs::parse_objective(
      objfun, candidate, /*need_gradient=*/false, lower_vec, upper_vec,
      /*allow_nonfinite=*/true
    );
    if (!std::isfinite(next.value)) {
      delta = std::max(delta / 4.0, 1e-12);
      continue;
    }
    double actual = direction * (value - next.value);
    double rho = -std::numeric_limits<double>::infinity();
    if (std::isfinite(predicted) && predicted > 0.0) {
      rho = actual / predicted;
    }
    if (std::isnan(rho)) {
      rho = -std::numeric_limits<double>::infinity();
    } else if (rho == std::numeric_limits<double>::infinity()) {
      rho = std::numeric_limits<double>::max();
    }

    if (rho < 0.25) {
      delta = std::max(delta / 4.0, 1e-12);
    } else if (rho > 0.75 && step_norm >= 0.9 * delta) {
      delta = std::min(2.0 * delta, rmax);
    }

    if (rho > eta) {
      if (!next.gradient_supplied) {
        next.gradient = vntrs::approximate_gradient(
          objfun, candidate, next.value, lower_vec, upper_vec
        );
      }
      if (!next.gradient.is_finite()) {
        stop("Function 'f' must return finite gradient values.");
      }
      arma::vec gradient_new = direction * next.gradient;
      if (next.hessian_supplied) {
        exact_hessian = true;
        supplied_curvature = direction * next.hessian;
        model = positive_definite_model(supplied_curvature);
      } else if (!exact_hessian) {
        model = damped_bfgs_update(
          model,
          step,
          gradient_new - gradient,
          accepted_updates == 0 && !components.gradient_supplied && !continuing
        );
      }
      x = candidate;
      value = next.value;
      gradient_original = next.gradient;
      gradient = gradient_new;
      ++accepted_updates;
    }

    if (delta <= 1e-12 && rho <= eta) {
      break;
    }
  }

  double gradient_measure = vntrs::gradient_measure(
    gradient, x, value, lower_vec, upper_vec, relative_scale
  );
  converged = converged || (
    std::isfinite(gradient_measure) && gradient_measure <= tol
  );
  if (converged && exact_hessian && !second_order_satisfied(
        supplied_curvature, gradient, x, lower_vec, upper_vec)) {
    converged = false;
  }

  double refinement_trigger = std::pow(
    std::max(tol, std::numeric_limits<double>::epsilon()), 0.25
  );
  if (components.gradient_supplied && !exact_hessian &&
      std::isfinite(gradient_measure) &&
      (converged || gradient_measure <= refinement_trigger)) {
    bool curvature_ok = true;
    for (int refinement = 0; refinement < 5; ++refinement) {
      arma::mat refinement_curvature = direction *
        vntrs::approximate_hessian_from_gradient(
          objfun, x, gradient_original, lower_vec, upper_vec
        );
      curvature_ok = second_order_satisfied(
        refinement_curvature, gradient, x, lower_vec, upper_vec
      );
      if (!curvature_ok) {
        break;
      }
      model = positive_definite_model(refinement_curvature);
      if (gradient_measure <= tol) {
        break;
      }

      arma::vec refinement_step;
      bool solved = arma::solve(
        refinement_step,
        model,
        -gradient,
        arma::solve_opts::likely_sympd
      );
      if (!solved || !refinement_step.is_finite()) {
        break;
      }
      double refinement_norm = arma::norm(refinement_step, 2);
      if (!std::isfinite(refinement_norm) || refinement_norm > rmax) {
        break;
      }

      arma::vec candidate = x + refinement_step;
      apply_bounds(candidate, lower_vec, upper_vec);
      vntrs::ObjectiveComponents refined = vntrs::parse_objective(
        objfun, candidate, /*need_gradient=*/true, lower_vec, upper_vec,
        /*allow_nonfinite=*/true
      );
      if (!std::isfinite(refined.value)) break;
      arma::vec refined_gradient = direction * refined.gradient;
      double refined_measure = vntrs::gradient_measure(
        refined_gradient, candidate, refined.value, lower_vec, upper_vec,
        relative_scale
      );
      double objective_change = direction * (value - refined.value);
      double value_slack = 10.0 * std::numeric_limits<double>::epsilon() *
        std::max(1.0, std::fabs(value));
      if (!std::isfinite(refined_measure) ||
          refined_measure >= gradient_measure ||
          objective_change < -value_slack) {
        break;
      }
      x = candidate;
      value = refined.value;
      gradient_original = refined.gradient;
      gradient = refined_gradient;
      gradient_measure = refined_measure;
    }
    converged = curvature_ok && gradient_measure <= tol;
  }

  if (converged && !components.gradient_supplied && !exact_hessian) {
    arma::mat numerical_curvature = direction *
      vntrs::approximate_hessian_from_values(
        objfun, x, value, gradient_original, lower_vec, upper_vec
      );
    if (!second_order_satisfied(
          numerical_curvature, gradient, x, lower_vec, upper_vec)) {
      converged = false;
    }
  }

  return List::create(
    Named("argument") = wrap(x),
    Named("value") = value,
    Named("converged") = converged,
    Named("iterations") = iter,
    Named("gradient") = wrap(gradient_original),
    Named("scaled_gradient") = gradient_measure,
    Named("model") = wrap(model),
    Named("radius") = delta
  );
}

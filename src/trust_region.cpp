#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include "objective_validation.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
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
    double low = lower_vec(i);
    double upp = upper_vec(i);
    if (std::isfinite(low) && std::isfinite(upp) && low > upp) {
      stop("Each lower bound must be less than or equal to its upper bound.");
    }
    if (std::isfinite(low) && x(i) < low) {
      x(i) = low;
    }
    if (std::isfinite(upp) && x(i) > upp) {
      x(i) = upp;
    }
  }

  double direction = minimize ? 1.0 : -1.0;
  double delta = std::max(rinit, 1e-8);
  delta = std::min(delta, rmax);

  int npar = static_cast<int>(x.n_elem);
  vntrs::ObjectiveComponents components = vntrs::parse_objective(
    objfun,
    x,
    /*require_finite_gradient=*/true,
    /*require_finite_hessian=*/false);
  double value = components.value;
  arma::vec grad_orig = components.gradient;
  arma::mat hess_orig = components.hessian;
  if (!hess_orig.is_finite()) {
    hess_orig.eye(npar, npar);
  }

  int iter = 0;

  while (true) {
    arma::vec grad = direction * grad_orig;
    arma::mat hess = direction * hess_orig;
    double grad_norm = std::sqrt(arma::dot(grad, grad));
    if (!std::isfinite(grad_norm)) {
      grad_norm = NA_REAL;
    }

    if (!std::isfinite(grad_norm) || grad_norm < tol || iter >= iterlim) {
      break;
    }

    ++iter;

    arma::mat chol_factor;
    arma::mat sym = 0.5 * (hess + hess.t());
    double jitter = 0.0;
    bool spd_ok = false;

    for (int attempt = 0; attempt < 7; ++attempt) {
      arma::mat candidate = arma::symmatu(sym);
      if (jitter > 0.0) {
        candidate.diag() += jitter;
      }
      if (arma::chol(chol_factor, candidate)) {
        hess = candidate;
        spd_ok = true;
        break;
      }
      jitter = (jitter == 0.0) ? 1e-8 : std::min(jitter * 10.0, 1e8);
    }

    if (!spd_ok) {
      double diag_shift = jitter;
      if (!sym.is_finite()) {
        diag_shift = std::max(1.0, jitter);
      } else {
        double min_diag = sym.diag().min();
        if (!std::isfinite(min_diag)) {
          min_diag = 1.0;
        }
        diag_shift = std::max(jitter, std::abs(min_diag) + 1.0);
      }

      arma::mat fallback = arma::symmatu(sym);
      fallback.diag() += diag_shift;
      if (arma::chol(chol_factor, fallback)) {
        hess = fallback;
        spd_ok = true;
      } else {
        hess = arma::symmatu(sym);
      }
    }

    arma::vec newton_step(x.n_elem, arma::fill::value(NA_REAL));
    if (spd_ok) {
      arma::vec rhs = -grad;
      arma::vec y = arma::solve(arma::trimatl(chol_factor.t()), rhs);
      arma::vec candidate = arma::solve(arma::trimatu(chol_factor), y);
      if (candidate.is_finite()) {
        newton_step = candidate;
      }
    }

    double grad_norm_eps = std::max(grad_norm, std::numeric_limits<double>::min());
    double gBg = arma::as_scalar(grad.t() * hess * grad);
    double alpha;
    if (!std::isfinite(gBg) || gBg <= 0.0) {
      alpha = delta / grad_norm_eps;
    } else {
      alpha = std::min((grad_norm * grad_norm) / gBg, delta / grad_norm_eps);
    }
    arma::vec cauchy_step = -alpha * grad;

    arma::vec step = cauchy_step;
    double newton_norm = std::sqrt(arma::dot(newton_step, newton_step));
    if (!std::isfinite(newton_norm)) {
      newton_norm = NA_REAL;
    }
    bool newton_valid = newton_step.is_finite() && std::isfinite(newton_norm) &&
      newton_norm <= delta;
    if (newton_valid) {
      step = newton_step;
    } else if (!newton_step.is_finite()) {
      step = cauchy_step;
    } else {
      arma::vec diff = newton_step - cauchy_step;
      double a = arma::dot(diff, diff);
      double b = 2.0 * arma::dot(cauchy_step, diff);
      double c = arma::dot(cauchy_step, cauchy_step) - delta * delta;
      double tau = 0.0;
      if (a > std::numeric_limits<double>::epsilon()) {
        double disc = b * b - 4.0 * a * c;
        disc = std::max(disc, 0.0);
        double sqrt_disc = std::sqrt(disc);
        tau = (-b + sqrt_disc) / (2.0 * a);
        tau = std::max(0.0, std::min(1.0, tau));
      }
      step = cauchy_step + tau * diff;
      double step_norm = std::sqrt(arma::dot(step, step));
      if (!std::isfinite(step_norm)) {
        step_norm = NA_REAL;
      }
      if (!std::isfinite(step_norm) || step_norm > delta * (1.0 + 1e-8)) {
        double cauchy_norm = std::sqrt(arma::dot(cauchy_step, cauchy_step));
        if (!std::isfinite(cauchy_norm)) {
          cauchy_norm = NA_REAL;
        }
        if (std::isfinite(cauchy_norm) && cauchy_norm > 0.0) {
          step = cauchy_step * (delta / cauchy_norm);
        } else {
          step.zeros();
        }
      }
    }

    double step_norm = std::sqrt(arma::dot(step, step));
    if (!std::isfinite(step_norm)) {
      step_norm = NA_REAL;
    }

    double quad_term = arma::as_scalar(step.t() * hess * step);
    double predicted = -(arma::dot(grad, step) + 0.5 * quad_term);
    if (!std::isfinite(predicted) || predicted <= 0.0) {
      predicted = std::numeric_limits<double>::min();
    }

    arma::vec candidate = x + step;
    bool projected = false;
    for (arma::uword i = 0; i < candidate.n_elem; ++i) {
      double low = lower_vec(i);
      double upp = upper_vec(i);
      if (std::isfinite(low) && candidate(i) < low) {
        candidate(i) = low;
        projected = true;
      }
      if (std::isfinite(upp) && candidate(i) > upp) {
        candidate(i) = upp;
        projected = true;
      }
    }
    if (projected) {
      step = candidate - x;
      step_norm = std::sqrt(arma::dot(step, step));
      if (!std::isfinite(step_norm)) {
        step_norm = NA_REAL;
      }
      quad_term = arma::as_scalar(step.t() * hess * step);
      predicted = -(arma::dot(grad, step) + 0.5 * quad_term);
      if (!std::isfinite(predicted) || predicted <= 0.0) {
        predicted = std::numeric_limits<double>::min();
      }
    }

    if (!std::isfinite(step_norm) || step_norm < std::numeric_limits<double>::epsilon()) {
      break;
    }

    vntrs::ObjectiveComponents next = vntrs::parse_objective(
      objfun,
      candidate,
      /*require_finite_gradient=*/true,
      /*require_finite_hessian=*/false);
    double value_new = next.value;
    arma::vec grad_new = next.gradient;
    arma::mat hess_new = next.hessian;
    if (!hess_new.is_finite()) {
      hess_new.eye(npar, npar);
    }

    double current_mod = direction * value;
    double candidate_mod = direction * value_new;
    double actual = current_mod - candidate_mod;
    double rho = actual / predicted;
    if (!std::isfinite(rho)) {
      rho = -std::numeric_limits<double>::infinity();
    }

    if (rho < 0.25) {
      delta = std::max(delta / 4.0, 1e-8);
    } else if (rho > 0.75 && std::abs(step_norm - delta) < 1e-8) {
      delta = std::min(2.0 * delta, rmax);
    }

    if (rho > eta) {
      x = candidate;
      value = value_new;
      grad_orig = grad_new;
      hess_orig = hess_new;
    }
  }

  arma::vec final_grad = direction * grad_orig;
  double final_norm = std::sqrt(arma::dot(final_grad, final_grad));
  if (!std::isfinite(final_norm)) {
    final_norm = NA_REAL;
  }

  arma::vec projected = final_grad;
  for (arma::uword i = 0; i < projected.n_elem; ++i) {
    double low = lower_vec(i);
    double upp = upper_vec(i);
    double xi = x(i);
    double lower_tol = 1e-8 * (1.0 + std::fabs(low));
    double upper_tol = 1e-8 * (1.0 + std::fabs(upp));
    bool at_lower = std::isfinite(low) && xi <= low + lower_tol;
    bool at_upper = std::isfinite(upp) && xi >= upp - upper_tol;
    if (at_lower && projected(i) > 0.0) {
      projected(i) = 0.0;
    }
    if (at_upper && projected(i) < 0.0) {
      projected(i) = 0.0;
    }
  }
  double projected_norm = std::sqrt(arma::dot(projected, projected));
  if (!std::isfinite(projected_norm)) {
    projected_norm = NA_REAL;
  }
  bool converged = std::isfinite(projected_norm) && projected_norm < tol;

  return List::create(
    Named("argument") = wrap(x),
    Named("value") = value,
    Named("converged") = converged,
    Named("iterations") = iter
  );
}

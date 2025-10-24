#ifndef VNTRS_OBJECTIVE_VALIDATION_H
#define VNTRS_OBJECTIVE_VALIDATION_H

#include <RcppArmadillo.h>
#include <limits>
#include <cmath>

namespace vntrs {

struct ObjectiveComponents {
  double value;
  arma::vec gradient;
  arma::mat hessian;
};

inline double extract_numeric_scalar(SEXP value_sexp, const char* field_name) {
  if (!Rf_isNumeric(value_sexp) || Rf_length(value_sexp) != 1) {
    Rcpp::stop("Function 'f' must return finite '%s'.", field_name);
  }
  double value = Rcpp::as<double>(value_sexp);
  if (!R_finite(value)) {
    Rcpp::stop("Function 'f' must return finite '%s'.", field_name);
  }
  return value;
}

inline double extract_value(SEXP result) {
  if (Rf_isNumeric(result)) {
    return extract_numeric_scalar(result, "value");
  }
  if (TYPEOF(result) != VECSXP) {
    Rcpp::stop("Function 'f' must return either a numeric value or a list with elements 'value', 'gradient', and 'hessian'.");
  }
  Rcpp::List out(result);
  if (!out.containsElementNamed("value")) {
    Rcpp::stop("Function 'f' must provide element 'value'.");
  }
  return extract_numeric_scalar(out["value"], "value");
}

inline arma::vec extract_gradient(const Rcpp::List& out,
                                  int npar,
                                  bool require_finite_gradient,
                                  bool& available) {
  available = out.containsElementNamed("gradient");
  if (!available) {
    return arma::vec(npar, arma::fill::value(NA_REAL));
  }
  SEXP gradient_sexp = out["gradient"];
  if (!Rf_isNumeric(gradient_sexp)) {
    Rcpp::stop("Function 'f' must return numeric 'gradient'.");
  }
  Rcpp::NumericVector gradient_vec(gradient_sexp);
  if (gradient_vec.size() != npar) {
    Rcpp::stop("Function 'f' must return gradient with length matching 'npar'.");
  }
  arma::vec gradient = Rcpp::as<arma::vec>(gradient_vec);
  if (require_finite_gradient && !gradient.is_finite()) {
    Rcpp::stop("Function 'f' must return finite gradient values.");
  }
  return gradient;
}

inline arma::mat extract_hessian(const Rcpp::List& out,
                                 int npar,
                                 bool require_finite_hessian,
                                 bool& available) {
  available = out.containsElementNamed("hessian");
  if (!available) {
    return arma::mat(npar, npar, arma::fill::value(NA_REAL));
  }
  SEXP hessian_sexp = out["hessian"];
  if (!Rf_isNumeric(hessian_sexp) || !Rf_isMatrix(hessian_sexp)) {
    Rcpp::stop("Function 'f' must return numeric matrix 'hessian'.");
  }
  Rcpp::NumericMatrix hessian_mat(hessian_sexp);
  if (hessian_mat.nrow() != npar || hessian_mat.ncol() != npar) {
    Rcpp::stop("Function 'f' must return Hessian with dimension 'npar' x 'npar'.");
  }
  arma::mat hessian = Rcpp::as<arma::mat>(hessian_mat);
  if (require_finite_hessian && !hessian.is_finite()) {
    Rcpp::stop("Function 'f' must return finite Hessian entries.");
  }
  return hessian;
}

inline arma::vec finite_difference_steps(const arma::vec& point) {
  arma::vec steps(point.n_elem);
  double eps = std::sqrt(std::numeric_limits<double>::epsilon());
  for (arma::uword i = 0; i < point.n_elem; ++i) {
    double scale = std::fabs(point(i));
    if (!std::isfinite(scale)) {
      scale = 1.0;
    }
    double step = eps * (scale + 1.0);
    if (!std::isfinite(step) || step <= std::numeric_limits<double>::min()) {
      step = eps;
    }
    steps(i) = step;
  }
  return steps;
}

inline double evaluate_value(Rcpp::Function f, const arma::vec& point) {
  SEXP result = f(Rcpp::wrap(point));
  return extract_value(result);
}

inline arma::vec approximate_gradient(Rcpp::Function f,
                                      const arma::vec& point,
                                      const arma::vec& steps,
                                      arma::vec& forward_values,
                                      arma::vec& backward_values) {
  arma::vec gradient(point.n_elem);
  forward_values.set_size(point.n_elem);
  backward_values.set_size(point.n_elem);

  for (arma::uword i = 0; i < point.n_elem; ++i) {
    arma::vec x_forward = point;
    arma::vec x_backward = point;
    double step = steps(i);
    x_forward(i) += step;
    x_backward(i) -= step;

    double value_forward = evaluate_value(f, x_forward);
    double value_backward = evaluate_value(f, x_backward);

    forward_values(i) = value_forward;
    backward_values(i) = value_backward;

    gradient(i) = (value_forward - value_backward) / (2.0 * step);
  }

  return gradient;
}

inline arma::mat approximate_hessian(Rcpp::Function f,
                                     const arma::vec& point,
                                     const arma::vec& steps,
                                     double base_value,
                                     const arma::vec& forward_values,
                                     const arma::vec& backward_values) {
  arma::uword n = point.n_elem;
  arma::mat hessian(n, n);
  hessian.zeros();

  for (arma::uword i = 0; i < n; ++i) {
    double step_i = steps(i);
    double diag = (forward_values(i) - 2.0 * base_value + backward_values(i)) /
      (step_i * step_i);
    hessian(i, i) = diag;
  }

  for (arma::uword i = 0; i < n; ++i) {
    for (arma::uword j = i + 1; j < n; ++j) {
      double step_i = steps(i);
      double step_j = steps(j);

      arma::vec x_pp = point;
      arma::vec x_pm = point;
      arma::vec x_mp = point;
      arma::vec x_mm = point;

      x_pp(i) += step_i; x_pp(j) += step_j;
      x_pm(i) += step_i; x_pm(j) -= step_j;
      x_mp(i) -= step_i; x_mp(j) += step_j;
      x_mm(i) -= step_i; x_mm(j) -= step_j;

      double f_pp = evaluate_value(f, x_pp);
      double f_pm = evaluate_value(f, x_pm);
      double f_mp = evaluate_value(f, x_mp);
      double f_mm = evaluate_value(f, x_mm);

      double mixed = (f_pp - f_pm - f_mp + f_mm) / (4.0 * step_i * step_j);
      hessian(i, j) = mixed;
      hessian(j, i) = mixed;
    }
  }

  return hessian;
}

inline ObjectiveComponents parse_objective(Rcpp::Function f,
                                           const arma::vec& point,
                                           bool require_finite_gradient,
                                           bool require_finite_hessian) {
  SEXP result = f(Rcpp::wrap(point));

  ObjectiveComponents components;
  components.value = extract_value(result);
  components.gradient.set_size(point.n_elem);
  components.hessian.set_size(point.n_elem, point.n_elem);

  bool gradient_available = false;
  bool hessian_available = false;
  arma::vec gradient(point.n_elem, arma::fill::value(NA_REAL));
  arma::mat hessian(point.n_elem, point.n_elem, arma::fill::value(NA_REAL));

  if (TYPEOF(result) == VECSXP) {
    Rcpp::List out(result);
    gradient = extract_gradient(out, point.n_elem, require_finite_gradient, gradient_available);
    hessian = extract_hessian(out, point.n_elem, require_finite_hessian, hessian_available);
    if (gradient_available && !gradient.is_finite()) {
      gradient_available = false;
    }
    if (hessian_available && !hessian.is_finite()) {
      hessian_available = false;
    }
  }

  arma::vec steps;
  arma::vec forward_values;
  arma::vec backward_values;

  if (!gradient_available || !hessian_available) {
    steps = finite_difference_steps(point);
  }

  if (!gradient_available) {
    gradient = approximate_gradient(f, point, steps, forward_values, backward_values);
  }

  if (!hessian_available) {
    if (forward_values.n_elem == 0 || backward_values.n_elem == 0) {
      forward_values.set_size(point.n_elem);
      backward_values.set_size(point.n_elem);
      forward_values.fill(NA_REAL);
      backward_values.fill(NA_REAL);
      for (arma::uword i = 0; i < point.n_elem; ++i) {
        arma::vec x_forward = point;
        arma::vec x_backward = point;
        double step = steps(i);
        x_forward(i) += step;
        x_backward(i) -= step;
        forward_values(i) = evaluate_value(f, x_forward);
        backward_values(i) = evaluate_value(f, x_backward);
      }
    }
    hessian = approximate_hessian(f, point, steps, components.value, forward_values, backward_values);
  }

  if (require_finite_gradient && !gradient.is_finite()) {
    Rcpp::stop("Function 'f' must return finite gradient values.");
  }
  if (require_finite_hessian && !hessian.is_finite()) {
    Rcpp::stop("Function 'f' must return finite Hessian entries.");
  }

  components.gradient = gradient;
  components.hessian = hessian;

  return components;
}

} // namespace vntrs

#endif // VNTRS_OBJECTIVE_VALIDATION_H

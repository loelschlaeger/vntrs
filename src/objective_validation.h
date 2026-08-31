#ifndef VNTRS_OBJECTIVE_VALIDATION_H
#define VNTRS_OBJECTIVE_VALIDATION_H

#include <RcppArmadillo.h>
#include <cmath>
#include <limits>

namespace vntrs {

struct ObjectiveComponents {
  double value;
  arma::vec gradient;
  arma::mat hessian;
  bool gradient_supplied;
  bool hessian_supplied;
};

inline double extract_scalar(SEXP x, const char* field,
                             bool require_finite = true) {
  if (!Rf_isNumeric(x) || Rf_length(x) != 1) {
    Rcpp::stop("Function 'f' must return finite '%s'.", field);
  }
  double value = Rcpp::as<double>(x);
  if (require_finite && !R_finite(value)) {
    Rcpp::stop("Function 'f' must return finite '%s'.", field);
  }
  return value;
}

inline double extract_value(SEXP result, bool require_finite = true) {
  if (Rf_isNumeric(result)) {
    return extract_scalar(result, "value", require_finite);
  }
  if (TYPEOF(result) != VECSXP) {
    Rcpp::stop(
      "Function 'f' must return a numeric value or a list with 'value' "
      "and optional 'gradient' and 'hessian'."
    );
  }
  Rcpp::List out(result);
  if (!out.containsElementNamed("value")) {
    Rcpp::stop("Function 'f' must provide element 'value'.");
  }
  return extract_scalar(out["value"], "value", require_finite);
}

inline arma::vec extract_gradient(const Rcpp::List& out, int n, bool& found) {
  found = out.containsElementNamed("gradient");
  if (!found) {
    return arma::vec(n, arma::fill::value(NA_REAL));
  }
  SEXP x = out["gradient"];
  if (!Rf_isNumeric(x) || Rf_length(x) != n) {
    Rcpp::stop("Function 'f' must return numeric 'gradient' of length 'npar'.");
  }
  arma::vec gradient = Rcpp::as<arma::vec>(x);
  if (!gradient.is_finite()) {
    Rcpp::stop("Function 'f' must return finite gradient values.");
  }
  return gradient;
}

inline arma::mat extract_hessian(const Rcpp::List& out, int n, bool& found) {
  found = out.containsElementNamed("hessian");
  if (!found) {
    return arma::mat(n, n, arma::fill::value(NA_REAL));
  }
  SEXP x = out["hessian"];
  if (!Rf_isNumeric(x) || !Rf_isMatrix(x)) {
    Rcpp::stop("Function 'f' must return numeric matrix 'hessian'.");
  }
  Rcpp::NumericMatrix matrix(x);
  if (matrix.nrow() != n || matrix.ncol() != n) {
    Rcpp::stop("Function 'f' must return Hessian with dimension 'npar' x 'npar'.");
  }
  arma::mat hessian = Rcpp::as<arma::mat>(matrix);
  if (!hessian.is_finite()) {
    Rcpp::stop("Function 'f' must return finite Hessian entries.");
  }
  return hessian;
}

inline arma::vec finite_difference_steps(const arma::vec& x,
                                         double power,
                                         double value_scale = 1.0) {
  double epsilon = std::numeric_limits<double>::epsilon();
  double base = std::pow(epsilon, power);
  double value_step = std::pow(
    epsilon * std::max(1.0, std::fabs(value_scale)), power
  );
  arma::vec steps(x.n_elem);
  for (arma::uword i = 0; i < x.n_elem; ++i) {
    double scale = std::isfinite(x(i)) ? std::max(1.0, std::fabs(x(i))) : 1.0;
    steps(i) = std::max(base * scale, value_step);
  }
  return steps;
}

struct Stencil { int direction; double step; };

inline Stencil select_stencil(double x, double step, double lower,
                              double upper, double points = 2.0) {
  double forward = std::isfinite(upper) ? std::max(0.0, upper - x) : INFINITY;
  double backward = std::isfinite(lower) ? std::max(0.0, x - lower) : INFINITY;
  if (forward >= step && backward >= step) {
    return Stencil{0, step};
  }
  int direction = forward >= backward ? 1 : -1;
  double room = direction > 0 ? forward : backward;
  step = std::isfinite(room) ? std::min(step, room / points) : step;
  double minimum = std::numeric_limits<double>::epsilon() *
    std::max(1.0, std::fabs(x));
  return Stencil{direction, step > minimum ? step : 0.0};
}

inline double zero_like(double) { return 0.0; }
inline arma::vec zero_like(const arma::vec& x) {
  return arma::vec(x.n_elem, arma::fill::zeros);
}

template <typename T, typename Evaluate>
T finite_difference(const arma::vec& x, const T& base, arma::uword i,
                    const Stencil& stencil, Evaluate evaluate) {
  if (stencil.step == 0.0) {
    return zero_like(base);
  }
  arma::vec first = x;
  if (stencil.direction == 0) {
    arma::vec second = x;
    first(i) += stencil.step;
    second(i) -= stencil.step;
    return (evaluate(first) - evaluate(second)) / (2.0 * stencil.step);
  }
  arma::vec second = x;
  first(i) += stencil.direction * stencil.step;
  second(i) += stencil.direction * 2.0 * stencil.step;
  T first_change = evaluate(first) - base;
  T second_change = evaluate(second) - base;
  return stencil.direction * (4.0 * first_change - second_change) /
    (2.0 * stencil.step);
}

inline double evaluate_value(Rcpp::Function f, const arma::vec& x) {
  return extract_value(f(Rcpp::wrap(x)));
}

inline arma::vec evaluate_supplied_gradient(Rcpp::Function f,
                                            const arma::vec& x) {
  SEXP result = f(Rcpp::wrap(x));
  extract_value(result);
  if (TYPEOF(result) != VECSXP) {
    Rcpp::stop("Function 'f' did not consistently provide 'gradient'.");
  }
  bool found = false;
  arma::vec gradient = extract_gradient(
    Rcpp::List(result), static_cast<int>(x.n_elem), found
  );
  if (!found) {
    Rcpp::stop("Function 'f' did not consistently provide 'gradient'.");
  }
  return gradient;
}

inline arma::vec approximate_gradient(Rcpp::Function f, const arma::vec& x,
                                      double value, const arma::vec& lower,
                                      const arma::vec& upper) {
  arma::vec gradient(x.n_elem);
  arma::vec steps = finite_difference_steps(x, 1.0 / 3.0, value);
  auto evaluate = [&](const arma::vec& point) { return evaluate_value(f, point); };
  for (arma::uword i = 0; i < x.n_elem; ++i) {
    gradient(i) = finite_difference(
      x, value, i, select_stencil(x(i), steps(i), lower(i), upper(i)), evaluate
    );
  }
  return gradient;
}

inline arma::mat approximate_hessian_from_gradient(
    Rcpp::Function f, const arma::vec& x, const arma::vec& gradient,
    const arma::vec& lower, const arma::vec& upper) {
  arma::mat hessian(x.n_elem, x.n_elem, arma::fill::zeros);
  arma::vec steps = finite_difference_steps(x, 1.0 / 3.0);
  for (arma::uword i = 0; i < x.n_elem; ++i) {
    Stencil stencil = select_stencil(x(i), steps(i), lower(i), upper(i), 1.0);
    if (stencil.direction == 0) {
      stencil.direction = 1;
    }
    if (stencil.step > 0.0) {
      arma::vec shifted = x;
      shifted(i) += stencil.direction * stencil.step;
      hessian.col(i) = stencil.direction *
        (evaluate_supplied_gradient(f, shifted) - gradient) / stencil.step;
    }
  }
  return 0.5 * (hessian + hessian.t());
}

inline arma::mat approximate_hessian_from_values(
    Rcpp::Function f, const arma::vec& x, double value,
    const arma::vec& lower, const arma::vec& upper) {
  arma::vec gradient = approximate_gradient(f, x, value, lower, upper);
  arma::vec steps = finite_difference_steps(x, 0.25, value);
  arma::mat hessian(x.n_elem, x.n_elem, arma::fill::zeros);
  auto evaluate = [&](const arma::vec& point) {
    double point_value = evaluate_value(f, point);
    return approximate_gradient(f, point, point_value, lower, upper);
  };
  for (arma::uword i = 0; i < x.n_elem; ++i) {
    hessian.col(i) = finite_difference(
      x, gradient, i, select_stencil(x(i), steps(i), lower(i), upper(i)), evaluate
    );
  }
  return 0.5 * (hessian + hessian.t());
}

inline ObjectiveComponents parse_objective(Rcpp::Function f,
                                           const arma::vec& x,
                                           bool need_gradient,
                                           const arma::vec& lower,
                                           const arma::vec& upper,
                                           bool allow_nonfinite = false) {
  SEXP result = f(Rcpp::wrap(x));
  ObjectiveComponents out{
    extract_value(result, !allow_nonfinite),
    arma::vec(x.n_elem, arma::fill::value(NA_REAL)),
    arma::mat(x.n_elem, x.n_elem, arma::fill::value(NA_REAL)),
    false,
    false
  };
  if (!std::isfinite(out.value)) return out;
  if (TYPEOF(result) == VECSXP) {
    Rcpp::List list(result);
    out.gradient = extract_gradient(list, x.n_elem, out.gradient_supplied);
    out.hessian = extract_hessian(list, x.n_elem, out.hessian_supplied);
  }
  if (need_gradient && !out.gradient_supplied) {
    out.gradient = approximate_gradient(f, x, out.value, lower, upper);
  }
  if (need_gradient && !out.gradient.is_finite()) {
    Rcpp::stop("Function 'f' must return finite gradient values.");
  }
  return out;
}

} // namespace vntrs

#endif

#ifndef VNTRS_OPTIMIZATION_METRICS_H
#define VNTRS_OPTIMIZATION_METRICS_H

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>

namespace vntrs {

inline arma::vec projected_gradient(const arma::vec& gradient,
                                    const arma::vec& x,
                                    const arma::vec& lower,
                                    const arma::vec& upper) {
  arma::vec projected = gradient;
  for (arma::uword i = 0; i < projected.n_elem; ++i) {
    double lower_tol = std::isfinite(lower(i))
      ? 1e-8 * (1.0 + std::fabs(lower(i))) : 0.0;
    double upper_tol = std::isfinite(upper(i))
      ? 1e-8 * (1.0 + std::fabs(upper(i))) : 0.0;
    bool at_lower = std::isfinite(lower(i)) &&
      x(i) <= lower(i) + lower_tol;
    bool at_upper = std::isfinite(upper(i)) &&
      x(i) >= upper(i) - upper_tol;
    if (at_lower && projected(i) > 0.0) {
      projected(i) = 0.0;
    }
    if (at_upper && projected(i) < 0.0) {
      projected(i) = 0.0;
    }
  }
  return projected;
}

inline double gradient_measure(const arma::vec& gradient,
                               const arma::vec& x,
                               double value,
                               const arma::vec& lower,
                               const arma::vec& upper,
                               bool relative_scale) {
  arma::vec projected = projected_gradient(gradient, x, lower, upper);
  if (!projected.is_finite()) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (!relative_scale) {
    return arma::norm(projected, 2);
  }
  double measure = 0.0;
  for (arma::uword i = 0; i < x.n_elem; ++i) {
    double parameter_scale = std::isfinite(x(i))
      ? std::max(1.0, std::fabs(x(i))) : 1.0;
    measure = std::max(measure, std::fabs(projected(i)) * parameter_scale);
  }
  return measure / std::max(1.0, std::fabs(value));
}

} // namespace vntrs

#endif

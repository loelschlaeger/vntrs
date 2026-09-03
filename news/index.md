# Changelog

## vntrs 0.3.0

- Reorganized optional arguments of
  [`vntrs()`](https://loelschlaeger.de/vntrs/reference/vntrs.md).

- Renamed argument `tolerance` to `identical_tolerance` to clarify that
  it identifies identical optima.

- Made duplicate-optimum detection scale-free by comparing relative
  component-wise parameter differences.

- Added `interruption_gradient_tolerance` to control the gradient-norm
  threshold used by premature interruption.

- Now samples initialization points across finite parameter bounds.
  `init_min` and `init_max` remain the fallback for parameters without
  finite bounds.

- Expanded the parameter documentation with explanations and recommended
  values.

- Improved initialization when the objective is non-finite at a randomly
  generated starting point.

- Now sorts returned optima from best to worst by objective value.

- Reduced objective-function evaluations by delaying numerical gradients
  until a trial step is accepted and by reusing available values and
  gradients. \# vntrs 0.2.2

- Replaced the absolute gradient convergence check with a scaled
  gradient criterion.

- Reworked missing-derivative handling.

- Added damped-BFGS trust-region models and reuse of local curvature
  across continued searches.

## vntrs 0.2.1

CRAN release: 2026-05-07

- Improved documentation.

- Removed renv.

## vntrs 0.2.0

CRAN release: 2025-10-25

- Added optional bound constraints to the trust-region optimizer. These
  are available through the `lower` and `upper` arguments.

- Improved documentation.

- Implemented the algorithm in Rcpp for faster execution.

## vntrs 0.1.1

CRAN release: 2023-12-21

- Small bug fixes.

## vntrs 0.1.0

CRAN release: 2021-10-18

- Initial release.

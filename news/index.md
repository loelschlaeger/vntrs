# Changelog

## vntrs 0.2.2

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

quadratic_function <- function(Q, b = NULL, c = 0) {
  if (is.null(b)) {
    b <- rep(0, nrow(Q))
  }
  function(x) {
    x <- as.numeric(x)
    value <- 0.5 * sum(x * (Q %*% x)) + sum(b * x) + c
    gradient <- Q %*% x + b
    list(
      value = as.numeric(value), gradient = as.numeric(gradient), hessian = Q
    )
  }
}

ill_conditioned_function <- function(scale = 1) {
  function(x) {
    value <- sum(x)
    gradient <- rep(scale, length(x))
    hessian <- matrix(0, nrow = length(x), ncol = length(x))
    list(value = value, gradient = gradient, hessian = hessian)
  }
}

flat_function <- function() {
  function(x) {
    list(value = 0, gradient = rep(0, length(x)), hessian = diag(length(x)))
  }
}

sleepy <- function(x) {
  Sys.sleep(0.01)
  quadratic_function(matrix(2))(x)
}

two_basin <- function(x) {
  stopifnot(is.numeric(x), length(x) == 1)
  if (x <= 1) {
    value <- x^2
    gradient <- 2 * x
  } else {
    value <- 1 + (x - 2)^2
    gradient <- 2 * (x - 2)
  }
  list(
    value = value,
    gradient = gradient,
    hessian = matrix(2, nrow = 1, ncol = 1)
  )
}

call_trust_region <- function(
    objfun, parinit, rinit = 1, rmax = 10, iterlim = 100,
    minimize = TRUE, tol = 1e-6, eta = 0.1,
    lower = rep(-Inf, length(parinit)), upper = rep(Inf, length(parinit))
) {
  .Call(
    `_vntrs_trust_region_cpp`, objfun, as.numeric(parinit), rinit, rmax,
    as.integer(iterlim), minimize, tol, eta, as.numeric(lower),
    as.numeric(upper)
  )
}

test_that("vntrs respects parameter bounds", {
  set.seed(1)
  res <- vntrs(
    f = quadratic_function(diag(2)),
    npar = 2,
    init_runs = 2,
    neighborhoods = 1,
    neighbors = 2,
    iterlim = 5,
    lower = c(-0.5, -0.5),
    upper = c(0.5, 0.5),
    quiet = TRUE
  )
  if (!is.null(res)) {
    expect_true(all(res$p1 >= -0.5 - 1e-8 & res$p1 <= 0.5 + 1e-8))
    expect_true(all(res$p2 >= -0.5 - 1e-8 & res$p2 <= 0.5 + 1e-8))
  }
})

test_that("vntrs marks global optimum", {
  set.seed(3)
  res <- vntrs(
    f = quadratic_function(diag(2)),
    npar = 2,
    init_runs = 1,
    neighborhoods = 1,
    neighbors = 1,
    iterlim = 5,
    quiet = TRUE
  )
  if (!is.null(res)) {
    best_idx <- which.min(res$value)
    expect_true(res$global[best_idx])
  }
})

test_that("vntrs finds minima and maxima", {
  set.seed(1)
  res_min <- vntrs(
    f = quadratic_function(diag(2) * 2),
    npar = 2,
    init_runs = 1,
    init_min = -1,
    init_max = 1,
    init_iterlim = 5,
    neighborhoods = 1,
    neighbors = 1,
    iterlim = 10,
    tolerance = 1e-6,
    quiet = TRUE
  )
  expect_s3_class(res_min, "data.frame")
  expect_equal(colnames(res_min), c("p1", "p2", "value", "global"))
  expect_true(res_min$global[1])
  res_max <- vntrs(
    f = quadratic_function(diag(2) * -2),
    npar = 2,
    init_runs = 1,
    init_min = -1,
    init_max = 1,
    init_iterlim = 5,
    neighborhoods = 1,
    neighbors = 1,
    iterlim = 10,
    tolerance = 1e-6,
    minimize = FALSE,
    quiet = TRUE
  )
  expect_s3_class(res_max, "data.frame")
  expect_equal(colnames(res_max), c("p1", "p2", "value", "global"))
  expect_true(res_max$global[1])
})

test_that("vntrs handles missing optima and time limits", {
  set.seed(1)
  expect_warning(
    result <- vntrs(
      f = ill_conditioned_function(),
      npar = 1,
      init_runs = 1,
      neighborhoods = 1,
      neighbors = 1,
      iterlim = 2,
      quiet = TRUE
    ),
    "No optima found"
  )
  expect_null(result)
  expect_warning(
    vntrs(
      f = sleepy,
      npar = 1,
      init_runs = 1,
      neighborhoods = 1,
      neighbors = 1,
      iterlim = 5,
      time_limit = 0.001,
      quiet = TRUE
    ),
    "time_limit"
  )
})

test_that("vntrs can retain all local optima", {
  set.seed(1)
  base_args <- list(
    f = two_basin,
    npar = 1,
    init_runs = 1,
    init_min = 0,
    init_max = 0,
    init_iterlim = 5,
    neighborhoods = 6,
    neighbors = 4,
    beta = 0.1,
    iterlim = 4,
    tolerance = 1e-6,
    inferior_tolerance = 1e-6,
    quiet = TRUE
  )
  full_res <- do.call(vntrs, c(base_args, list(collect_all = TRUE)))
  expect_s3_class(full_res, "data.frame")
  expect_equal(nrow(full_res), 2L)
  expect_setequal(round(full_res$value, 8), c(0, 1))
  expect_true(any(full_res$global))
})

test_that("trust_region converges for convex quadratic", {
  f <- quadratic_function(diag(2) * 2)
  res <- call_trust_region(
    objfun = f, parinit = c(5, -3), rinit = 1, rmax = 5, iterlim = 100,
    minimize = TRUE
  )
  expect_true(res$converged)
  expect_true(norm(matrix(res$argument), "F") < 1e-3)
})

test_that("trust_region expands and shrinks radius appropriately", {
  f <- quadratic_function(diag(2) * 2)
  res <- call_trust_region(
    objfun = f, parinit = c(5, 5), rinit = 0.05, rmax = 1, iterlim = 50
  )
  expect_true(res$converged)

  misleading <- local({
    count <- 0
    base <- quadratic_function(diag(2) * 2)
    function(x) {
      count <<- count + 1
      out <- base(x)
      if (count == 2) {
        out$value <- out$value + 10
      }
      out
    }
  })
  res_misleading <- call_trust_region(
    objfun = misleading, parinit = c(1, 1), rinit = 0.1, rmax = 1, iterlim = 5
  )
  expect_false(res_misleading$converged)
})

test_that("trust_region validates parinit", {
  expect_error(
    call_trust_region(
      objfun = quadratic_function(matrix(2)), parinit = numeric(0)
    ),
    "parinit"
  )
})

test_that("trust_region respects bound constraints", {
  f <- quadratic_function(diag(2) * 2)
  lower <- c(1, -Inf)
  upper <- c(2, -0.5)
  res <- call_trust_region(
    objfun = f,
    parinit = c(10, 10),
    lower = lower,
    upper = upper,
    iterlim = 50
  )
  expect_true(res$converged)
  expect_equal(res$argument[1], lower[1], tolerance = 1e-8)
  expect_equal(res$argument[2], upper[2], tolerance = 1e-8)
})


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

two_optima <- function(x) {
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
  vntrs:::trust_region_cpp(
    objfun = objfun,
    parinit = as.numeric(parinit),
    rinit = rinit,
    rmax = rmax,
    iterlim = as.integer(iterlim),
    minimize = minimize,
    tol = tol,
    eta = eta,
    lower = as.numeric(lower),
    upper = as.numeric(upper)
  )
}

call_vntrs_cpp <- function(
    f,
    npar = 1,
    minimize = TRUE,
    init_runs = 1,
    init_min = -1,
    init_max = 1,
    init_iterlim = 5,
    neighborhoods = 1,
    neighbors = 1,
    beta = 0.05,
    iterlim = 5,
    tolerance = 1e-6,
    inferior_tolerance = 1e-6,
    interruption_gradient_tolerance = 1e-3,
    has_time_limit = FALSE,
    time_limit = 0,
    lower = rep(-Inf, npar),
    upper = rep(Inf, npar),
    quiet = TRUE,
    collect_all_optima = FALSE,
    gradient_tolerance = 1e-6
) {
  vntrs:::vntrs_cpp(
    f = f,
    npar = as.integer(npar),
    minimize = minimize,
    init_runs = as.integer(init_runs),
    init_min = init_min,
    init_max = init_max,
    init_iterlim = as.integer(init_iterlim),
    neighborhoods = as.integer(neighborhoods),
    neighbors = as.integer(neighbors),
    beta = beta,
    iterlim = as.integer(iterlim),
    tolerance = tolerance,
    inferior_tolerance = inferior_tolerance,
    interruption_gradient_tolerance = interruption_gradient_tolerance,
    has_time_limit = has_time_limit,
    time_limit = time_limit,
    lower = as.numeric(lower),
    upper = as.numeric(upper),
    quiet = quiet,
    collect_all_optima = collect_all_optima,
    gradient_tolerance = gradient_tolerance
  )
}

valid_vntrs_args <- function() {
  list(
    f = quadratic_function(matrix(2)),
    npar = 1,
    init_runs = 1,
    init_min = 0,
    init_max = 0,
    init_iterlim = 1,
    neighborhoods = 1,
    neighbors = 1,
    beta = 0,
    iterlim = 1,
    identical_tolerance = 1e-6,
    inferior_tolerance = 1e-6,
    interruption_gradient_tolerance = 1e-3,
    quiet = TRUE
  )
}

expect_vntrs_error <- function(...) {
  expect_error(do.call(vntrs, utils::modifyList(valid_vntrs_args(), list(...))))
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
  expect_s3_class(res, "data.frame")
  expect_true(all(res$p1 >= -0.5 - 1e-8 & res$p1 <= 0.5 + 1e-8))
  expect_true(all(res$p2 >= -0.5 - 1e-8 & res$p2 <= 0.5 + 1e-8))
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
  expect_s3_class(res, "data.frame")
  best_idx <- which.min(res$value)
  expect_true(res$global[best_idx])
})

test_that("documented examples return the expected optima", {
  set.seed(1)
  rosenbrock <- function(x) 100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
  rosenbrock_result <- vntrs(f = rosenbrock, npar = 2)

  expect_equal(nrow(rosenbrock_result), 1L)
  expect_equal(sum(rosenbrock_result$global), 1L)

  set.seed(1)
  camel <- function(x) {
    (4 - 2.1 * x[1]^2 + x[1]^4 / 3) * x[1]^2 +
      x[1] * x[2] + (-4 + 4 * x[2]^2) * x[2]^2
  }
  camel_result <- vntrs(
    f = camel,
    npar = 2,
    lower = c(-3, -2),
    upper = c(3, 2),
    collect_all = TRUE,
    neighborhoods = 10
  )

  expect_equal(nrow(camel_result), 6L)
  expect_equal(sum(camel_result$global), 2L)
  expect_equal(sum(!camel_result$global), 4L)
})

test_that("vntrs finds minima and maxima", {
  set.seed(1)
  res_min <- vntrs(
    f = quadratic_function(diag(2) * 2),
    npar = 2,
    init_runs = 1,
    init_min = -1,
    init_max = 1,
    neighborhoods = 1,
    neighbors = 1,
    iterlim = 10,
    identical_tolerance = 1e-6,
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
    neighborhoods = 1,
    neighbors = 1,
    iterlim = 10,
    identical_tolerance = 1e-6,
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
    f = two_optima,
    npar = 1,
    init_runs = 1,
    init_min = 0,
    init_max = 0,
    neighborhoods = 6,
    neighbors = 4,
    beta = 0.1,
    iterlim = 4,
    identical_tolerance = 1e-6,
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


test_that("vntrs validates public arguments", {
  expect_vntrs_error(f = 1)
  expect_vntrs_error(npar = 0)
  expect_vntrs_error(minimize = NA)
  expect_vntrs_error(init_runs = 0)
  expect_vntrs_error(init_min = "low")
  expect_vntrs_error(init_max = -1)
  expect_vntrs_error(init_iterlim = 0)
  expect_vntrs_error(neighborhoods = 0)
  expect_vntrs_error(neighbors = 0)
  expect_vntrs_error(beta = -1)
  expect_vntrs_error(iterlim = 0)
  expect_vntrs_error(identical_tolerance = -1)
  expect_vntrs_error(inferior_tolerance = -1)
  expect_vntrs_error(interruption_gradient_tolerance = -1)
  expect_vntrs_error(gradient_tolerance = 0)
  expect_vntrs_error(time_limit = -1)
  expect_vntrs_error(time_limit = 0)
  expect_vntrs_error(lower = c(0, 0))
  expect_vntrs_error(lower = NA_real_)
  expect_vntrs_error(upper = c(0, 0))
  expect_vntrs_error(upper = NA_real_)
  expect_vntrs_error(lower = 1, upper = 0)
  expect_vntrs_error(collect_all = NA)
  expect_vntrs_error(quiet = NA)
})

test_that("objective return values are validated", {
  expect_vntrs_error(f = function(x) NULL)
  expect_vntrs_error(f = function(x) list(gradient = 0, hessian = matrix(1)))
  expect_vntrs_error(f = function(x) c(1, 2))
  expect_vntrs_error(f = function(x) Inf)
  expect_vntrs_error(
    f = function(x) list(value = 1, gradient = "bad", hessian = matrix(1))
  )
  expect_vntrs_error(
    f = function(x) {
      list(value = 1, gradient = numeric(), hessian = matrix(1))
    }
  )
  expect_vntrs_error(
    f = function(x) {
      list(value = 1, gradient = NA_real_, hessian = matrix(1))
    }
  )
  expect_vntrs_error(
    f = function(x) list(value = 1, gradient = 0, hessian = 1)
  )
  expect_vntrs_error(
    f = function(x) {
      list(value = 1, gradient = 0, hessian = matrix(1, 2, 2))
    }
  )
  expect_vntrs_error(
    f = function(x) {
      list(value = 1, gradient = 0, hessian = matrix(NA_real_))
    }
  )
})

test_that("objectives may omit derivatives", {
  numeric_only <- function(x) sum(x^2)
  value_only <- function(x) list(value = sum(x^2))
  gradient_only <- function(x) list(value = sum(x^2), gradient = 2 * x)
  hessian_only <- function(x) {
    list(value = sum(x^2), hessian = matrix(2, nrow = 1, ncol = 1))
  }

  objectives <- list(
    numeric_only,
    value_only,
    gradient_only,
    hessian_only
  )
  for (objective in objectives) {
    res <- call_trust_region(
      objfun = objective,
      parinit = 1,
      rinit = 1,
      rmax = 2,
      iterlim = 50
    )
    expect_true(res$converged)
    expect_equal(res$value, 0, tolerance = 1e-5)
  }
})

test_that("generated Rcpp wrappers call compiled routines", {
  set.seed(11)
  res <- call_vntrs_cpp(
    f = quadratic_function(matrix(2)),
    npar = 1,
    init_min = 0,
    init_max = 0,
    iterlim = 10
  )
  expect_s3_class(res, "data.frame")
  expect_equal(colnames(res), c("p1", "value", "global"))
  expect_true(res$global[1])
})

test_that("trust_region validates bounds and supports maximization", {
  f <- quadratic_function(matrix(2))
  expect_error(call_trust_region(f, parinit = 1, lower = c(0, 0)))
  expect_error(call_trust_region(f, parinit = 1, lower = 1, upper = 0))

  concave <- quadratic_function(matrix(-2))
  res <- call_trust_region(
    objfun = concave,
    parinit = 1,
    minimize = FALSE
  )
  expect_true(res$converged)
  expect_equal(res$value, 0, tolerance = 1e-5)
})

test_that("finite differences handle non-finite starting scales", {
  constant <- function(x) 0
  res <- call_trust_region(
    objfun = constant,
    parinit = Inf,
    iterlim = 1
  )
  expect_true(res$converged)
  expect_equal(res$value, 0)
})

#' Check \code{f}.
#' @description
#' Function that checks \code{f} for \code{\link{search}}.
#' @param f
#' A function that computes value, gradient, and Hessian of the function to be
#' optimized and returns them as a named list with elements \code{value},
#' \code{gradient}, and \code{hessian}. It is required that \code{hessian} is
#' of class \code{matrix}.
#' @param npar
#' The number of parameters of \code{f}.
#' @inheritParams check_controls
#' @return
#' Invisibly a boolean, if \code{TRUE} than \code{f} is proper, \code{FALSE}
#' otherwise.

check_f = function(f, npar, controls) {

  if(!is.function(f))
    stop("'f' must be a function.")

  ### draw random test points
  test_runs = 10
  y = runif(test_runs*npar, controls$init_min, controls$init_max)
  y = matrix(y, nrow = test_runs, ncol = npar)
  y = round(y,1)

  ### perform test
  for(run in seq_len(test_runs)){
    call = paste0("f(",paste(y[run,],collapse=","),")")
    out = try(f(y[run,]))
    if(class(out) == "try-error")
      stop("Could not compute ",call,".")
    if(!is.list(out))
      stop(call," does not return a list.")
    if(is.null(out[["value"]]))
      stop(call, " does not return a list with element 'value'.")
    if(!(is.numeric(out[["value"]]) && length(out[["value"]])==1))
      stop("The element 'value' in the output list of ",call,
           " is not a single numeric value.")
    if(is.null(out[["gradient"]]))
      stop(call, " does not return a list with element 'gradient'.")
    if(!(is.numeric(out[["gradient"]]) && length(out[["gradient"]])==npar))
      stop("The element 'gradient' in the output list of ",call,
           " is not a numeric vector of length 'npar'.")
    if(is.null(out[["hessian"]]))
      stop(call, " does not return a list with element 'hessian'.")
    if(!(is.numeric(out[["hessian"]]) &&
         is.matrix(out[["hessian"]]) &&
         all(dim(out[["gradient"]])==c(npar,npar))))
      stop("The element 'hessian' in the output list of ",call,
           " is not a numeric matrix of dimension 'npar' x 'npar'.")
  }

  return(invisible(TRUE))
}
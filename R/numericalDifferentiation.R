#' @title Calculate gradient or Hessian matrix
#'
#' @description \code{calc_gradient} or \code{calc_hessian} calculates the gradient or Hessian matrix
#'  of the given function at the given point using central difference numerical approximation.
#'  \code{get_gradient} or \code{get_hessian} calculates the gradient or Hessian matrix of the
#'  log-likelihood function at the parameter values of class \code{'gsmar'} object. \code{get_soc}
#'  returns eigenvalues of the Hessian matrix.
#'
#' @inheritParams simulateGSMAR
#' @param x a numeric vector specifying the point where the gradient or Hessian should be calculated.
#' @param fn a function that takes in argument \code{x} as the \strong{first} argument.
#' @param h difference used to approximate the derivatives.
#' @param ... other arguments passed to \code{fn}.
#' @details Especially the functions \code{get_foc} or \code{get_soc} can be used to check whether
#'  the found estimates denote a (local) maximum point, a saddle point or something else.
#' @return Gradient functions return numerical approximation of the gradient, and Hessian functions return
#'  numerical approximation of the Hessian. \code{get_soc} returns eigenvalues of the Hessian matrix, \code{get_foc}
#'  is the same as \code{get_gradient} but named conveniently.
#' @section Warning:
#'  No argument checks!
#' @examples
#'  # Simple function
#'  foo <- function(x) x^2 + x
#'  calc_gradient(x=1, fn=foo)
#'  calc_gradient(x=-0.5, fn=foo)
#'
#'  # More complicated function
#'  foo <- function(x, a, b) a*x[1]^2 - b*x[2]^2
#'  calc_gradient(x=c(1, 2), fn=foo, a=0.3, b=0.1)
#'
#'  # GMAR model:
#'  params12 <- c(0.18281409, 0.92657275, 0.00214552,
#'   0.85725129, 0.68210294, 0.01900299, 0.88342018)
#'  gmar12 <- GSMAR(logVIX, 1, 2, params12)
#'  get_gradient(gmar12)
#'  get_foc(gmar12)
#'  get_hessian(gmar12)
#'  get_soc(gmar12)
#' @export

calc_gradient <- function(x, fn, h=6e-06, ...) {
  fn <- match.fun(fn)
  n <- length(x)
  I <- diag(1, nrow=n, ncol=n)
  vapply(1:n, function(i1) (fn(x + h*I[i1,], ...) - fn(x - h*I[i1,], ...))/(2*h), numeric(1))
}


#' @rdname calc_gradient
#' @export
calc_hessian <- function(x, fn, h=6e-06, ...) {
  fn <- match.fun(fn)
  n <- length(x)
  I <- diag(1, nrow=n, ncol=n)
  Hess <- matrix(ncol=n, nrow=n)
  for(i1 in 1:n) {
    for(i2 in i1:n) {
      dr1 <- (fn(x + h*I[i1,] + h*I[i2,], ...) - fn(x - h*I[i1,] + h*I[i2,], ...))/(2*h)
      dr2 <- (fn(x + h*I[i1,] - h*I[i2,], ...) - fn(x - h*I[i1,] - h*I[i2,], ...))/(2*h)
      Hess[i1, i2] <- (dr1 - dr2)/(2*h)
      Hess[i2, i1] <- Hess[i1, i2] # Take use of symmetry
    }
  }
  Hess
}



#' @rdname calc_gradient
#' @export
get_gradient <- function(gsmar, h=6e-06) {
  check_gsmar(gsmar)
  warn_dfs(gsmar, warn_about="derivs")
  foo <- function(x) {
    loglikelihood(data=gsmar$data, p=gsmar$model$p, M=gsmar$model$M, params=x, model=gsmar$model$model,
                  restricted=gsmar$model$restricted, constraints=gsmar$model$constraints,
                  conditional=gsmar$model$conditional, parametrization=gsmar$model$parametrization,
                  minval = NA)
  }
  calc_gradient(x=gsmar$params, fn=foo, h=h)
}

#' @rdname calc_gradient
#' @export
get_foc <- function(gsmar, h=6e-06) {
  get_gradient(gsmar=gsmar, h=h)
}

#' @rdname calc_gradient
#' @export
get_hessian <- function(gsmar, h=6e-06) {
  check_gsmar(gsmar)
  warn_dfs(gsmar, warn_about="derivs")
  foo <- function(x) {
    loglikelihood(data=gsmar$data, p=gsmar$model$p, M=gsmar$model$M, params=x, model=gsmar$model$model,
                  restricted=gsmar$model$restricted, constraints=gsmar$model$constraints,
                  conditional=gsmar$model$conditional, parametrization=gsmar$model$parametrization,
                  minval = NA)
  }
  calc_hessian(x=gsmar$params, fn=foo, h=h)
}

#' @rdname calc_gradient
#' @export
get_soc <- function(gsmar, h=6e-06) {
  hess <- get_hessian(gsmar, h)
  if(anyNA(hess)) stop("Missing values in the Hessian matrix. Are estimates at the border of the parameter space?")
  eigen(get_hessian(gsmar, h))$value
}

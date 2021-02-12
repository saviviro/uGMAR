#' @title Calculate gradient or Hessian matrix
#'
#' @description \code{calc_gradient} or \code{calc_hessian} calculates the gradient or Hessian matrix
#'  of the given function at the given point using central difference numerical approximation.
#'  \code{get_gradient} (and \code{get_foc}) or \code{get_hessian} calculates the gradient or Hessian matrix of the
#'  log-likelihood function at the parameter values of a class \code{'gsmar'} object.
#'  \code{get_soc} returns eigenvalues of the Hessian matrix.
#'
#' @inheritParams simulateGSMAR
#' @param x a numeric vector specifying the point at which the gradient or Hessian should be evaluated.
#' @param fn a function that takes in the argument \code{x} as the \strong{first} argument.
#' @param h the difference used to approximate the derivatives.
#' @param varying_h a numeric vector with the same length as \code{x} specifying the difference \code{h}
#'  for each dimension separately. If \code{NULL} (default), then the difference given as parameter \code{h}
#'  will be used for all dimensions.
#' @param custom_h same as \code{varying_h} but if \code{NULL} (default), then the difference \code{h} used for differentiating
#'   overly large degrees of freedom parameters is adjusted to avoid numerical problems, and the difference is \code{6e-6} for the other
#'   parameters.
#' @param ... other arguments passed to \code{fn}.
#' @details In particular, the functions \code{get_foc} and \code{get_soc} can be used to check whether
#'  the found estimates denote a (local) maximum point, a saddle point, or something else.
#' @return The gradient functions return numerical approximation of the gradient, and the Hessian functions return
#'  numerical approximation of the Hessian. \code{get_soc} returns eigenvalues of the Hessian matrix, \code{get_foc}
#'  is the same as \code{get_gradient} but named conveniently.
#' @seealso \code{\link{profile_logliks}}
#' @section Warning:
#'  No argument checks!
#' @examples
#' # Simple function
#' foo <- function(x) x^2 + x
#' calc_gradient(x=1, fn=foo)
#' calc_gradient(x=-0.5, fn=foo)
#' calc_hessian(x=2, fn=foo)
#'
#' # More complicated function
#' foo <- function(x, a, b) a*x[1]^2 - b*x[2]^2
#' calc_gradient(x=c(1, 2), fn=foo, a=0.3, b=0.1)
#' calc_hessian(x=c(1, 2), fn=foo, a=0.3, b=0.1)
#'
#' # GMAR model
#' params12 <- c(1.70, 0.85, 0.30, 4.12, 0.73, 1.98, 0.63)
#' gmar12 <- GSMAR(data=simudata, p=1, M=2, params=params12, model="GMAR")
#' get_gradient(gmar12)
#' get_foc(gmar12)
#' get_hessian(gmar12)
#' get_soc(gmar12)
#' @export

calc_gradient <- function(x, fn, h=6e-06, varying_h=NULL, ...) {
  fn <- match.fun(fn)
  n <- length(x)
  I <- diag(1, nrow=n, ncol=n)
  if(is.null(varying_h)) {
    h <- rep(h, times=n)
  } else {
    stopifnot(length(varying_h) == length(x))
    h <- varying_h
  }
  vapply(1:n, function(i1) (fn(x + h[i1]*I[i1,], ...) - fn(x - h[i1]*I[i1,], ...))/(2*h[i1]), numeric(1))
}


#' @rdname calc_gradient
#' @export
calc_hessian <- function(x, fn, h=6e-06, varying_h=NULL, ...) {
  fn <- match.fun(fn)
  n <- length(x)
  I <- diag(1, nrow=n, ncol=n)
  if(is.null(varying_h)) {
    h <- rep(h, times=n)
  } else {
    stopifnot(length(varying_h) == length(x))
    h <- varying_h
  }
  Hess <- matrix(ncol=n, nrow=n)
  for(i1 in 1:n) {
    for(i2 in i1:n) {
      dr1 <- (fn(x + h[i1]*I[i1,] + h[i2]*I[i2,], ...) - fn(x - h[i1]*I[i1,] + h[i2]*I[i2,], ...))/(2*h[i1])
      dr2 <- (fn(x + h[i1]*I[i1,] - h[i2]*I[i2,], ...) - fn(x - h[i1]*I[i1,] - h[i2]*I[i2,], ...))/(2*h[i1])
      Hess[i1, i2] <- (dr1 - dr2)/(2*h[i2])
      Hess[i2, i1] <- Hess[i1, i2] # Take use of symmetry
    }
  }
  Hess
}


#' @rdname calc_gradient
#' @export
get_gradient <- function(gsmar, custom_h=NULL) {
  check_gsmar(gsmar)
  if(is.null(custom_h)) {
    varying_h <- get_varying_h(p=gsmar$model$p, M=gsmar$model$M, params=gsmar$params, model=gsmar$model$model)
  } else {
    stopifnot(length(custom_h) == length(gsmar$params))
    varying_h <- custom_h
  }

  foo <- function(x) {
    loglikelihood(data=gsmar$data, p=gsmar$model$p, M=gsmar$model$M, params=x, model=gsmar$model$model,
                  restricted=gsmar$model$restricted, constraints=gsmar$model$constraints,
                  conditional=gsmar$model$conditional, parametrization=gsmar$model$parametrization,
                  minval = NA)
  }
  calc_gradient(x=gsmar$params, fn=foo, varying_h=varying_h)
}

#' @rdname calc_gradient
#' @export
get_foc <- function(gsmar, custom_h=NULL) {
  get_gradient(gsmar=gsmar, custom_h=custom_h)
}

#' @rdname calc_gradient
#' @export
get_hessian <- function(gsmar, custom_h=NULL) {
  check_gsmar(gsmar)
  if(is.null(custom_h)) {
    varying_h <- get_varying_h(p=gsmar$model$p, M=gsmar$model$M, params=gsmar$params, model=gsmar$model$model)
  } else {
    stopifnot(length(custom_h) == length(gsmar$params))
    varying_h <- custom_h
  }

  foo <- function(x) {
    loglikelihood(data=gsmar$data, p=gsmar$model$p, M=gsmar$model$M, params=x, model=gsmar$model$model,
                  restricted=gsmar$model$restricted, constraints=gsmar$model$constraints,
                  conditional=gsmar$model$conditional, parametrization=gsmar$model$parametrization,
                  minval = NA)
  }
  calc_hessian(x=gsmar$params, fn=foo, varying_h=varying_h)
}

#' @rdname calc_gradient
#' @export
get_soc <- function(gsmar, custom_h=NULL) {
  hess <- get_hessian(gsmar, custom_h=custom_h)
  if(anyNA(hess)) stop("Missing values in the Hessian matrix. Are estimates at the border of the parameter space?")
  eigen(hess)$value
}


#' @title Get differences 'h' which are adjusted for overly large degrees of freedom parameters
#'
#' @description \code{get_varying_h} adjusts differences for overly large degrees of freedom parameters
#'   for finite difference approximation of the derivatives of the log-likelihood function. StMAR and
#'   G-StMAR models are supported.
#'
#' @inheritParams loglikelihood_int
#' @details This function is used for approximating gradient and Hessian of a StMAR or G-StMAR model. Large
#'   degrees of freedom parameters cause significant numerical error if too small differences are used.
#' @return Returns a vector with the same length as \code{params}. For other parameters than degrees
#'   of freedom parameters larger than 100, the differences will be \code{6e-6}. For the large degrees of
#'   freedom parameters, the difference will be \code{signif(df/1000, digits=2)}.

get_varying_h <- function(p, M, params, model) {
  if(model != "GMAR") {
    dfs <- pick_dfs(p=p, M=M, params=params, model=model)
    adj_diffs <- numeric(length(dfs))
    adj_diffs[dfs <= 100] <- 6e-6
    adj_diffs[dfs > 100] <- signif(dfs[dfs > 100]/1000, digits=2)
    varying_h <- c(rep(6e-6, times=length(params) - length(dfs)), adj_diffs)
  } else {
    varying_h <- rep(6e-6, times=length(params))
  }
  varying_h
}

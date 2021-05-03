#' @title Calculate standard errors for estimates of a GMAR, StMAR, or G-StMAR model
#'
#' @description \code{standard_errors} numerically approximates standard errors for the given estimates of GMAR, StMAR, or GStMAR model.
#'
#' @inheritParams loglikelihood_int
#' @param custom_h a numeric vector with the same length as \code{params} specifying the difference 'h' used in finite difference approximation
#'   for each parameter separately. If \code{NULL} (default), then the difference used for differentiating overly large degrees of freedom
#'   parameters is adjusted to avoid numerical problems, and the difference is \code{6e-6} for the other parameters.
#' @inheritParams fitGSMAR
#' @return Returns approximate standard errors of the parameter values in a numeric vector.

standard_errors <- function(data, p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL, conditional=TRUE,
                           parametrization=c("intercept", "mean"), custom_h=NULL, minval) {
  if(missing(minval)) minval <- get_minval(data)

  # Function to differenciate
  fn <- function(params) {
    tryCatch(loglikelihood_int(data=data, p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints,
                      boundaries=TRUE, conditional=conditional, parametrization=parametrization, checks=FALSE, minval=minval),
             error=function(e) minval)
  }

  # Set up the differences used for finite difference approximation
  if(is.null(custom_h)) {
    varying_h <- get_varying_h(p=p, M=M, params=params, model=model)
  } else {
    stopifnot(length(custom_h) == length(params))
    varying_h <- custom_h
  }

  # Numerically approximated Hessian matrix
  Hess <- calc_hessian(x=params, fn=fn, varying_h=varying_h)

  # Inverse of the observed information matrix
  inv_obs_inf <- tryCatch(solve(-Hess), error=function(cond) return(matrix(NA, ncol=length(params), nrow=length(params))))

  # The diagonal of the inverted observed information matrix evaluated at the estimates
  diag_inv_obs_inf <- diag(inv_obs_inf)

  # Standard errors: NA if can't be calculated (because of numerical issues)
  unlist(lapply(diag_inv_obs_inf, function(x) ifelse(is.na(x) | x < 0, NA, sqrt(x))))
}

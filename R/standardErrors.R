#' @title Calculate standard errors for estimates of GMAR, StMAR or GStMAR model
#'
#' @description \code{standardErrors} numerically approximates standard errors for the given estimates of GMAR, StMAR or GStMAR model
#'
#' @inheritParams loglikelihood_int
#' @param custom_h a numeric vector with the same length as \code{params} specifying the difference 'h' used in finite difference approximation
#'   for each parameter separetely. If \code{NULL} (default), then the difference used for differentiating overly large degrees of freedom
#'   parameters is adjusted to avoid numerical problems, and the difference is \code{6e-6} for the other parameters.
#' @inheritParams fitGSMAR
#' @return Approximate standard errors of the parameter values

standardErrors <- function(data, p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL, conditional=TRUE,
                           parametrization=c("intercept", "mean"), custom_h=NULL, minval) {
  if(missing(minval)) minval <- -(10^(ceiling(log10(length(data))) + 1) - 1)

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
  inv_obs_inf <- tryCatch(solve(-Hess), error=function(cond) return(matrix(NA, ncol=d, nrow=d)))

  # Calculate the standard errors if possible: break loop if all calculated and change the difference if not
  diag_inv_obs_inf <- diag(inv_obs_inf)

  # Standard errors: NA if can't be calculated
  unlist(lapply(diag_inv_obs_inf, function(x) ifelse(is.na(x) | x < 0, NA, sqrt(x))))
}

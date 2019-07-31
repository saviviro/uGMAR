#' @title Calculate standard errors for estimates of GMAR, StMAR or GStMAR model
#'
#' @description \code{standardErrors} numerically approximates standard errors for the given estimates of GMAR, StMAR or GStMAR model
#'
#' @inheritParams loglikelihood_int
#' @inheritParams fitGSMAR
#' @return Approximate standard errors of the parameter values

standardErrors <- function(data, p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL, conditional=TRUE,
                           parametrization=c("intercept", "mean"), h=6e-6, minval) {
  if(missing(minval)) minval <- -(10^(ceiling(log10(length(data))) + 1) - 1)

  # Function to differenciate
  fn <- function(params) {
    tryCatch(loglikelihood_int(data=data, p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints,
                      boundaries=TRUE, conditional=conditional, parametrization=parametrization, checks=FALSE, minval=minval),
             error=function(e) minval)
  }
  #differences <- c(6e-06, 1e-04, 0.001)
  #d <- length(params)
  #I <- diag(1, ncol=d, nrow=d) # Indicates which parameter is derivated

  # Numerically approximated Hessian matrix
  Hess <- calc_hessian(x=params, fn=fn, h=h)

  # Inverse of the observed information matrix
  inv_obs_inf <- tryCatch(solve(-Hess), error=function(cond) return(matrix(NA, ncol=d, nrow=d)))

  # Calculate the standard errors if possible: break loop if all calculated and change the difference if not
  diag_inv_obs_inf <- diag(inv_obs_inf)

  # Standard errors: NA if can't be calculated
  unlist(lapply(diag_inv_obs_inf, function(x) ifelse(is.na(x) | x < 0, NA, sqrt(x))))

#  for(h in differences) {
#    Hess <- calc_hessian(x=params, fn=fn, h=h)
#
    # Inverse of the observed information matrix
#    inv_obs_inf <- tryCatch(solve(-Hess), error=function(cond) return(matrix(NA, ncol=d, nrow=d)))

    # Calculate the standard errors if possible: break loop if all calculated and change the difference if not
#    diag_inv_obs_inf <- diag(inv_obs_inf)
#    if((all(!is.na(diag_inv_obs_inf)) & all(diag_inv_obs_inf >= 0)) | h == differences[length(differences)]) {
#      std_errors <- unlist(lapply(diag_inv_obs_inf, function(x) ifelse(is.na(x) | x < 0, NA, sqrt(x))))
#      break;
#    }
#  }
#  std_errors
}

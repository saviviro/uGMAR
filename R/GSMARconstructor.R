#' @title Create object of class 'gsmar' defining a GMAR, StMAR, or G-StMAR model
#'
#' @description \code{GSMAR} creates an S3 object of class \code{'gsmar'} that defines a GMAR, StMAR, or G-StMAR model.
#'
#' @inheritParams fitGSMAR
#' @inheritParams loglikelihood_int
#' @param calc_qresiduals should quantile residuals be calculated? Default is \code{TRUE} iff the model contains data.
#' @param calc_cond_moments should conditional means and variances be calculated? Default is \code{TRUE} iff the model contains data.
#' @param calc_std_errors should approximate standard errors be calculated?
#' @param custom_h A numeric vector with same the length as the parameter vector: i:th element of custom_h is the difference
#'  used in central difference approximation for partial differentials of the log-likelihood function for the i:th parameter.
#'  If \code{NULL} (default), then the difference used for differentiating overly large degrees of freedom parameters
#'  is adjusted to avoid numerical problems, and the difference is \code{6e-6} for the other parameters.
#' @details Models can be built without data, e.q., in order to simulate from the process, but some things such as quantile
#'  residuals and conditional moments can't be calculated without data.
#' @return Returns an object of class \code{'gsmar'} defining the specified GMAR, StMAR, or G-StMAR model. If data is supplied,
#'  the returned object contains (by default) empirical mixing weights, some conditional and unconditional moments, and quantile
#'  residuals. Note that the first p observations are taken as the initial values so the mixing weights, conditional moments, and
#'  quantile residuals start from the p+1:th observation (interpreted as t=1).
#' @seealso \code{\link{fitGSMAR}}, \code{\link{iterate_more}}, \code{\link{add_data}}, \code{\link{stmar_to_gstmar}},
#'  \code{\link{swap_parametrization}}, \code{\link{get_gradient}}, \code{\link{simulate.gsmar}},
#'  \code{\link{predict.gsmar}}, \code{\link{cond_moments}}, \code{\link{uncond_moments}}, \code{\link{LR_test}}, \code{\link{Wald_test}}
#' @inherit is_stationary references
#' @examples
#' # GMAR model without data
#' params22 <- c(0.9, 0.4, 0.2, 0.5, 0.7, 0.5, -0.2, 0.7, 0.7)
#' gmar22 <- GSMAR(p=2, M=2, params=params22, model="GMAR")
#' gmar22
#'
#' # StMAR model, without data
#' params12t <- c(1.38, 0.88, 0.27, 3.8, 0.74, 3.15, 0.8, 300, 3.6)
#' stmar12t <- GSMAR(p=1, M=2, params=params12t, model="StMAR")
#' stmar12t
#'
#' # G-StMAR model with data
#' params42gs <- c(0.04, 1.34, -0.59, 0.54, -0.36, 0.01, 0.06, 1.28, -0.36,
#'                 0.2, -0.15, 0.04, 0.19, 9.75)
#' gstmar42 <- GSMAR(data=M10Y1Y, p=4, M=c(1, 1), params=params42gs,
#'                   model="G-StMAR")
#' gstmar42
#'
#' # Restricted G-StMAR model with data
#' params42gsr <- c(0.13, 0.03, 1.29, -0.4, 0.25, -0.2, 0.03, 0.05, 0.51, 2.76)
#' gstmar42r <- GSMAR(data=M10Y1Y, p=4, M=c(1, 1), params=params42gsr,
#'                    model="G-StMAR", restricted=TRUE)
#' gstmar42r
#' @export

GSMAR <- function(data, p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL, conditional=TRUE,
                  parametrization=c("intercept", "mean"), calc_qresiduals, calc_cond_moments, calc_std_errors=FALSE, custom_h=NULL) {
  # Checks etc
  model <- match.arg(model)
  parametrization <- match.arg(parametrization)
  check_model(model)
  stopifnot(parametrization %in% c("intercept", "mean"))
  check_pM(p=p, M=M, model=model)
  check_constraint_mat(p=p, M=M, restricted=restricted, constraints=constraints)
  parameter_checks(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints)
  npars <- length(params)
  if(!is.null(custom_h)) stopifnot(length(custom_h) == npars)

  # Initialize some statistics (possibly to be calculated)
  qresiduals <- NULL
  regime_cmeans <- NA
  regime_cvars <- NA
  total_cmeans <- NA
  total_cvars <- NA

  # Calculate some of the statistics based on the arguments
  if(missing(calc_qresiduals)) calc_qresiduals <- ifelse(missing(data), FALSE, TRUE)
  if(missing(calc_cond_moments)) calc_cond_moments <- ifelse(missing(data), FALSE, TRUE)
  if(missing(data) || is.null(data)) { # No data -> some statistics cannot be calculated
    if(calc_qresiduals) warning("Quantile residuals can't be calculated without data")
    if(calc_cond_moments) warning("Conditional moments can't be calculated without data")
    data <- NULL
    lok_and_mw <- list(loglik=NA, mw=NA)
    IC <- data.frame(AIC=NA, HQIC=NA, BIC=NA)
  } else { # Data is provided -> calculate log-likelihood, quantile residuals, conditional moments etc
    data <- check_and_correct_data(data=data, p=p)
    lok_and_mw <- loglikelihood_int(data, p, M, params, model=model, restricted=restricted, constraints=constraints,
                                    conditional=conditional, parametrization=parametrization, boundaries=FALSE,
                                    checks=TRUE, to_return="loglik_and_mw", minval=NA)
    if(calc_qresiduals) {
      qresiduals <- quantile_residuals_int(data, p, M, params, model=model, restricted=restricted, constraints=constraints,
                                          parametrization=parametrization)
    }

    if(calc_cond_moments) {
      # Function to calculate conditional moments
      get_cm <- function(to_return) loglikelihood_int(data, p, M, params, model=model, restricted=restricted, constraints=constraints,
                                                      conditional=conditional, parametrization=parametrization, boundaries=FALSE,
                                                      checks=TRUE, to_return=to_return, minval=NA)
      regime_cmeans <- get_cm("regime_cmeans")
      regime_cvars <- get_cm("regime_cvars")
      total_cmeans <- get_cm("total_cmeans")
      total_cvars <- get_cm("total_cvars")
    }

    obs <- ifelse(conditional, length(data) - p, length(data)) # With conditional log-likelihood function, p obs are lost as initial vals
    IC <- get_IC(loglik=lok_and_mw$loglik, npars=npars, obs=obs) # Information criteria
  }

  # Calculate approximate standard errors
  if(calc_std_errors) {
    if(is.null(data)) {
      warning("Approximate standard errors can't be calculated without data")
      std_errors <- rep(NA, npars)
    } else {
      warn_dfs(p=p, M=M, params=params, model=model)
      std_errors <- tryCatch(standard_errors(data=data, p=p, M=M, params=params, model=model, restricted=restricted,
                                            constraints=constraints, parametrization=parametrization, conditional=conditional,
                                            custom_h=custom_h, minval=-(10^(ceiling(log10(length(data))) + 1) - 1)),
                             error=function(e) {
                               warning("Approximate standard errors can't be calculated")
                               rep(NA, npars)
                             })
    }
  } else {
    std_errors <- rep(NA, npars)
  }

  # Returns the GSMAR model
  structure(list(data=data,
                 model=list(p=p,
                            M=M,
                            model=model,
                            restricted=restricted,
                            constraints=constraints,
                            conditional=conditional,
                            parametrization=parametrization),
                 params=params,
                 std_errors=std_errors,
                 mixing_weights=lok_and_mw$mw,
                 regime_cmeans=regime_cmeans,
                 regime_cvars=regime_cvars,
                 total_cmeans=total_cmeans,
                 total_cvars=total_cvars,
                 quantile_residuals=qresiduals,
                 loglik=structure(lok_and_mw$loglik,
                                  class="logLik",
                                  df=npars),
                 IC=IC,
                 uncond_moments=uncond_moments_int(p=p, M=M, params=params, model=model, restricted=restricted,
                                                  constraints=constraints, parametrization=parametrization),
                 all_estimates=NULL,
                 all_logliks=NULL,
                 which_converged=NULL,
                 which_round=NULL),
            class="gsmar")
}


#' @title Add data to object of class 'gsmar' defining a GMAR, StMAR, or G-StMAR model
#'
#' @description \code{add_data} adds or updates data to object of class '\code{gsmar}' that defines a GMAR, StMAR,
#'  or G-StMAR model. Also calculates empirical mixing weights, conditional moments, and quantile residuals accordingly.
#'
#' @inheritParams loglikelihood_int
#' @inheritParams GSMAR
#' @param gsmar a class 'gsmar' object, typically generated by \code{fitGSMAR} or \code{GSMAR}.
#' @return Returns an object of class 'gsmar' defining the GMAR, StMAR, or G-StMAR model with the data added to the model.
#'   If the object already contained data, the data will be updated. Does not modify the 'gsmar' object given as argument!
#' @seealso \code{\link{fitGSMAR}}, \code{\link{GSMAR}}, \code{\link{iterate_more}}, \code{\link{get_gradient}},
#'  \code{\link{get_regime_means}}, \code{\link{swap_parametrization}}, \code{\link{stmar_to_gstmar}}
#' @inherit is_stationary references
#' @examples
#' # G-StMAR model without data
#' params42gs <- c(0.04, 1.34, -0.59, 0.54, -0.36, 0.01, 0.06, 1.28, -0.36,
#'                 0.2, -0.15, 0.04, 0.19, 9.75)
#' gstmar42 <- GSMAR(p=4, M=c(1, 1), params=params42gs, model="G-StMAR")
#' gstmar42
#'
#' # Add data to the model
#' gstmar42 <- add_data(data=M10Y1Y, gsmar=gstmar42)
#' gstmar42
#' @export

add_data <- function(data, gsmar, calc_qresiduals=TRUE, calc_cond_moments=TRUE, calc_std_errors=FALSE, custom_h=NULL) {
  # Check the data and update it the model with the constructor function GSMAR
  check_gsmar(gsmar)
  check_and_correct_data(data=data, p=gsmar$model$p)
  GSMAR(data=data, p=gsmar$model$p, M=gsmar$model$M, params=gsmar$params,
        model=gsmar$model$model, restricted=gsmar$model$restricted, constraints=gsmar$model$constraints,
        conditional=gsmar$model$conditional, parametrization=gsmar$model$parametrization,
        calc_qresiduals=calc_qresiduals, calc_cond_moments=calc_cond_moments,
        calc_std_errors=calc_std_errors, custom_h=custom_h)
}



#' @title Swap the parametrization of object of class 'gsmar' defining a GMAR, StMAR, or G-StMAR model
#'
#' @description \code{swap_parametrization} swaps the parametrization of object of class '\code{gsmar}'
#'  to \code{"mean"} if the current parametrization is \code{"intercept"}, and vice versa.
#'
#' @inheritParams add_data
#' @details \code{swap_parametrization} is a convenient tool if you have estimated the model in
#'  "intercept"-parametrization but wish to work with "mean"-parametrization in the future,
#'  or vice versa. For example, approximate standard errors are readily available for
#'  parametrized parameters only.
#' @inherit GSMAR references return
#' @inherit add_data seealso
#' @examples
#' # G-StMAR model with intercept parametrization
#' params42gs <- c(0.04, 1.34, -0.59, 0.54, -0.36, 0.01, 0.06, 1.28, -0.36,
#'                 0.2, -0.15, 0.04, 0.19, 9.75)
#' gstmar42 <- GSMAR(data=M10Y1Y, p=4, M=c(1, 1), params=params42gs,
#'                   model="G-StMAR")
#' summary(gstmar42)
#'
#' # Swap to mean parametrization
#' gstmar42 <- swap_parametrization(gstmar42)
#' summary(gstmar42)
#' @export

swap_parametrization <- function(gsmar, calc_std_errors=TRUE, custom_h=NULL) {
  check_gsmar(gsmar)

  # Swap the parametrization to the parameter vector
  change_to <- ifelse(gsmar$model$parametrization == "intercept", "mean", "intercept")
  new_params <- change_parametrization(p=gsmar$model$p, M=gsmar$model$M, params=gsmar$params,
                                       restricted=gsmar$model$restricted, constraints=gsmar$model$constraints,
                                       change_to=change_to)

  # Build the model with the new parametrization
  GSMAR(data=gsmar$data, p=gsmar$model$p, M=gsmar$model$M, params=new_params,
        model=gsmar$model$model, restricted=gsmar$model$restricted,
        constraints=gsmar$model$constraints, conditional=gsmar$model$conditional,
        parametrization=change_to, calc_std_errors=calc_std_errors, custom_h=custom_h)
}


#' @title Estimate a G-StMAR model based on a StMAR model with large degrees of freedom parameters
#'
#' @description \code{stmar_to_gstmar} estimates a G-StMAR model based on a StMAR model with large degree
#'  of freedom parameters.
#'
#' @inheritParams GSMAR
#' @inheritParams add_data
#' @inheritParams stmarpars_to_gstmar
#' @param estimate set \code{TRUE} if the new model should be estimated with a variable metric algorithm
#'  using the StMAR model parameter value as the initial value. By default \code{TRUE} iff the model
#'  contains data.
#' @param calc_std_errors set \code{TRUE} if the approximate standard errors should be calculated.
#'  By default \code{TRUE} iff the model contains data.
#'@param maxit the maximum number of iterations for the variable metric algorithm. Ignored if \code{estimate==FALSE}.
#' @details If a StMAR model contains large estimates for the degrees of freedom parameters,
#'   one should consider switching to the corresponding G-StMAR model that lets the corresponding
#'   regimes to be GMAR type. \code{stmar_to_gstmar} does this switch conveniently.
#' @inherit GSMAR references return
#' @inherit add_data seealso
#' @examples
#' \donttest{
#'  # These are long running example that take approximately 15 seconds to run.
#'  fit42t <- fitGSMAR(data=M10Y1Y, p=4, M=2, model="StMAR", ncalls=1, seeds=6)
#'  fit42t # Overly large degrees of freedom estimate!
#'
#'  # Switch to the appropriate G-StMAR model:
#'  fit42gs <- stmar_to_gstmar(fit42t)
#'  fit42gs
#' }
#' @export

stmar_to_gstmar <- function(gsmar, maxdf=100, estimate, calc_std_errors, maxit=100, custom_h=NULL) {
  # Checks etc
  if(!gsmar$model$model %in% c("StMAR", "G-StMAR")) stop("Only StMAR and G-StMAR models are supported")
  if(missing(estimate)) estimate <- ifelse(is.null(gsmar$data), FALSE, TRUE)
  if(missing(calc_std_errors)) calc_std_errors <- ifelse(is.null(gsmar$data), FALSE, TRUE)
  if(estimate & is.null(gsmar$data)) stop("Can't estimate the model without data")

  # Delete the degrees of freedom parameters of regimes with >maxdf values and obtain the G-StMAR parameter vector accordingly
  new_params <- stmarpars_to_gstmar(p=gsmar$model$p, M=gsmar$model$M, params=gsmar$params,
                                    model=gsmar$model$model, restricted=gsmar$model$restricted,
                                    constraints=gsmar$model$constraints, maxdf=maxdf)
  if(is.null(gsmar$model$constraints) || gsmar$model$restricted) {
    new_constraints <- gsmar$model$constraints
  } else {
    new_constraints <- gsmar$model$constraints[new_params$reg_order] # Reorder the constrnt mats if the ordering of the regimes changed
  }

  # The order M of the new model
  if(new_params$M[2] == 0) {
    new_model <- "GMAR"
    new_M <- new_params$M[1]
  } else if(new_params$M[1] == 0) {
    new_model <- "StMAR"
    new_M <- new_params$M[2]
  } else {
    new_model <- "G-StMAR"
    new_M <- new_params$M
  }

  if(estimate) {
    # Build and estimate the G-StMAR model with a variable metric algorithm using the modified parameter vector as the initial value
    tmp_mod <- GSMAR(data=gsmar$data, p=gsmar$model$p, M=new_M, params=new_params$params,
                     model=new_model, restricted=gsmar$model$restricted,
                     constraints=new_constraints, conditional=gsmar$model$conditional,
                     parametrization=gsmar$model$parametrization, calc_qresiduals=FALSE,
                     calc_cond_moments=FALSE, calc_std_errors=FALSE, custom_h=custom_h)
    new_mod <- iterate_more(tmp_mod, maxit=maxit, custom_h=custom_h)
  } else { # Build the G-StMAR model without estimation
    new_mod <- GSMAR(data=gsmar$data, p=gsmar$model$p, M=new_M, params=new_params$params,
                     model=new_model, restricted=gsmar$model$restricted,
                     constraints=new_constraints, conditional=gsmar$model$conditional,
                     parametrization=gsmar$model$parametrization, calc_qresiduals=FALSE,
                     calc_cond_moments=FALSE, calc_std_errors=calc_std_errors,
                     custom_h=custom_h)
  }
  new_mod
}


#' @title Construct a GSMAR model based on results from an arbitrary estimation round of \code{fitGSMAR}
#'
#' @description \code{alt_gsmar} constructs a GSMAR model based on results from an arbitrary estimation round of \code{fitGSMAR}.
#'
#' @inheritParams add_data
#' @inheritParams GSMAR
#' @param which_round based on which estimation round should the model be constructed? An integer value in 1,...,\code{ncalls}.
#' @param which_largest based on estimation round with which largest log-likelihood should the model be constructed?
#'   An integer value in 1,...,\code{ncalls}. For example, \code{which_largest=2} would take the second largest log-likelihood
#'   and construct the model based on the corresponding estimates. If specified, then \code{which_round} is ignored.
#' @details It's sometimes useful to examine other estimates than the one with the highest log-likelihood value. This function
#'   is just a simple wrapper to \code{GSMAR} that picks the correct estimates from an object returned by \code{fitGSMAR}.
#'
#'   In addition to the S3 methods listed under the topic "Methods (by generic)", the \code{predict} and \code{simulate} methods
#'   are also available for the class 'gsmar' objects (see \code{?predict.gsmar} and \code{?simulate.gsmar}).
#' @inherit GSMAR references return
#' @inherit add_data seealso
#' @examples
#' \donttest{
#'  # These are long running examples that take approximately ...
#'  fit42t <- fitGSMAR(data=M10Y1Y, p=4, M=2, model="StMAR", ncalls=2,
#'                     seeds=c(1, 6))
#'  fit42t # Bad estimate in the boundary of the stationarity region!
#'
#'  # So we build a model based on the next-best local maximum point:
#'  fit42t_alt <- alt_gsmar(fit42t, which_largest=2)
#'  fit42t_alt # Overly large degrees of freedom paramter estimate
#'
#'  # Switch to the appropriate G-StMAR model:
#'  fit42gs <- stmar_to_gstmar(fit42t_alt)
#'  fit42gs
#' }
#' @export

alt_gsmar <- function(gsmar, which_round=1, which_largest, calc_qresiduals=TRUE, calc_cond_moments=TRUE, calc_std_errors=TRUE,
                      custom_h=NULL) {
  stopifnot(!is.null(gsmar$all_estimates))
  stopifnot(which_round >= 1 && which_round <= length(gsmar$all_estimates))
  if(!missing(which_largest)) { # If which_largest is used, select the correct estimation round
    stopifnot(which_largest >= 1 && which_largest <= length(gsmar$all_estimates))
    which_round <- order(gsmar$all_logliks, decreasing=TRUE)[which_largest]
  }

  # Build the model based on the given estimation round
  ret <- GSMAR(data=gsmar$data, p=gsmar$model$p, M=gsmar$model$M, params=gsmar$all_estimates[[which_round]],
               model=gsmar$model$model, restricted=gsmar$model$restricted, constraints=gsmar$model$constraints,
               conditional=gsmar$model$conditional, parametrization=gsmar$model$parametrization,
               calc_qresiduals=calc_qresiduals, calc_cond_moments=calc_cond_moments, calc_std_errors=calc_std_errors,
               custom_h=custom_h)

  # Obtain the rest of the statistics related to the estimation
  ret$all_estimates <- gsmar$all_estimates
  ret$all_logliks <- gsmar$all_logliks
  ret$which_converged <- gsmar$which_converged
  if(!is.null(gsmar$which_round)) {
    ret$which_round <- which_round
  }
  warn_ar_roots(ret)
  ret
}

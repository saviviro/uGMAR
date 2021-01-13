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
#'  \code{\link{swap_parametrization}}, \code{\link{get_gradient}}, \code{\link{simulateGSMAR}},
#'  \code{\link{predict.gsmar}}, \code{\link{cond_moments}}, \code{\link{uncond_moments}}, \code{\link{LR_test}}, \code{\link{Wald_test}}
#' @inherit is_stationary references
#' @examples
#' # GMAR model without data
#' params12 <- c(0.18, 0.93, 0.01, 0.86, 0.68, 0.02, 0.88)
#' gmar12 <- GSMAR(p=1, M=2, params=params12, model="GMAR")
#' gmar12
#'
#' # StMAR model, without data
#' params12t <- c(1.38, 0.88, 0.27, 3.8, 0.74, 3.15, 0.8, 300, 3.6)
#' stmar12t <- GSMAR(p=1, M=2, params=params12t, model="StMAR")
#' stmar12t
#'
#' # Restricted G-StMAR-model
#' params42gsr <- c(0.11, 0.03, 1.27, -0.39, 0.24, -0.17, 0.03, 1.01, 0.3, 2.03)
#' gstmar42r <- GSMAR(data=T10Y1Y, p=4, M=c(1, 1), params=params42gsr,
#'  model="G-StMAR", restricted=TRUE)
#' gstmar42r
#'
#' # Two-regime GMAR p=2 model with the second AR coeffiecient of
#' # of the second regime contrained to zero.
#' constraints <- list(diag(1, ncol=2, nrow=2), as.matrix(c(1, 0)))
#' params22c <- c(0.03, 1.27, -0.29, 0.03, -0.01, 0.91, 0.34, 0.88)
#' gmar22c <- GSMAR(T10Y1Y, p=2, M=2, params=params22c,
#'  model="GMAR", constraints=constraints)
#' gmar22c
#' @export

GSMAR <- function(data, p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL, conditional=TRUE,
                  parametrization=c("intercept", "mean"), calc_qresiduals, calc_cond_moments, calc_std_errors=FALSE, custom_h=NULL) {
  model <- match.arg(model)
  parametrization <- match.arg(parametrization)
  check_model(model)
  stopifnot(parametrization %in% c("intercept", "mean"))
  check_pM(p=p, M=M, model=model)
  check_constraint_mat(p=p, M=M, restricted=restricted, constraints=constraints)
  parameter_checks(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints)
  npars <- length(params)
  if(!is.null(custom_h)) stopifnot(length(custom_h) == npars)

  qresiduals <- NULL
  regime_cmeans <- NA
  regime_cvars <- NA
  total_cmeans <- NA
  total_cvars <- NA

  if(missing(calc_qresiduals)) calc_qresiduals <- ifelse(missing(data), FALSE, TRUE)
  if(missing(calc_cond_moments)) calc_cond_moments <- ifelse(missing(data), FALSE, TRUE)
  if(missing(data) || is.null(data)) {
    if(calc_qresiduals) warning("Quantile residuals can't be calculated without data")
    if(calc_cond_moments) warning("Conditional moments can't be calculated without data")
    data <- NULL
    lok_and_mw <- list(loglik=NA, mw=NA)
    IC <- data.frame(AIC=NA, HQIC=NA, BIC=NA)
  } else {
    data <- check_and_correct_data(data=data, p=p)
    lok_and_mw <- loglikelihood_int(data, p, M, params, model=model, restricted=restricted, constraints=constraints,
                                    conditional=conditional, parametrization=parametrization, boundaries=FALSE,
                                    checks=TRUE, to_return="loglik_and_mw", minval=NA)
    if(calc_qresiduals) {
      qresiduals <- quantile_residuals_int(data, p, M, params, model=model, restricted=restricted, constraints=constraints,
                                          parametrization=parametrization)
    }

    if(calc_cond_moments) {
      get_cm <- function(to_return) loglikelihood_int(data, p, M, params, model=model, restricted=restricted, constraints=constraints,
                                                      conditional=conditional, parametrization=parametrization, boundaries=FALSE,
                                                      checks=TRUE, to_return=to_return, minval=NA)
      regime_cmeans <- get_cm("regime_cmeans")
      regime_cvars <- get_cm("regime_cvars")
      total_cmeans <- get_cm("total_cmeans")
      total_cvars <- get_cm("total_cvars")
    }

    obs <- ifelse(conditional, length(data) - p, length(data))
    IC <- get_IC(loglik=lok_and_mw$loglik, npars=npars, obs=obs)
  }

  if(calc_std_errors) {
    if(is.null(data)) {
      warning("Approximate standard errors can't be calculated without data")
      std_errors <- rep(NA, npars)
    } else {
      warn_dfs(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints)
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
                 qr_tests=NULL),
            class="gsmar")
}


#' @title Add data to object of class 'gsmar' defining a GMAR, StMAR, or G-StMAR model
#'
#' @description \code{add_data} adds or updates data to object of class '\code{gsmar}' that defines a GMAR, StMAR,
#'  or G-StMAR model. Also calculates empirical mixing weights, conditional moments, and quantile residuals accordingly.
#'
#' @inheritParams simulateGSMAR
#' @inheritParams loglikelihood_int
#' @inheritParams GSMAR
#' @return Returns an object of class 'gsmar' defining the GMAR, StMAR, or G-StMAR model with the data added to the model.
#'   If the object already contained data, the data will be updated. Does not modify the 'gsmar' object given as argument!
#' @seealso \code{\link{fitGSMAR}}, \code{\link{GSMAR}}, \code{\link{iterate_more}}, \code{\link{get_gradient}},
#'  \code{\link{get_regime_means}}, \code{\link{swap_parametrization}}, \code{\link{stmar_to_gstmar}}
#' @inherit is_stationary references
#' @examples
#' # Restricted G-StMAR-model without data
#' params42gsr <- c(0.11, 0.03, 1.27, -0.39, 0.24, -0.17, 0.03, 1.01, 0.3, 2.03)
#' gstmar42r <- GSMAR(p=4, M=c(1, 1), params=params42gsr,
#'  model="G-StMAR", restricted=TRUE)
#' gstmar42r
#'
#' # Add data to the model
#' gstmar42r <- add_data(data=T10Y1Y, gstmar42r)
#' gstmar42r
#' @export

add_data <- function(data, gsmar, calc_qresiduals=TRUE, calc_cond_moments=TRUE, calc_std_errors=FALSE, custom_h=NULL) {
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
#' # Restricted G-StMAR-model with intercept paarametrization
#' params42gsr <- c(0.11, 0.03, 1.27, -0.39, 0.24, -0.17, 0.03, 1.01, 0.3, 2.03)
#' gstmar42r <- GSMAR(data=T10Y1Y, p=4, M=c(1, 1), params=params42gsr,
#'  model="G-StMAR", restricted=TRUE)
#' summary(gstmar42r)
#'
#' # Swap to mean parametrization
#' gstmar42r <- swap_parametrization(gstmar42r)
#' summary(gstmar42r)
#' @export

swap_parametrization <- function(gsmar, calc_std_errors=TRUE, custom_h=NULL) {
  check_gsmar(gsmar)
  change_to <- ifelse(gsmar$model$parametrization == "intercept", "mean", "intercept")
  new_params <- change_parametrization(p=gsmar$model$p, M=gsmar$model$M, params=gsmar$params,
                                       restricted=gsmar$model$restricted, constraints=gsmar$model$constraints,
                                       change_to=change_to)
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
#'  # These are long running examples and use parallel computing
#'  fit43t <- fitGSMAR(T10Y1Y, 4, 3, model="StMAR", ncalls=2, seeds=1:2)
#'  fit43t
#'  fit43gst <- stmar_to_gstmar(fit43t)
#'  fit43gst
#' }
#' @export

stmar_to_gstmar <- function(gsmar, maxdf=100, estimate, calc_std_errors, maxit=100, custom_h=NULL) {
  if(gsmar$model$model != "StMAR") stop("Only StMAR models are supported as input")
  if(missing(estimate)) estimate <- ifelse(is.null(gsmar$data), FALSE, TRUE)
  if(missing(calc_std_errors)) calc_std_errors <- ifelse(is.null(gsmar$data), FALSE, TRUE)
  if(estimate & is.null(gsmar$data)) stop("Can't estimate the model without data")

  new_params <- stmarpars_to_gstmar(p=gsmar$model$p, M=gsmar$model$M, params=gsmar$params,
                                    restricted=gsmar$model$restricted, constraints=gsmar$model$constraints,
                                    maxdf=maxdf)
  if(is.null(gsmar$model$constraints) || gsmar$model$restricted) {
    new_constraints <- gsmar$model$constraints
  } else {
    new_constraints <- gsmar$model$constraints[new_params$reg_order]
  }

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
    tmp_mod <- GSMAR(data=gsmar$data, p=gsmar$model$p, M=new_M, params=new_params$params,
                     model=new_model, restricted=gsmar$model$restricted,
                     constraints=new_constraints, conditional=gsmar$model$conditional,
                     parametrization=gsmar$model$parametrization, calc_qresiduals=FALSE,
                     calc_cond_moments=FALSE, calc_std_errors=FALSE, custom_h=custom_h)
    new_mod <- iterate_more(tmp_mod, maxit=maxit, custom_h=custom_h)
  } else {
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
#' @inheritParams simulateGSMAR
#' @inheritParams GSMAR
#' @param which_round based on which estimation round should the model be constructed? An integer value in 1,...,\code{ncalls}.
#' @param which_largest based on estination round with which largest log-likelihood should the model be constructed?
#'   An integer value in 1,...,\code{ncalls}. For example, \code{which_largest=2} would take the second largest log-likelihood
#'   and construct the model based on the corresponding estimates. If used, then \code{which_round} is ignored.
#' @details It's sometimes useful to examine other estimates than the one with the highest log-likelihood value. This function
#'   is just a simple wrapper to \code{GSMAR} that picks the correct estimates from an object returned by \code{fitGSMAR}.
#' @inherit GSMAR references return
#' @inherit add_data seealso
#' @examples
#' \donttest{
#'  # These are long running examples and use parallel computing
#'  fit43t <- fitGSMAR(T10Y1Y, 4, 3, model="StMAR", ncalls=2, seeds=1:2)
#'  fit43t
#'  fit43t2 <- alt_gsmar(fit43t, which_largest=2)
#'  fit43t2
#' }
#' @export

alt_gsmar <- function(gsmar, which_round=1, which_largest, calc_qresiduals=TRUE, calc_cond_moments=TRUE, calc_std_errors=TRUE,
                      custom_h=NULL) {
  stopifnot(!is.null(gsmar$all_estimates))
  stopifnot(which_round >= 1 || which_round <= length(gsmar$all_estimates))
  if(!missing(which_largest)) {
    stopifnot(which_largest >= 1 || which_largest <= length(gsmar$all_estimates))
    which_round <- order(gsmar$all_logliks, decreasing=TRUE)[which_largest]
  }
  GSMAR(data=gsmar$data, p=gsmar$model$p, M=gsmar$model$M, params=gsmar$all_estimates[[which_round]],
        model=gsmar$model$model, restricted=gsmar$model$restricted, constraints=gsmar$model$constraints,
        conditional=gsmar$model$conditional, parametrization=gsmar$model$parametrization,
        calc_qresiduals=calc_qresiduals, calc_cond_moments=calc_cond_moments, calc_std_errors=calc_std_errors, custom_h=custom_h)
}

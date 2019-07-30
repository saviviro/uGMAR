#' @title Create object of class 'gsmar' defining a GMAR, StMAR or G-StMAR model
#'
#' @description \code{GSMAR} creates an S3 object of class \code{'gsmar'} that defines a GMAR, StMAR or G-StMAR model.
#'
#' @inheritParams loglikelihood_int
#' @param calc_qresiduals should quantile residuals be calculated? Default is \code{TRUE} iff the model contains data.
#' @param calc_cond_moments should conditional means and variances be calculated? Default is \code{TRUE} iff the model contains data.
#' @param calc_std_errors should approximate standard errors be calculated?
#' @details Models can be built without data, e.q., in order to simulate from the process, but some elements such as quantile
#'  residuals and conditional moments can't be calculated without data.
#' @return Returns an object of class \code{'gsmar'} defining the specified GMAR, StMAR or G-StMAR model. If data is suplied, the returned object
#'   contains (by default) empirical mixing weights, conditional means and variances and quantile residuals. Note that the first p observations are
#'   taken as the initial values so mixing weights, conditional moments and qresiduals start from the p+1:th observation (interpreted as t=1).
#' @seealso \code{\link{fitGSMAR}}, \code{\link{iterate_more}}, \code{\link{add_data}}, \code{\link{stmar_to_gstmar}},
#'  \code{\link{swap_parametrization}}, \code{\link{get_gradient}}, \code{\link{simulateGSMAR}},
#'  \code{\link{predict.gsmar}}, \code{\link{condMoments}}, \code{\link{uncondMoments}}
#' @inherit isStationary references
#' @examples
#' # GMAR model
#' params12 <- c(0.18, 0.93, 0.01, 0.86, 0.68, 0.02, 0.88)
#' gmar12 <- GSMAR(data=logVIX, p=1, M=2, params=params12, model="GMAR")
#' gmar12
#'
#' # Restricted GMAR model
#' params12r <- c(0.21, 0.23, 0.92, 0.01, 0.02, 0.86)
#' gmar12r <- GSMAR(data=logVIX, p=1, M=2, params=params12r, model="GMAR",
#'  restricted=TRUE)
#' gmar12r
#'
#' # StMAR model, without data
#' params12t <- c(1.38, 0.88, 0.27, 3.8, 0.74, 3.15, 0.8, 300, 3.6)
#' stmar12t <- GSMAR(p=1, M=2, params=params12t, model="StMAR")
#' stmar12t
#'
#' # G-StMAR model (similar to the StMAR model above), without data
#' params12gs <- c(1.38, 0.88, 0.27, 3.8, 0.74, 3.15, 0.8, 3.6)
#' gstmar12 <- GSMAR(p=1, M=c(1, 1), params=params12gs, model="G-StMAR")
#' gstmar12
#'
#' # Restricted G-StMAR-model
#' params12gsr <- c(0.31, 0.33, 0.88, 0.01, 0.02, 0.77, 2.72)
#' gstmar12r <- GSMAR(data=logVIX, p=1, M=c(1, 1), params=params12gsr,
#'  model="G-StMAR", restricted=TRUE)
#' gstmar12r
#'
#' # GMAR(p=2, M=2) model such that the second AR coefficient of the
#' # second regime is constrained to zero.
#' constraints <- list(diag(1, ncol=2, nrow=2), as.matrix(c(1, 0)))
#' params22c <- c(0.61, 0.83, -0.06, 0.02, 0.21, 0.91, 0.01, 0.16)
#' gmar22c <- GSMAR(logVIX, p=2, M=2, params=params22c,
#'  model="GMAR", constraints=constraints)
#' gmar22c
#'
#' # Such StMAR(3,2) that the AR coefficients are restricted to be
#' # the same for both regimes and that the second AR coefficients are
#' # constrained to zero.
#' params32trc <- c(0.35, 0.33, 0.88, -0.02, 0.01, 0.01, 0.36, 4.53, 1000)
#' stmar32rc <- GSMAR(logVIX, p=3, M=2, params=params32trc, model="StMAR",
#'  restricted=TRUE, constraints=matrix(c(1, 0, 0, 0, 0, 1), ncol=2))
#' stmar32rc
#'
#' # Mixture version of Heterogenuous autoregressive (HAR) model (without data)
#' paramsHAR <- c(1, 0.1, 0.2, 0.3, 1, 2, 0.15, 0.25, 0.35, 2, 0.55)
#' r1 = c(1, rep(0, 21)); r2 = c(rep(0.2, 5), rep(0, 17)); r3 = rep(1/22, 22)
#' R0 = cbind(r1, r2, r3)
#' mixhar <- GSMAR(p=22, M=2, params=paramsHAR, model="GMAR", constraints=list(R0, R0))
#' mixhar
#' @export

GSMAR <- function(data, p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL, conditional=TRUE,
                  parametrization=c("intercept", "mean"), calc_qresiduals, calc_cond_moments, calc_std_errors=FALSE) {
  model <- match.arg(model)
  parametrization <- match.arg(parametrization)
  check_model(model)
  stopifnot(parametrization %in% c("intercept", "mean"))
  checkPM(p=p, M=M, model=model)
  checkConstraintMat(p=p, M=M, restricted=restricted, constraints=constraints)
  parameterChecks(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints)
  npars <- length(params)

  qresiduals <- NULL
  regime_cmeans <- NA
  regime_cvars <- NA
  total_cmeans <- NA
  total_cvars <- NA

  if(missing(calc_qresiduals)) calc_qresiduals <- ifelse(missing(data), FALSE, TRUE)
  if(missing(calc_cond_moments)) calc_cond_moments <- ifelse(missing(data), FALSE, TRUE)
  if(missing(data) || is.null(data)) {
    if(calc_qresiduals == TRUE) warning("Quantile residuals can't be calculated without data")
    if(calc_cond_moments == TRUE) warning("Conditional moments can't be calculated without data")
    data <- NULL
    lok_and_mw <- list(loglik=NA, mw=NA)
    IC <- data.frame(AIC=NA, HQIC=NA, BIC=NA)
  } else {
    data <- checkAndCorrectData(data=data, p=p)
    lok_and_mw <- loglikelihood_int(data, p, M, params, model=model, restricted=restricted, constraints=constraints,
                                    conditional=conditional, parametrization=parametrization, boundaries=FALSE,
                                    checks=TRUE, to_return="loglik_and_mw", minval=NA)
    if(calc_qresiduals == TRUE) {
      qresiduals <- quantileResiduals_int(data, p, M, params, model=model, restricted=restricted, constraints=constraints,
                                          parametrization=parametrization)
    }

    if(calc_cond_moments == TRUE) {
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

  if(calc_std_errors == TRUE) {
    if(is.null(data)) {
      warning("Approximate standard errors can't be calculated without data")
      std_errors <- rep(NA, npars)
    } else {
      warn_dfs(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints, warn_about="errors")
      std_errors <- tryCatch(standardErrors(data=data, p=p, M=M, params=params, model=model, restricted=restricted,
                                            constraints=constraints, parametrization=parametrization, conditional=conditional,
                                            minval=-(10^(ceiling(log10(length(data))) + 1) - 1)),
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
                 uncond_moments=uncondMoments_int(p=p, M=M, params=params, model=model, restricted=restricted,
                                                  constraints=constraints, parametrization=parametrization),
                 all_estimates=NULL,
                 all_logliks=NULL,
                 which_converged=NULL,
                 qr_tests=NULL),
            class="gsmar")
}


#' @title Add data to object of class 'gsmar' defining a GMAR, StMAR or G-StMAR model
#'
#' @description \code{add_data} adds or updates data to object of class '\code{gsmar}' that defines a GMAR, StMAR or G-StMAR
#'  model. Also calculates mixing weights, conditional moments and quantile residuals accordingly.
#'
#' @inheritParams simulateGSMAR
#' @inheritParams loglikelihood_int
#' @inheritParams GSMAR
#' @return Returns an object of class 'gsmar' defining the GMAR, StMAR or G-StMAR model with the data added to the model.
#'   If the object already contained data, the data will be updated. Does not modify the 'gsmar' object given as argument!
#' @seealso \code{\link{fitGSMAR}}, \code{\link{GSMAR}}, \code{\link{iterate_more}}, \code{\link{get_gradient}},
#'  \code{\link{get_regime_means}}, \code{\link{swap_parametrization}}, \code{\link{stmar_to_gstmar}}
#' @inherit isStationary references
#' @examples
#' # GMAR model without data
#' params12 <- c(0.18, 0.93, 0.01, 0.86, 0.68, 0.02, 0.88)
#' gmar12 <- GSMAR(p=1, M=2, params=params12, model="GMAR")
#' gmar12
#'
#' # Add data to the model
#' gmar12 <- add_data(data=logVIX, gmar12)
#' gmar12
#' @export

add_data <- function(data, gsmar, calc_qresiduals=TRUE, calc_cond_moments=TRUE, calc_std_errors=FALSE) {
  check_gsmar(gsmar)
  checkAndCorrectData(data=data, p=gsmar$model$p)
  GSMAR(data=data, p=gsmar$model$p, M=gsmar$model$M, params=gsmar$params,
        restricted=gsmar$model$restricted, constraints=gsmar$model$constraints,
        conditional=gsmar$model$conditional, parametrization=gsmar$model$parametrization,
        calc_qresiduals=calc_qresiduals, calc_cond_moments=calc_cond_moments,
        calc_std_errors=calc_std_errors)
}



#' @title Swap the parametrization of object of class 'gsmar' defining a gsmar model
#'
#' @description \code{swap_parametrization} swaps the parametrization of object of class '\code{gsmar}'
#'  to \code{"mean"} if the current parametrization is \code{"intercept"}, and vice versa.
#'
#' @inheritParams add_data
#' @details \code{swap_parametrization} is convenient tool if you have estimated the model in
#'  "intercept"-parametrization, but wish to work with "mean"-parametrization in the future, or vice versa.
#'  In \code{gsmarkit}, for example the approximate standard errors are only available for
#'  parametrized parameters.
#' @inherit GSMAR references return
#' @inherit add_data seealso
#' @examples
#' # GMAR model with intercept parametrization
#' params12 <- c(0.18, 0.93, 0.01, 0.86, 0.68, 0.02, 0.88)
#' gmar12 <- GSMAR(data=logVIX, p=1, M=2, params=params12, model="GMAR")
#' gmar12
#'
#' # Swap to mean parametrization
#' gmar12 <- swap_parametrization(gmar12)
#' gmar12
#' @export

swap_parametrization <- function(gsmar, calc_std_errors=TRUE) {
  check_gsmar(gsmar)
  change_to <- ifelse(gsmar$model$parametrization == "intercept", "mean", "intercept")
  new_params <- change_parametrization(p=gsmar$model$p, M=gsmar$model$M, params=gsmar$params,
                                       restricted=gsmar$model$restricted, constraints=gsmar$model$constraints,
                                       change_to=change_to)
  GSMAR(data=gsmar$data, p=gsmar$model$p, M=gsmar$model$M, params=new_params,
        model=gsmar$model$model, restricted=gsmar$model$restricted,
        constraints=gsmar$model$constraints, conditional=gsmar$model$conditional,
        parametrization=change_to, calc_std_errors=calc_std_errors)
}


#' @title Estimate a G-StMAR model based on StMAR model with large degrees of freedom parameters
#'
#' @description \code{stmar_to_gstmar} estimates a G-StMAR model based on StMAR model with large degrees
#'  of freedom parameterss
#'
#' @inheritParams add_data
#' @inheritParams stmarpars_to_gstmar
#' @param estimate set \code{TRUE} if the new model should be estimated with variable metric algorithm using the StMAR model
#'   parameters as the initial values. By default \code{TRUE} iff the model contains data.
#' @param calc_std_errors set \code{TRUE} if the approximate standard errors should be calculated.
#'  By default \code{TRUE} iff the model contains data.
#'@param maxit the maximum number of iterations for the variable metric algorithm. Ignored if \code{estimate == FALSE}.
#' @details If a StMAR model contains large estimates for the degrees of freedom parameters
#'   one should consider switching to the corresponding G-StMAR model that lets the corresponding regimes to be GMAR type.
#'   \code{stmar_to_gstmar}  makes it convenient to do this switch.
#' @inherit GSMAR references return
#' @inherit add_data seealso
#' @examples
#' \donttest{
#'  # These are long running examples and use parallel computing
#'  fit13tr <- fitGSMAR(logVIX, 1, 3, model="StMAR", restricted=TRUE)
#'  fit13tr
#'  fit13gsr <- stmar_to_gstmar(fit13tr)
#'  fit13gsr
#' }
#' @export

stmar_to_gstmar <- function(gsmar, maxdf=100, estimate, calc_std_errors, maxit=100) {
  if(gsmar$model$model != "StMAR") stop("Only StMAR models are supported as input")
  if(missing(estimate)) estimate <- ifelse(is.null(gsmar$data), FALSE, TRUE)
  if(missing(calc_std_errors)) calc_std_errors <- ifelse(is.null(gsmar$data), FALSE, TRUE)
  if(estimate == TRUE & is.null(gsmar$data)) stop("Can't estimate the model without data")

  new_params <- stmarpars_to_gstmar(p=gsmar$model$p, M=gsmar$model$M, params=gsmar$params,
                                    restricted=gsmar$model$restricted, constraints=gsmar$model$constraints,
                                    maxdf=maxdf)
  if(is.null(gsmar$model$constraints) || gsmar$model$restricted == TRUE) {
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

  if(estimate == TRUE) {
    tmp_mod <- GSMAR(data=gsmar$data, p=gsmar$model$p, M=new_M, params=new_params$params,
                     model=new_model, restricted=gsmar$model$restricted,
                     constraints=new_constraints, conditional=gsmar$model$conditional,
                     parametrization=gsmar$model$parametrization, calc_qresiduals=FALSE,
                     calc_cond_moments=FALSE, calc_std_errors=FALSE)
    new_mod <- iterate_more(tmp_mod, maxit=maxit)
  } else {
    new_mod <- GSMAR(data=gsmar$data, p=gsmar$model$p, M=new_M, params=new_params$params,
                     model=new_model, restricted=gsmar$model$restricted,
                     constraints=new_constraints, conditional=gsmar$model$conditional,
                     parametrization=gsmar$model$parametrization, calc_std_errors=calc_std_errors)
  }
  new_mod
}


#' @title Construct a GSMAR model based on results from an arbitrary estimation round of \code{fitGSMAR}
#'
#' @description \code{alt_gsmar} constructs a GSMAR model based on results from an arbitrary estimation round of \code{fitGSMAR}.
#'
#' @inheritParams simulateGSMAR
#' @param which_round based on which estimation round should the model be constructed? An integer value in 1,...,\code{ncalls}.
#' @details It's sometimes useful to examine other estimates than the one with the highest log-likelihood value. This function
#'   is just a simple wrapper to \code{GSMAR} that picks the correct estimates from an returned by \code{fitGSMAR}.
#' @inherit GSMAR references return
#' @inherit add_data seealso
#' @examples
#' \donttest{
#'  # These are long running examples and use parallel computing
#'  fit12t <- fitGSMAR(IE, 1, 2, model="StMAR", ncalls=2, seeds=1:2)
#'  fit12t
#'  fit12t2 <- alt_gsmar(fit12t, which_round=2)
#'  fit12t2
#' }
#' @export

alt_gsmar <- function(gsmar, which_round=1) {
  stopifnot(!is.null(gsmar$all_estimates))
  stopifnot(which_round >= 1 || which_round <= length(gsmar$all_estimates))
  GSMAR(data=gsmar$data, p=gsmar$model$p, M=gsmar$model$M, params=gsmar$all_estimates[[which_round]],
        model=gsmar$model$model, restricted=gsmar$model$restricted, constraints=gsmar$model$constraints,
        conditional=gsmar$model$conditional, parametrization=gsmar$model$parametrization,
        calc_qresiduals=TRUE, calc_cond_moments=TRUE, calc_std_errors=TRUE)
}

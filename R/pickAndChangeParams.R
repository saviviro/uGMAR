#' @title Pick phi0 or mean parameters from parameter vector
#'
#' @description \code{pick_phi0} picks and returns the phi0 or mean parameters from parameter vector.
#'
#' @inheritParams loglikelihood_int
#' @return Returns a vector of length \code{M} containing the phi0 or mean parameters depending
#'  parametrization.

pick_phi0 <- function(p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL) {
  if(!is.null(constraints)) {
    params <- reformConstrainedPars(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints)
  }
  M <- sum(M)
  if(restricted == FALSE) {
    return(matrix(params[1:(M*(p + 2))], ncol=M)[1,])
  } else {
    return(params[1:M])
  }
}


#' @title Pick degrees of freedom parameters from parameter vector
#'
#' @description \code{pick_dfs} picks and returns the degrees of freedom parameters from parameter vector.
#'
#' @inheritParams loglikelihood_int
#' @return Returns a vector of length \code{M} or \code{M2} containing the  degrees of freedom parameters

pick_dfs <- function(p, M, params, model=c("GMAR", "StMAR", "G-StMAR")) {
  if(model == "GMAR") {
    return(numeric(0))
  } else if(model == "G-StMAR") {
    M2 <- M[2]
  } else {
    M2 <- M
  }
  d <- length(params)
  params[(d - M2 + 1):d]
}


#' @title Pick mixing weights parameters from parameter vector
#'
#' @description \code{pick_alphas} picks and returns the mixing weights parameters
#'  (including the non-parametrized one for the last component) from parameter vector.
#'
#' @inheritParams loglikelihood_int
#' @return Returns a vector of length \code{M} containing the mixing weights parameters \eqn{\alpha_m}.

pick_alphas <- function(p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL) {
  model <- match.arg(model)
  if(sum(M) == 1) {
    return(1)
  } else {
    if(!is.null(constraints)) {
      params <- reformConstrainedPars(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints)
    }
    M <- sum(M)
    if(restricted == FALSE) {
      alphas <- params[(M*(p + 2) + 1):(M*(p + 3) - 1)]
    } else {
      alphas <- params[(p + 2*M + 1):(3*M + p - 1)]
    }
  }
  c(alphas, 1-sum(alphas))
}


#' @title Pick \eqn{\phi_0}/\eqn{\mu}, AR-coefficients and variance parameters from parameter vector
#'
#' @description \code{pick_pars} picks \eqn{\phi_0}/\eqn{\mu}, ar-coefficient and variance parameters from parameter vector
#'
#' @inheritParams loglikelihood_int
#' @return Returns a \eqn{(Mx(p+2))} matrix containing the parameters, column for each component. First row for \eqn{\phi_0}/\eqn{\mu} depending on the parametrization,
#'   second row for \eqn{\phi_1}, second last row for \eqn{\phi_p} and last row for \eqn{\sigma^2}.

pick_pars <- function(p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL) {
  params <- removeAllConstraints(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints)
  matrix(params[1:(sum(M)*(p + 2))], ncol=sum(M))
}


#' @title Change parametrization of the parameter vector
#'
#' @description \code{change_parametrization} changes the parametrization of the given parameter
#'   vector to \code{change_to}.
#'
#' @inheritParams loglikelihood_int
#' @param change_to either "intercept" or "mean" specifying to which parametrization it should be switched to.
#'   If set to \code{"intercept"}, it's assumed that \code{params} is mean-parametrized, and if set to \code{"mean"}
#'   it's assumed that \code{params} is intercept-parametrized.
#' @return Returns parameter vector described in \code{params}, but with parametrization changed from intercept to mean
#'   (when \code{change_to==mean}) or from mean to intercept (when \code{change_to==intercept}).
#' @section Warning:
#'  No argument checks!
#' @inherit isStationary references

change_parametrization <- function(p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE,
                                   constraints=NULL, change_to=c("intercept", "mean")) {
  model <- match.arg(model)
  change_to <- match.arg(change_to)
  stopifnot(change_to %in% c("intercept", "mean"))
  params_orig <- params
  params <- reformConstrainedPars(p=p, M=M, params=params, model=model, restricted=restricted,
                                  constraints=constraints)
  pars <- reformParameters(p=p, M=M, params=params, model=model, restricted=restricted)$pars
  M <- sum(M)

  if(change_to == "intercept") { # Current parametrization is "mean"
    new_pars <- vapply(1:M, function(m) pars[1, m]*(1 - sum(pars[(2:(nrow(pars) - 1)), m])), numeric(1))
  } else {  # change_to == "mean" and current parametrization is "intercept"
    new_pars <- vapply(1:M, function(m) pars[1, m]/(1 - sum(pars[(2:(nrow(pars) - 1)), m])), numeric(1))
  }

  if(restricted == FALSE) {
    if(is.null(constraints)) {
      all_q <- rep(p, times=M)
    } else {
      all_q <- vapply(1:M, function(m) ncol(as.matrix(constraints[[m]])), numeric(1))
    }
    j <- 0
    for(m in 1:M) {
      params_orig[j + 1] <- new_pars[m]
      j <- j + all_q[m] + 2
    }
  } else {
    params_orig[1:M] <- new_pars
  }
  params_orig
}


#' @title Calculate and return regime means \eqn{\mu_{m}}
#'
#' @description \code{get_regime_means} calculates regime means \eqn{\mu_{m} =  \phi_{m,0}/(1-\sum\phi_{i,m})}
#'   for the given GMAR, StMAR or G-StMAR model
#'
#' @inheritParams simulateGSMAR
#' @return Returns a length \code{M} vector containing regime mean \eqn{\mu_{m}} in the m:th column, \eqn{m=1,..,M}.
#' @inherit isStationary references
#' @examples
#' # GMAR model
#' params13 <- c(1.4, 0.88, 0.26, 2.46, 0.82, 0.74, 5.0, 0.68, 5.2, 0.72, 0.2)
#' gmar13 <- GSMAR(data=VIX, p=1, M=3, params=params13, model="GMAR")
#' get_regime_means(gmar13)
#'
#' # StMAR model
#' params12t <- c(1.38, 0.88, 0.27, 3.8, 0.74, 3.15, 0.8, 100, 3.6)
#' stmar12t <- GSMAR(data=VIX, p=1, M=2, params=params12t, model="StMAR")
#' get_regime_means(stmar12t)
#'
#' # G-StMAR model (similar to the StMAR model above)
#' params12gs <- c(1.38, 0.88, 0.27, 3.8, 0.74, 3.15, 0.8, 3.6)
#' gstmar12 <- GSMAR(data=VIX, p=1, M=c(1, 1), params=params12gs,
#'  model="G-StMAR")
#' get_regime_means(gstmar12)
#' @export

get_regime_means <- function(gsmar) {
  check_gsmar(gsmar)
  p <- gsmar$model$p
  M <- gsmar$model$M
  params <- gsmar$params
  model <- gsmar$model$model
  restricted <- gsmar$model$restricted
  constraints <- gsmar$model$constraints

  if(gsmar$model$parametrization == "intercept") {
    params <- change_parametrization(p=p, M=M, params=params, model=model, restricted=restricted,
                                     constraints=constraints, change_to="mean")
  }

  pick_phi0(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints)
}


#' @title Calculate absolute values of the roots of the AR characteristic polynomials
#'
#' @description \code{get_ar_roots} calculates absolute values of the roots of the AR characteristic polynomials
#'   for each component.
#'
#' @inheritParams simulateGSMAR
#' @return Returns a list with \code{M} elements each containing the absolute values of the roots
#'  of the AR characteristic polynomial corresponding to each mixture component.
#' @inherit isStationary references
#' @examples
#' params13 <- c(1.4, 0.88, 0.26, 2.46, 0.82, 0.74, 5.0, 0.68, 5.2, 0.72, 0.2)
#' gmar13 <- GSMAR(data=VIX, p=1, M=3, params=params13, model="GMAR")
#' get_ar_roots(gmar13)
#' @export

get_ar_roots <- function(gsmar) {
  check_gsmar(gsmar)
  p <- gsmar$model$p
  M <- gsmar$model$M
  params <- removeAllConstraints(p=p, M=M, params=gsmar$params, model=gsmar$model$model,
                                 restricted=gsmar$model$restricted, constraints=gsmar$model$constraints)
  M <- sum(M)
  pars <- matrix(params[1:(M*(p + 2))], ncol=M)
  lapply(1:M, function(i1) abs(polyroot(c(1, -pars[2:(p + 1), i1]))))
}



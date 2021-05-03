#' @title Transform a StMAR or G-StMAR model parameter vector to a corresponding G-StMAR model parameter vector
#'  with large dfs parameters reduced.
#'
#' @description \code{stmarpars_to_gstmar} transforms a StMAR model parameter vector to a corresponding
#'  G-StMAR model parameter vector with large dfs parameters reduced by switching the related regimes
#'  to be GMAR type.
#'
#' @inheritParams loglikelihood_int
#' @param maxdf regimes with degrees of freedom parameter value larger than this will be turned into
#'  GMAR type.
#' @return Returns a list with three elements: \code{$params} contains the corresponding G-StMAR model
#'  parameter vector, \code{$reg_order} contains the permutation that was applied to the regimes
#'  (GMAR type components first, and decreasing ordering by mixing weight parameters), and
#'  \code{$M} a vector of length two containing the number of GMAR type regimes in the first element
#'  and the number of StMAR type regimes in the second.
#' @examples
#'  params12 <- c(2, 0.9, 0.1, 0.8, 0.5, 0.5, 0.4, 12, 300)
#'  stmarpars_to_gstmar(p=1, M=2, params=params12, model="StMAR", maxdf=100)
#' @export

stmarpars_to_gstmar <- function(p, M, params, model=c("GMAR", "StMAR", "G-StMAR"),
                                restricted=FALSE, constraints=NULL, maxdf=100) {
  if(!model %in% c("StMAR", "G-StMAR")) stop("Only StMAR and G-StMAR models are supported!")
  check_pM(p, M, model=model)
  check_params_length(p, M, params, model=model, restricted=restricted, constraints=constraints)
  check_constraint_mat(p, M, restricted=restricted, constraints=constraints)
  M_orig <- M # Length two vector for G-StMAR model
  M <- sum(M) # The number of regimes
  if(model == "StMAR") {
    M1 <- 0
    M2 <- M
  } else { # model == "G-StMAR"
    M1 <- M_orig[1]
    M2 <- M_orig[2]
  }

  dfs <- pick_dfs(p=p, M=M_orig, params=params, model=model)
  if(!any(dfs > maxdf)) {
    warning("No degrees of freedom parameter is larger than 'maxdf'. The original model is returned.")
    if(model == "StMAR") {
      ret_M <- c(0, M)
    } else {
      ret_M <- M_orig
    }
    return(list(params=params,
                reg_order=1:M,
                M=ret_M))
  }
  regs_to_change <- which(dfs > maxdf) + M1
  if(length(regs_to_change) == M2) message("All regimes are changed to GMAR type. The result is therefore a GMAR model and not a G-StMAR model.")
  alphas <- pick_alphas(p=p, M=M_orig, params=params, model=model, restricted=restricted, constraints=constraints)
  all_regs <- lapply(1:M, function(i1) extract_regime(p=p, M=M_orig, params, model=model, restricted=restricted,
                                                      constraints=constraints, regime=i1, with_dfs=FALSE))
  if(model == "StMAR") {
    gmar_regs <- regs_to_change
  } else { # model == "G-StMAR
    gmar_regs <- c(1:M1, regs_to_change)
  }
  reg_order <- c(gmar_regs[order(alphas[gmar_regs], decreasing=TRUE)], # GMAR type regimes
                 (1:M)[-gmar_regs][order(alphas[-gmar_regs], decreasing=TRUE)]) # StMAR type regimes
  tmp_pars <- unlist(all_regs[reg_order])

  if(!is.null(constraints) && any(reg_order != 1:M)) {
    message(paste0("The order of the constraint matrices for was changed to ", toString(reg_order), "."))
  }
  if(restricted) { # Add the missing AR parameters
    q <- ifelse(is.null(constraints), p, ncol(constraints))
    tmp_pars <- c(tmp_pars[seq(from=1, to=2*M, by=2)],
                  params[(M + 1):(M + q)],
                  tmp_pars[seq(from=2, to=2*M, by=2)])
  }
  if(length(regs_to_change) == M2) {
    new_dfs <- numeric(0)
  } else {
    new_dfs <- dfs[-(regs_to_change - M1)][order(reg_order[(length(regs_to_change) + 1 + M1):M], decreasing=FALSE)]
  }
  list(params=c(tmp_pars, alphas[reg_order][-M], new_dfs),
       reg_order=reg_order,
       M=c(length(regs_to_change) + M1, M2 - length(regs_to_change)))
}


#' @title Pick phi0 or mean parameters from parameter vector
#'
#' @description \code{pick_phi0} picks and returns the phi0 or mean parameters from parameter vector.
#'
#' @inheritParams loglikelihood_int
#' @return Returns a vector of length \code{M} containing the phi0 or mean parameters depending
#'  parametrization.

pick_phi0 <- function(p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL) {
  if(!is.null(constraints)) {
    params <- reform_constrained_pars(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints)
  }
  M <- sum(M)
  if(restricted == FALSE) {
    return(matrix(params[1:(M*(p + 2))], ncol=M)[1,])
  } else {
    return(params[1:M])
  }
}


#' @title Pick degrees of freedom parameters from a parameter vector
#'
#' @description \code{pick_dfs} picks and returns the degrees of freedom parameters from
#'   the given parameter vector.
#'
#' @inheritParams loglikelihood_int
#' @return Returns a vector of length \code{M} or \code{M2} containing the degrees of freedom parameters.

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
#'   (including the non-parametrized one for the last component) from the given
#'   parameter vector.
#'
#' @inheritParams loglikelihood_int
#' @return Returns a vector of length \code{M} containing the mixing weight parameters \eqn{\alpha_m}.

pick_alphas <- function(p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL) {
  model <- match.arg(model)
  if(sum(M) == 1) {
    return(1)
  } else {
    if(!is.null(constraints)) {
      params <- reform_constrained_pars(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints)
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


#' @title Pick \eqn{\phi_0} (or \eqn{\mu}), AR-coefficients, and variance parameters from a parameter vector
#'
#' @description \code{pick_pars} picks \eqn{\phi_0}/\eqn{\mu}, AR-coefficients, and variance parameters from
#'  the given parameter vector.
#'
#' @inheritParams loglikelihood_int
#' @return Returns a \eqn{(Mx(p+2))} matrix containing the parameters, column for each component.
#'  The first row for \eqn{\phi_0} or \eqn{\mu} depending on the parametrization, the second row
#'  for \eqn{\phi_1}, ..., the second to last row for \eqn{\phi_p}, and the last row for \eqn{\sigma^2}.

pick_pars <- function(p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL) {
  params <- remove_all_constraints(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints)
  matrix(params[1:(sum(M)*(p + 2))], ncol=sum(M))
}


#' @title Change parametrization of a parameter vector
#'
#' @description \code{change_parametrization} changes the parametrization of the given parameter
#'   vector to \code{change_to}.
#'
#' @inheritParams loglikelihood_int
#' @param change_to either "intercept" or "mean" specifying to which parametrization it should be switched to.
#'   If set to \code{"intercept"}, it's assumed that \code{params} is mean-parametrized, and if set to \code{"mean"}
#'   it's assumed that \code{params} is intercept-parametrized.
#' @return Returns parameter vector described in \code{params} but with parametrization changed from intercept to mean
#'   (when \code{change_to==mean}) or from mean to intercept (when \code{change_to==intercept}).
#' @section Warning:
#'  No argument checks!
#' @inherit is_stationary references

change_parametrization <- function(p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE,
                                   constraints=NULL, change_to=c("intercept", "mean")) {
  model <- match.arg(model)
  change_to <- match.arg(change_to)
  stopifnot(change_to == "intercept" | change_to == "mean")
  params_orig <- params
  params <- reform_constrained_pars(p=p, M=M, params=params, model=model, restricted=restricted,
                                  constraints=constraints)
  pars <- reform_parameters(p=p, M=M, params=params, model=model, restricted=restricted)$pars
  M <- sum(M)

  if(change_to == "intercept") { # Current parametrization is "mean"
    new_pars <- pars[1, ]*(1 - colSums(pars[2:(p + 1), , drop=FALSE]))
  } else {  # change_to == "mean" and current parametrization is "intercept"
    new_pars <- pars[1, ]/(1 - colSums(pars[2:(p + 1), , drop=FALSE]))
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



#' @title Calculate absolute values of the roots of the AR characteristic polynomials
#'
#' @description \code{get_ar_roots} calculates the absolute values of the roots of the AR
#'   characteristic polynomials for each mixture component.
#'
#' @inheritParams simulateGSMAR
#' @return Returns a list with \code{M} elements each containing the absolute values of the roots
#'  of the AR characteristic polynomial corresponding to each mixture component.
#' @inherit is_stationary references
#' @examples
#' params12 <- c(1.70, 0.85, 0.30, 4.12, 0.73, 1.98, 0.63)
#' gmar12 <- GSMAR(data=simudata, p=1, M=2, params=params12, model="GMAR")
#' get_ar_roots(gmar12)
#' @export

get_ar_roots <- function(gsmar) {
  check_gsmar(gsmar)
  p <- gsmar$model$p
  M <- gsmar$model$M
  params <- remove_all_constraints(p=p, M=M, params=gsmar$params, model=gsmar$model$model,
                                 restricted=gsmar$model$restricted, constraints=gsmar$model$constraints)
  M <- sum(M)
  pars <- matrix(params[1:(M*(p + 2))], ncol=M)
  lapply(1:M, function(i1) abs(polyroot(c(1, -pars[2:(p + 1), i1]))))
}



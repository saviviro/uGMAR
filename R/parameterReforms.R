#' @title Reform any parameter vector into standard form.
#'
#' @description \code{reform_parameters} takes a parameter vector of any (non-constrained) GMAR, StMAR, or G-StMAR model
#'  and returns a list with the parameter vector in the standard form, parameter matrix containing AR coefficients and
#'  component variances, mixing weights alphas, and in case of StMAR or G-StMAR model also degrees of freedom parameters.
#'
#' @inheritParams is_stationary_int
#' @details This function does not support models imposing linear constraints. No argument checks in this function.
#' @return Returns a list with...
#' \describe{
#'   \item{\code{$params}}{parameter vector in the standard form.}
#'   \item{\code{$pars}}{corresponding parameter matrix containing AR coefficients and
#'     component variances. First row for phi0 or means depending on the parametrization. Column for each component.}
#'   \item{\code{$alphas}}{numeric vector containing mixing weight parameters for all of the components (also for the last one).}
#'   \item{\code{$dfs}}{numeric vector containing degrees of freedom parameters for all of components.
#'     Returned only if \code{model == "StMAR"} or \code{model == "G-StMAR"}.}
#'  }
#'  @keywords internal

reform_parameters <- function(p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE) {
  model <- match.arg(model)
  params <- reform_restricted_pars(p=p, M=M, params=params, model=model, restricted=restricted) # Parameters on the non-restricted form
  list(params=params,
       pars=pick_pars(p=p, M=M, params=params, model=model, restricted=FALSE, constraints=NULL),
       alphas=pick_alphas(p=p, M=M, params=params, model=model, restricted=FALSE, constraints=NULL),
       dfs=pick_dfs(p=p, M=M, params=params, model=model))
}


#' @title Reform parameter vector with linear constraints to correspond non-constrained parameter vector.
#'
#' @description \code{reform_constrained_pars} reforms the parameter vector of a model with linear constrains
#'   to the "standard form" so that it's comparable with non-constrained models.
#'
#' @inheritParams loglikelihood_int
#' @return Returns such parameter vector corresponding to the input vector that is the form described in \code{params}
#' for non-restricted or restricted models (for non-constrained models), and can hence be used just as the
#' parameter vectors of non-constrained models.
#' @keywords internal

reform_constrained_pars <- function(p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL) {
  if(is.null(constraints)) {
    return(params)
  } else {
    model <- match.arg(model)
    M_orig <- M
    M <- sum(M)
    if(restricted == FALSE) {
      params0 <- numeric(0)
      j <- 0 # Controls where we at in the parameter vector
      for(i1 in 1:M) { # Go through the regimes
        C_m <- as.matrix(constraints[[i1]])
        q_m <- ncol(C_m) # Number of AR parameters in the regime
        psi_m <- params[(j + 2):(j + q_m + 1)]
        params0 <- c(params0, params[j + 1], C_m%*%psi_m, params[j + q_m + 2]) # Expand the constraints to the AR coefficients
        j <- j + q_m + 2 # Update the counter
      }
      if(M > 1) {
        params0 <- c(params0, params[(j + 1):(j + M - 1)]) # add alphas
      }
    } else { # If restricted==TRUE
      q <- ncol(as.matrix(constraints)) # Number of AR parameters
      psi <- params[(M + 1):(M + q)]
      params0 <- c(params[1:M], constraints%*%psi, params[(M + q + 1):(2*M + q)]) # Expand the constraints to the AR coefficients
      if(M > 1) {
        params0 <- c(params0, params[(2*M + q + 1):(3*M + q - 1)]) # add alphas
      }
    }
  }
  c(params0, pick_dfs(p=p, M=M_orig, params=params, model=model)) # add dfs
}


#' @title Reform parameter vector with restricted autoregressive parameters to correspond non-restricted parameter vector.
#'
#' @description \code{reform_restricted_pars} reforms parameter vector with restricted autoregressive parameters to correspond
#'  non-restricted parameter vector.
#'
#' @inheritParams loglikelihood_int
#' @return Returns such parameter vector corresponding to the input vector that is the form described in \code{params}
#' for non-restricted models (for non-constrained models). Linear constraints are not supported.
#' @keywords internal

reform_restricted_pars <- function(p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE) {
  if(restricted == FALSE) {
    return(params)
  }
  model <- match.arg(model)

  # Pick the parameter values
  M_orig <- M
  M <- sum(M)
  phi0 <- params[1:M]
  arcoefs <- matrix(rep(params[(M + 1):(M + p)], M), ncol=M)
  sigmas <- params[(M + p + 1):(p + 2*M)]
  pars <- rbind(phi0, arcoefs, sigmas)
  alphas <- params[(p + 2*M + 1):(3*M + p - 1)]
  dfs <- pick_dfs(p=p, M=M_orig, params=params, model=model)

  # Collect parameters together in the "standard" form
  c(as.vector(pars), alphas, dfs)
}


#' @title Transform constrained and restricted parameter vector into the regular form
#'
#' @description \code{remove_all_constraints} transforms constrained and restricted parameter vector into the regular form.
#'
#' @inheritParams loglikelihood_int
#' @return Returns such parameter vector corresponding to the input vector that is the form described in \code{params}
#' for non-restricted and non-constrained models.
#' @keywords internal

remove_all_constraints <- function(p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL) {
  model <- match.arg(model)
  params <- reform_constrained_pars(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints) # Remove constraints
  reform_restricted_pars(p=p, M=M, params=params, model=model, restricted=restricted) # Remove "restricted" -type constraints
}



#' @title Sort the mixture components of a GMAR, StMAR, or G-StMAR model
#'
#' @description \code{sort_components} sorts mixture components of the specified GMAR, StMAR, or G-StMAR model
#'   according to the mixing weight parameters when the parameter vector has the "standard/regular form" for
#'   restricted or non-restricted models.
#'
#' @inheritParams is_stationary
#' @details This function does not support models imposing linear constraints.
#' @return Returns a parameter vector sorted according to its mixing weight parameters,
#'   described in \code{params}.
#' @keywords internal

sort_components <- function(p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE) {
  model <- match.arg(model)
  M_orig <- M
  M <- sum(M)
  if(M == 1) return(params)
  alphas <- pick_alphas(p=p, M=M_orig, params=params, model=model, restricted=restricted, constraints=NULL)
  dfs <- pick_dfs(p=p, M=M_orig, params=params, model=model)

  # Obtain the new ordering of the regimes
  if(model != "G-StMAR") {
    M2 <- ifelse(model == "StMAR", M, 0)
    ord <- order(alphas, decreasing=TRUE) # The new ordering of the regimes
    if(all(ord == 1:M)) return(params) # Already in the correct order
    if(model == "StMAR") dfs <- dfs[ord]
  } else { # G-StMAR model, sort the M1 GMAR components and M2 StMAR components eparately
    M1 <- M_orig[1]
    M2 <- M_orig[2]
    if(M1 == 1 && M2 == 1) return(params) # Only one component of each type - nothing to sort
    ord_M1 <- order(alphas[1:M1], decreasing=TRUE)
    ord_M2 <- order(alphas[(M1 + 1):M], decreasing=TRUE)
    ord <- c(ord_M1, M1 + ord_M2) # Overall ordering
    if(all(ord == 1:M)) return(params) # Already in the correct order
    dfs <- dfs[ord_M2]
  }

  # Sort the regimes
  alphas <- alphas[ord][-M]
  pars <- pick_pars(p=p, M=M_orig, params=params, model=model, restricted=restricted, constraints=NULL)
  pars <- pars[,ord]

  if(restricted == FALSE) {
    return(c(as.vector(pars), alphas, dfs))
  } else { # If restricted == TRUE
    return(c(pars[1,], pars[2:(p + 1), 1], pars[p + 2,], alphas, dfs))
  }
}



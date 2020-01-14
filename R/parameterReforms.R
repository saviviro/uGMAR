#' @title Reform any parameter vector into standard form.
#'
#' @description \code{reformParameters} takes a parameter vector of any (non-constrained) GMAR, StMAR, or G-StMAR model
#'  and returns a list with the parameter vector in the standard form, parameter matrix containing AR coefficients and
#'  component variances, mixing weights alphas, and in case of StMAR or G-StMAR model also degrees of freedom parameters.
#'
#' @inheritParams isStationary_int
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

reformParameters <- function(p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE) {
  model <- match.arg(model)
  params <- reformRestrictedPars(p=p, M=M, params=params, model=model, restricted=restricted)
  list(params=params,
       pars=pick_pars(p=p, M=M, params=params, model=model, restricted=FALSE, constraints=NULL),
       alphas=pick_alphas(p=p, M=M, params=params, model=model, restricted=FALSE, constraints=NULL),
       dfs=pick_dfs(p=p, M=M, params=params, model=model))
}


#' @title Reform parameter vector with linear constraints to correspond non-constrained parameter vector.
#'
#' @description \code{reformConstrainedPars} reforms the parameter vector of a model with linear constrains
#'   to the "standard form" so that it's comparable with non-constrained models.
#'
#' @inheritParams loglikelihood_int
#' @return Returns such parameter vector corresponding to the input vector that is the form described in \code{params}
#' for non-restricted or restricted models (for non-constrained models), and can hence be used just as the
#' parameter vectors of non-constrained models.

reformConstrainedPars <- function(p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL) {
  if(is.null(constraints)) {
    return(params)
  } else {
    model <- match.arg(model)
    M_orig <- M
    M <- sum(M)
    if(restricted == FALSE) {
      params0 <- numeric(0)
      j <- 0
      for(i1 in 1:M) {
        C_m <- as.matrix(constraints[[i1]])
        q_m <- ncol(C_m)
        psi_m <- params[(j + 2):(j + q_m + 1)]
        params0 <- c(params0, params[j + 1], C_m%*%psi_m, params[j + q_m + 2])
        j <- j + q_m + 2
      }
      if(M > 1) {
        params0 <- c(params0, params[(j + 1):(j + M - 1)]) # add alphas
      }
    } else { # If restricted==TRUE
      q <- ncol(as.matrix(constraints))
      psi <- params[(M + 1):(M + q)]
      params0 <- c(params[1:M], constraints%*%psi, params[(M + q + 1):(2*M + q)])
      if(M > 1) {
        params0 <- c(params0, params[(2*M + q + 1):(3*M + q - 1)]) # add alphas
      }
    }
  }
  c(params0, pick_dfs(p=p, M=M_orig, params=params, model=model)) # add dfs
}


#' @title Reform parameter vector with restricted autoregressive parameters to correspond non-restricted parameter vector.
#'
#' @description \code{reformRestrictedPars} reforms parameter vector with restricted autoregressive parameters to correspond
#'  non-restricted parameter vector.
#'
#' @inheritParams loglikelihood_int
#' @return Returns such parameter vector corresponding to the input vector that is the form describted in \code{params}
#' for non-restricted models (for non-constrained models). Linear constraints are not supported.

reformRestrictedPars <- function(p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE) {
  if(restricted == FALSE) {
    return(params)
  } else {
    model <- match.arg(model)
    M_orig <- M
    M <- sum(M)
    phi0 <- params[1:M]
    arcoefs <- matrix(rep(params[(M + 1):(M + p)], M), ncol=M)
    sigmas <- params[(M + p + 1):(p + 2*M)]
    pars <- rbind(phi0, arcoefs, sigmas)
    alphas <- params[(p + 2*M + 1):(3*M + p - 1)]
    dfs <- pick_dfs(p=p, M=M_orig, params=params, model=model)
    return(c(as.vector(pars), alphas, dfs))
  }
}


#' @title Transform constrained and restricted parameter vector into the regular form
#'
#' @description \code{removeAllConstraints} transforms constrained and restricted parameter vector into the regular form.
#'
#' @inheritParams loglikelihood_int
#' @return Returns such parameter vector corresponding to the input vector that is the form described in \code{params}
#' for non-restricted and non-constrained models.

removeAllConstraints <- function(p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL) {
  model <- match.arg(model)
  params <- reformConstrainedPars(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints)
  reformRestrictedPars(p=p, M=M, params=params, model=model, restricted=restricted)
}



#' @title Sort the mixture components of a GMAR, StMAR, or G-StMAR model
#'
#' @description \code{sortComponents} sorts mixture components of the specified GMAR, StMAR, or G-StMAR model
#'   according to the mixing weight parameters when the parameter vector has the "standard/regular form" for
#'   restricted or non-restricted models.
#'
#' @inheritParams isStationary
#' @details This function does not support models imposing linear constraints.
#' @return Returns a parameter vector sorted according to its mixing weight parameters,
#'   described in \code{params}.

sortComponents <- function(p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE) {
  model <- match.arg(model)
  if(model != "G-StMAR") {
    if(M == 1) {
      return(params)
    }
    alphas <- pick_alphas(p=p, M=M, params=params, model=model, restricted=restricted, constraints=NULL)
    sortedByAlphas <- order(alphas, decreasing=TRUE)
    alphas <- alphas[sortedByAlphas]
    pars <- pick_pars(p=p, M=M, params=params, model=model, restricted=restricted, constraints=NULL)
    pars <- pars[,sortedByAlphas]
    dfs <- pick_dfs(p=p, M=M, params=params, model=model)
    if(model == "StMAR") {
      dfs <- dfs[sortedByAlphas]
    }

    if(restricted == FALSE) {
      return(c(as.vector(pars), alphas[-M], dfs))
    } else { # If restricted == TRUE
      return(c(pars[1,], pars[2:(p + 1), 1], pars[p + 2,], alphas[-M], dfs))
    }
  } else { # If model == "G-StMAR": Sort the M1 and M2 components separately (parameter picks are non-standard)
    M1 <- M[1]
    M2 <- M[2]
    Msum <- sum(M)
    if(restricted == FALSE) {
      pars0 <- numeric(0) # Collect pars in here
      alphas0 <- numeric(0) # Collect alphas in here
      for(i1 in 1:2) { # Go through GMAR and StMAR parts
        if(i1 == 1) {
          pars <- matrix(params[1:(M1*(p + 2))], ncol=M1)
        } else {
          pars <- matrix(params[((M1*(p + 2)) + 1):(Msum*(p + 2))], ncol=M2)
        }
        if(M[i1] > 1) { # If only one component of a given type, no need to sort
          if(i1 == 1) {
            alphas <- params[(Msum*(p + 2) + 1):(Msum*(p + 2) + M1)]
          } else {
            alphas <- params[(Msum*(p + 2) + M1 + 1):(Msum*(p + 3) - 1)]
            dfs <- pick_dfs(p=p, M=M, params=params, model=model)
          }
          if(i1 == 2) { # If i1 == 2, add the non-parametrized alpha
            alphas <- c(alphas, 1 - (sum(alphas) + sum(alphas0)))
          }
          sortedByAlphas <- order(alphas, decreasing=TRUE)
          pars <- pars[,sortedByAlphas]
          alphas <- alphas[sortedByAlphas]
          if(i1 == 2) {
            dfs <- dfs[sortedByAlphas]
            alphas <- alphas[-M2] # Delete the last alpha
          }
        } else { # If only one component, still need to collect alphas from GMAR and dfs for StMAR
          if(i1 == 1) {
            alphas <- params[Msum*(p + 2) + 1]
          } else {
            dfs <- params[Msum*(p + 3)]
            alphas <- numeric(0) # No alphas for StMAR if only one StMAR component
          }
        }
        pars0 <- cbind(pars0, pars)
        alphas0 <- c(alphas0, alphas)
      }
      return(c(as.vector(pars0), alphas0, dfs))
    } else { # If restricted == TRUE & model == "G-StMAR"
      phi00 <- numeric(0)
      sigmas0 <- numeric(0)
      alphas0 <- numeric(0)
      arcoefs <- params[(Msum + 1):(Msum + p)]
      for(i1 in 1:2) {
        if(i1 == 1) {
          phi0 <- params[1:M1]
          sigmas <- params[(Msum + p + 1):(Msum + p + M1)]
        } else {
          phi0 <- params[(M1 + 1):Msum]
          sigmas <- params[(Msum + p + M1 + 1):(2*Msum + p)]
          dfs <- params[(3*Msum + p):(3*Msum + M2 + p - 1)]
        }
        if(M[i1] > 1) { # if M[i1]==1, no need to sort
          if(i1 == 1) {
            alphas <- params[(2*Msum + p + 1):(2*Msum + p + M1)]
          } else {
            alphas <- params[(2*Msum + p + M1 + 1):(3*Msum + p - 1)]
          }
          if(i1 == 2) { # If i1 == 2, add the non-parametrized alpha
            alphas <- c(alphas, 1 - (sum(alphas) + sum(alphas0)))
          }
          sortedByAlphas <- order(alphas, decreasing=TRUE)
          phi0 <- phi0[sortedByAlphas]
          sigmas <- sigmas[sortedByAlphas]
          alphas <- alphas[sortedByAlphas]
          if(i1 == 2) {
            dfs <- dfs[sortedByAlphas]
            alphas <- alphas[-M2] # Delete the last alpha
          }
          # }
        } else { # If only one component, no sort but collect alphas
          if(i1 == 1) {
            alphas <- params[2*Msum + p + 1]
          } else {
            alphas <- NULL # No alphas if only one StMAR-component
          }
        }
        phi00 <- c(phi00, phi0)
        sigmas0 <- c(sigmas0, sigmas)
        alphas0 <- c(alphas0, alphas)
      }
      return(c(phi00, arcoefs, sigmas0, alphas0, dfs))
    }
  }
}



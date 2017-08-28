#' @title Reform any parameter vector into standard form.
#'
#' @description \code{reformParameters} takes a parameter vector of any GMAR or StMAR model and returns a list with the
#'  parameter vector in the standard form, parameter matrix containing AR coefficients and component
#'  variances, mixing weights alphas and in case of StMAR also degrees of freedom parameters.
#'
#' @inheritParams isStationary
#' @details This function does not support models parametrized with general linear constraints! Nor does it have any argument checks.
#' @return Returns a list with...
#' \describe{
#'   \item{\code{$params}}{parameter vector in the standard form.}
#'   \item{\code{$pars}}{corresponding parameter matrix containing AR coefficients and
#'     component variances. Column for each component.}
#'   \item{\code{$alphas}}{numeric vector containing mixing weights for all components (also for the last one).}
#'   \item{\code{$dfs}}{numeric vector containing degrees of freedom parameters for all components.
#'     Returned only if \code{StMAR==TRUE}.}
#'  }

reformParameters <- function(p, M, params, StMAR=FALSE, restricted=FALSE) {
  if(M==1) {
    pars = matrix(params[1:(M*(p+2))], ncol=M)
    alphas = c(1)
    if(StMAR==TRUE) {
      dfs = params[p+3]
    }
  } else {
    if(restricted==FALSE) {
      pars = matrix(params[1:(M*(p+2))], ncol=M) # Component parameters by column (except alphas and df in case of StMAR)
      alphas = params[(M*(p+2)+1):(M*(p+3)-1)]
      alphas = c(alphas, 1-sum(alphas)) # last alpha
      if(StMAR==TRUE) {
        dfs = params[(M*(p+3)):(M*(p+4)-1)] # degrees of freedom
      }
    } else {
      # If restricted TRUE: transform the restricted parameter vector into the standard form.
      phi0 = params[1:M]
      arcoefs = matrix(rep(params[(M+1):(M+p)], M), ncol=M)
      variances = params[(M+p+1):(p+2*M)]
      pars = rbind(phi0, arcoefs, variances)
      alphas = params[(p+2*M+1):(3*M+p-1)]
      if(StMAR==TRUE) {
        dfs = params[(3*M+p):(4*M+p-1)] # degrees of freedom
        params = c(as.vector(pars), alphas, dfs)
      } else {
        params = c(as.vector(pars), alphas)
      }
      alphas = c(alphas, 1-sum(alphas))
    }
  }
  rownames(pars) = NULL
  ret = list(params, pars, alphas)
  names(ret) = c("params", "pars", "alphas")
  if(StMAR==TRUE) {
    ret[[length(ret)+1]] = dfs
    names(ret)[length(ret)] = "dfs"
  }
  return(ret)
}


#' @title Reform parameter vector with linear constraints to correspond non-constrained parameter vector.
#'
#' @description \code{reformConstrainedParameters} reforms parameter vector of a model with linear constrains
#'   to the "standard form" so that it's comparable with non-constrained models.
#'
#' @inheritParams loglikelihood
#' @return Returns such parameter vector corresponding to the input vector that is the form describted in \code{params}
#' for non-restricted or restricted models (for non-constrained models), and can hence be used just as the
#' parameter vectors of non-constrained models.

reformConstrainedPars <- function(p, M, params, StMAR=FALSE, restricted=FALSE, R) {
  if(restricted==FALSE) {
    params0 = c()
    j = 0
    for(i1 in 1:M) {
      R_m = as.matrix(R[[i1]])
      q_m = ncol(R_m)
      psi_m = params[(j+2):(j+q_m+1)]
      params0 = c(params0, params[j+1], R_m%*%psi_m, params[j+q_m+2])
      j = j + q_m + 2
    }
    if(M>1) {
      params0 = c(params0, params[(j+1):(j+M-1)]) # add alphas
    }
    if(StMAR==TRUE) {
      params0 = c(params0, params[(j+M):(j+2*M-1)]) # add dfs
    }
  } else { # If restricted==TRUE
    q = ncol(as.matrix(R))
    psi = params[(M+1):(M+q)]
    params0 = c(params[1:M], R%*%psi, params[(M+q+1):(2*M+q)])
    if(M>1) {
      params0 = c(params0, params[(2*M+q+1):(3*M+q-1)]) # add alphas
    }
    if(StMAR==TRUE) {
      params0 = c(params0, params[(3*M+q):(4*M+q-1)]) # add dfs
    }
  }
  return(params0)
}



#' @title Sort the mixture components of GMAR or StMAR model
#'
#' @description \code{sortComponents} sorts mixture components of the specified GMAR or StMAR model by the mixing weights.
#'
#' @inheritParams isStationary
#' @details This function does not support models parametrized with general linear constraints!
#' @return Returns a parameter vector sorted by it's mixing weights.
#'  \describe{
#'    \item{For \strong{non-restricted} models:}{
#'      \describe{
#'        \item{For \strong{GMAR} model:}{Size \eqn{(M(p+3)-1x1)} vector \strong{\eqn{\theta}}\eqn{=}(\strong{\eqn{\upsilon_{1}}},...,\strong{\eqn{\upsilon_{M}}},
#'          \eqn{\alpha_{1},...,\alpha_{M-1}}), where \strong{\eqn{\upsilon_{m}}}\eqn{=(\phi_{m,0},}\strong{\eqn{\phi_{m}}}\eqn{,
#'          \sigma_{m}^2)} and \strong{\eqn{\phi_{m}}}=\eqn{(\phi_{m,1},...,\phi_{m,p}), m=1,...,M}.}
#'        \item{For \strong{StMAR} model:}{Size \eqn{(M(p+4)-1x1)} vector (\strong{\eqn{\theta, v}})\eqn{=}(\strong{\eqn{\upsilon_{1}}},...,\strong{\eqn{\upsilon_{M}}},
#'          \eqn{\alpha_{1},...,\alpha_{M-1}, v_{1},...,v_{M}}).}
#'      }
#'    }
#'    \item{For \strong{restricted} models:}{
#'      \describe{
#'        \item{For \strong{GMAR} model:}{Size \eqn{(3M+p-1x1)} vector \strong{\eqn{\theta}}\eqn{=(\phi_{1,0},...,\phi_{M,0},}\strong{\eqn{\phi}}\eqn{,
#'          \sigma_{1}^2,...,\sigma_{M}^2,\alpha_{1},...,\alpha_{M-1})}, where \strong{\eqn{\phi}}=\eqn{(\phi_{1},...,\phi_{M})}.}
#'        \item{For \strong{StMAR} model:}{Size \eqn{(4M+p-1x1)} vector (\strong{\eqn{\theta, v}})\eqn{=(\phi_{1,0},...,\phi_{M,0},}\strong{\eqn{\phi}}\eqn{,
#'          \sigma_{1}^2,...,\sigma_{M}^2,\alpha_{1},...,\alpha_{M-1}, v_{1},...,v_{M})}.}
#'      }
#'    }
#'  }
#' @section Warning:
#'  This function doesn't have any argument checks.

sortComponents <- function(p, M, params, StMAR=FALSE, restricted=FALSE) {
  if(M==1) {
    return(params)
  }
  if(restricted==FALSE) {
    pars = matrix(params[1:(M*(p+2))], ncol=M) # Component parameters by column (except alphas and dfs)
    alphas = params[(M*(p+2)+1):(M*(p+3)-1)]
    if(length(alphas)==1) {
      if(alphas<0.5) {
        pars = pars[,c(2,1)]
        if(StMAR==TRUE) {
          dfs = params[(M*(p+3)):(M*(p+4)-1)] # degrees of freedom
          return(c(as.vector(pars), 1-alphas, rev(dfs)))
        } else {
          return(c(as.vector(pars), 1-alphas))
        }
      } else {
        return(params)
      }
    } else {
      alphas0 = c(alphas, 1-sum(alphas))
      sortedByAlphas = order(alphas0, decreasing=TRUE)
      pars = pars[,sortedByAlphas]
      alphas = alphas0[sortedByAlphas]
      if(StMAR==TRUE) {
        dfs = params[(M*(p+3)):(M*(p+4)-1)] # degrees of freedom
        dfs = dfs[sortedByAlphas]
        return(c(as.vector(pars), alphas[-M], dfs))
      } else {
        return(c(as.vector(pars), alphas[-M]))
      }
    }
  } else { # If restricted==TRUE
    phi0 = params[1:M]
    arcoefs = params[(M+1):(M+p)]
    armat = matrix(rep(arcoefs, M), ncol=M)
    variances = params[(M+p+1):(p+2*M)]
    pars = rbind(phi0, armat, variances)
    alphas = params[(p+2*M+1):(3*M+p-1)]
    if(length(alphas)==1) {
      if(alphas<0.5) {
        pars = pars[,c(2,1)]
        if(StMAR==TRUE) {
          dfs = params[(3*M+p):(4*M+p-1)] # degrees of freedom
          pars = pars[,c(2,1)]
          return(c(pars[1,], arcoefs, pars[p+2,], 1-alphas, rev(dfs)))
        } else {
          return(c(pars[1,], arcoefs, pars[p+2,], 1-alphas))
        }
      } else {
        return(params)
      }
    } else {
      alphas0 = c(alphas, 1-sum(alphas))
      sortedByAlphas = order(alphas0, decreasing=TRUE)
      pars = pars[,sortedByAlphas]
      alphas = alphas0[sortedByAlphas]
      if(StMAR==TRUE) {
        dfs = params[(3*M+p):(4*M+p-1)] # degrees of freedom
        dfs = dfs[sortedByAlphas]
        return(c(pars[1,], arcoefs, pars[p+2,], alphas[-M], dfs))
      } else {
        return(c(pars[1,], arcoefs, pars[p+2,], alphas[-M]))
      }
    }
  }
}



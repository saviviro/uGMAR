#' @title Check the stationary and identification conditions of specified GMAR or StMAR model.
#'
#' @description FOR INTERNAL USE. \code{isStationary_int} checks the stationary condition and \code{isIdentifiable} checks the identification conditions
#'  of the specified GMAR or StMAR model.
#'
#' @inheritParams loglikelihood
#' @param params a real valued parameter vector specifying the model.
#'  \describe{
#'    \item{For \strong{non-restricted} models:}{
#'      \describe{
#'        \item{For \strong{GMAR} model:}{Size \eqn{(M(p+3)-1x1)} vector \strong{\eqn{\theta}}\eqn{=}(\strong{\eqn{\upsilon_{1}}},...,\strong{\eqn{\upsilon_{M}}},
#'          \eqn{\alpha_{1},...,\alpha_{M-1}}), where \strong{\eqn{\upsilon_{m}}}\eqn{=(\phi_{m,0},}\strong{\eqn{\phi_{m}}}\eqn{,
#'          \sigma_{m}^2)} and \strong{\eqn{\phi_{m}}}=\eqn{(\phi_{m,1},...,\phi_{m,p}), m=1,...,M}.}
#'        \item{For \strong{StMAR} model:}{Size \eqn{(M(p+4)-1x1)} vector (\strong{\eqn{\theta, \nu}})\eqn{=}(\strong{\eqn{\upsilon_{1}}},...,\strong{\eqn{\upsilon_{M}}},
#'          \eqn{\alpha_{1},...,\alpha_{M-1}, \nu_{1},...,\nu_{M}}).}
#'      }
#'    }
#'    \item{For \strong{restricted} models:}{
#'      \describe{
#'        \item{For \strong{GMAR} model:}{Size \eqn{(3M+p-1x1)} vector \strong{\eqn{\theta}}\eqn{=(\phi_{1,0},...,\phi_{M,0},}\strong{\eqn{\phi}}\eqn{,
#'          \sigma_{1}^2,...,\sigma_{M}^2,\alpha_{1},...,\alpha_{M-1})}, where \strong{\eqn{\phi}}=\eqn{(\phi_{1},...,\phi_{M})}.}
#'        \item{For \strong{StMAR} model:}{Size \eqn{(4M+p-1x1)} vector (\strong{\eqn{\theta, \nu}})\eqn{=(\phi_{1,0},...,\phi_{M,0},}\strong{\eqn{\phi}}\eqn{,
#'          \sigma_{1}^2,...,\sigma_{M}^2,\alpha_{1},...,\alpha_{M-1}, \nu_{1},...,\nu_{M})}.}
#'      }
#'    }
#'  }
#'  Symbol \eqn{\phi} denotes an AR coefficient, \eqn{\sigma^2} a variance, \eqn{\alpha} a mixing weight and \eqn{v} a degrees of
#'  freedom parameter.
#'  Note that in the case \strong{M=1} the parameter \eqn{\alpha} is dropped, and in the case of \strong{StMAR} model
#'  the degrees of freedom parameters \eqn{\nu_{m}} have to be larger than \eqn{2}.
#' @details These functions don't support models parametrized with general linear constraints.
#' @return Returns \code{TRUE} or \code{FALSE} accordingly.
#' @section Warning:
#'  These functions don't have any argument checks!
#' @references
#'  \itemize{
#'    \item Kalliovirta L., Meitz M. and Saikkonen P. (2015) Gaussian Mixture Autoregressive model for univariate time series.
#'          \emph{Journal of Time Series Analysis}, \strong{36}, 247-266.
#'    \item References regarding the StMAR and G-StMAR models will be updated when they are published.
#'  }

isStationary_int <- function(p, M, params, restricted=FALSE) {
  M = sum(M) # For G-StMAR models
  if(restricted==FALSE) {
    pars = matrix(params[1:(M*(p+2))], ncol=M)
    for(i1 in 1:M) {
      if(any(abs(polyroot(c(1, -pars[2:(p+1), i1])))<=1+1e-10)) {
        return(FALSE)
      }
    }
  } else {
    absroots = abs(polyroot(c(1, -params[(M+1):(M+p)])))
    if(any(absroots<=1+1e-10)) {
      return(FALSE)
    }
  }
  return(TRUE)
}


#' @title Check the stationary condition of specified GMAR, StMAR or G-StMAR model.
#'
#' @description \code{isStationary} checks the stationarity condition of the specified GMAR, StMAR or G-StMAR model.
#'
#' @inheritParams loglikelihood
#' @details This function uses numerical approximations and it will falsely return \code{FALSE} for stationary models
#'   when the stationarity condition is really close to break.
#' @inherit isStationary_int return references
#' @examples
#' # GMAR model
#' params22 <- c(0.4, 0.39, 0.6, 0.3, 0.4, 0.1, 0.6, 0.3, 0.8)
#' isStationary(2, 2, params22)
#'
#' # StMAR model
#' params12t <- c(-0.3, 1, 0.9, 0.1, 0.8, 0.6, 0.7, 10, 12)
#' isStationary(1, 2, params12t, StMAR=TRUE)
#'
#' # G-StMAR model
#' params12gs <- c(1, 0.1, 1, 2, 0.2, 2, 0.8, 20)
#' isStationary(1, c(1,1), params12gs, GStMAR=TRUE)
#'
#' # Restricted GMAR model
#' params13r <- c(0.1, 0.2, 0.3, -0.99, 0.1, 0.2, 0.3, 0.5, 0.3)
#' isStationary(1, 3, params13r, restricted=TRUE)
#'
#' # Restricted StMAR model
#' params22tr <- c(-0.1, -0.2, 0.01, 0.99, 0.3, 0.4, 0.9, 3, 13)
#' isStationary(2, 2, params22tr, StMAR=TRUE, restricted=TRUE)
#'
#' # Restricted G-StMAR model
#' params13gsr <- c(1, 2, 3, -0.99, 1, 2, 3, 0.5, 0.4, 20, 30)
#' isStationary(1, c(1,2), params13gsr, GStMAR=TRUE, restricted=TRUE)
#'
#' # GMAR model as a mixture of AR(2) and AR(1) models
#' R <- list(diag(1, ncol=2, nrow=2), as.matrix(c(1, 0)))
#' params22c <- c(1.2, 0.8, 0.2, 0.3, 3.3, 0.7, 3, 0.8)
#' isStationary(2, 2, params22c, constraints=TRUE, R=R)
#'
#' # Such StMAR(3,2) that the AR coefficients are restricted to be the
#' # same for both regimes and that the second AR coefficients are
#' # constrained to zero.
#' params32trc <- c(1, 2, 0.8, -0.3, 1, 2, 0.7, 11, 12)
#' isStationary(3, 2, params32trc, StMAR=TRUE, restricted=TRUE,
#'              constraints=TRUE, R=matrix(c(1, 0, 0, 0, 0, 1), ncol=2))
#' @export

isStationary <- function(p, M, params, StMAR=FALSE, GStMAR=FALSE, restricted=FALSE, constraints=FALSE, R) {
  checkLogicals(StMAR=StMAR, GStMAR=GStMAR)
  checkPM(p=p, M=M, GStMAR=GStMAR)
  if(length(params)!=nParams(p=p, M=M, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted, constraints=constraints, R=R)) {
    stop("The parameter vector has wrong dimension")
  }
  if(constraints==TRUE) {
    checkConstraintMat(p=p, M=M, R, restricted=restricted)
    params = reformConstrainedPars(p, M, params, R=R, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted)
  }
  return(isStationary_int(p=p, M=M, params=params, restricted=restricted))
}


#' @rdname isStationary_int

isIdentifiable <- function(p, M, params, restricted=FALSE, StMAR=FALSE, GStMAR=FALSE) {
  if(GStMAR==TRUE) {
    M1 = M[1]
    M2 = M[2]
    M = sum(M)
    if(M1==1 & M2==1) {
      return(TRUE)
    }
  }
  if(M==1) {
    return(TRUE)
  }
  if(restricted==FALSE) {
    pars = matrix(params[1:(M*(p+2))], ncol=M)
    alphas = params[(M*(p+2)+1):(M*(p+3)-1)]
    alphas = c(alphas, 1-sum(alphas))
    if(StMAR==TRUE) {
      pars = rbind(pars, params[(M*(p+3)):(M*(p+4)-1)])
    } else if(GStMAR==TRUE) {
      pars1 = as.matrix(pars[,1:M1])
      pars2 = as.matrix(pars[,(M1+1):M])
      pars2 = rbind(pars2, params[(M*(p+3)):(M*(p+3)+M2-1)])
      alphas1 = alphas[1:M1]
      alphas2 = alphas[(M1+1):M]
    }
  } else { # If restricted==TRUE
    pars = rbind(params[1:M], matrix(rep(params[(M+1):(M+p)], M), ncol=M), params[(M+p+1):(p+2*M)])
    alphas = params[(p+2*M+1):(3*M+p-1)]
    alphas = c(alphas, 1-sum(alphas))
    if(StMAR==TRUE) {
      pars = rbind(pars, params[(3*M+p):(4*M+p-1)])
    } else if(GStMAR==TRUE) {
      pars1 = as.matrix(pars[,1:M1])
      pars2 = as.matrix(pars[,(M1+1):M])
      pars2 = rbind(pars2, params[(3*M+p):(3*M+p+M2-1)])
      alphas1 = alphas[1:M1]
      alphas2 = alphas[(M1+1):M]
    }
  }

  if(GStMAR==TRUE) { # Check GMAR parameters and StMAR-parameters separately
    if(M1>1) {
      if(is.unsorted(rev(alphas1), strictly=TRUE)) {
        return(FALSE)
      } else if(anyDuplicated(t(pars1))!=0) {
        return(FALSE)
      }
    }
    if(M2>1) {
      if(is.unsorted(rev(alphas2), strictly=TRUE)) {
        return(FALSE)
      } else if(anyDuplicated(t(pars2))!=0) {
        return(FALSE)
      }
    }
  } else { # if GStMAR is NOT true
    if(is.unsorted(rev(alphas), strictly=TRUE)) {
      return(FALSE)
    } else if(anyDuplicated(t(pars))!=0) {
      return(FALSE)
    }
  }
  return(TRUE)
}


#' @title Check the data is set correctly and correct if not
#'
#' @description \code{checkAndCorrectData} checks that the data is set correctly and corrects it if not.
#'  Throws an error if it can't convert the data to the correct form.
#'
#' @inheritParams loglikelihood
#' @return Returns a numeric column matrix containing the data.

checkAndCorrectData <- function(data, p) {
  if(!is.matrix(data)) {
    data = as.matrix(data)
  }
  if(ncol(data)>nrow(data)) {
    data=t(data)   # The data matrix should contain observation per row (not per column)
  }
  if(ncol(data)>1) {
    stop("The data has more than one columns")
  } else if(anyNA(data)) {
    stop("The data contains NA values, which is not supported")
  } else if(length(data)<p+2) {
    stop("The data must contain at least p+2 values")
  }
  return(data)
}


#' @title Check the parameter vector is specified correctly
#'
#' @description \code{parameterChecks} checks dimension, restrictions, stationarity and identifibility of the given parameters
#'   of GMAR, StMAR or G-StMAR model. Throws an error if any check fails.
#'
#' @inheritParams loglikelihood
#' @param params a real valued parameter vector specifying the model.
#'      \describe{
#'        \item{For \strong{GMAR} model:}{Size \eqn{(M(p+3)-1x1)} vector \strong{\eqn{\theta}}\eqn{=}(\strong{\eqn{\upsilon_{1}}},...,\strong{\eqn{\upsilon_{M}}},
#'          \eqn{\alpha_{1},...,\alpha_{M-1}}), where \strong{\eqn{\upsilon_{m}}}\eqn{=(\phi_{m,0},}\strong{\eqn{\phi_{m}}}\eqn{,
#'          \sigma_{m}^2)} and \strong{\eqn{\phi_{m}}}=\eqn{(\phi_{m,1},...,\phi_{m,p}), m=1,...,M}.}
#'        \item{For \strong{StMAR} model:}{Size \eqn{(M(p+4)-1x1)} vector (\strong{\eqn{\theta, \nu}})\eqn{=}(\strong{\eqn{\upsilon_{1}}},...,\strong{\eqn{\upsilon_{M}}},
#'          \eqn{\alpha_{1},...,\alpha_{M-1}, \nu_{1},...,\nu_{M}}).}
#'        \item{For \strong{G-StMAR} model:}{Size \eqn{(M(p+3)+M2-1x1)} vector (\strong{\eqn{\theta, \nu}})\eqn{=}(\strong{\eqn{\upsilon_{1}}},...,\strong{\eqn{\upsilon_{M}}},
#'          \eqn{\alpha_{1},...,\alpha_{M-1}, \nu_{M1+1},...,\nu_{M}}).}
#'      }
#' @param pars a parameter matrix containing parameters \eqn{(\phi_{i,0},...,\phi_{i,p}, \sigma_{i}^2)}
#'  so that i:th column denotes i:th component.
#' @param alphas size Mx1 parameter vector containing mixing weights for all components.
#' @details Note that the "params" -parameter vector is assumed to be in the "standard" form for restricted models as well.
#' @return Throws an error if any check fails. Doesn't return anything.

parameterChecks <- function(p, M, params, pars, alphas, StMAR=FALSE, GStMAR=FALSE, constraints=FALSE) {
  if(StMAR==TRUE) {
    dfs = params[(M*(p+3)):(M*(p+4)-1)]
    if(length(params)!=(M*(p+4)-1)) {
      stop("The parameter vector has wrong dimension")
    } else if(any(dfs<=2)) {
      stop("The degrees of freedom parameters have to be larger than 2")
    }
  } else if(GStMAR==TRUE) {
    M1 = M[1]
    M2 = M[2]
    M = sum(M)
    dfs = params[(M*(p+3)):(M*(p+3)+M2-1)]
    if(length(params)!=M*(p+3)-1+M2) {
      stop("The parameter vector has wrong dimension")
    } else if(any(dfs<=2)) {
      stop("The degrees of freedom parameters have to be larger than 2")
    }
  } else {
    if(length(params)!=(M*(p+3)-1)) {
      stop("The parameter vector has wrong dimension")
    }
  }
  if(M>=2 & sum(alphas[-M])>=1) {
    stop("The mixing weights don't sum to one")
  } else if(any(pars[p+2,]<=0)) {
    stop("Component variances have to be larger than zero")
  } else if(!isStationary_int(p, M, params, restricted=FALSE)) {
    stop("The model doesn't satisfy the stationary condition")
  }
  if(GStMAR==TRUE) { M = c(M1, M2)}
  if(constraints==FALSE & !isIdentifiable(p, M, params, restricted=FALSE, StMAR=StMAR, GStMAR=GStMAR)) {
    stop("The model doesn't satisfy the identification conditions")
  }
}


#' @title Check constraint matrices R
#'
#' @description \code{checkConstraintMat} ckecks for some parts that the constraint matrices R are correctly set.
#' @inheritParams loglikelihood
#' @return Doesn't return anything, but throws an error if finds out that something is wrong.

checkConstraintMat <- function(p, M, R, restricted=FALSE) {
  M = sum(M) # For G-StMAR
  if(restricted==TRUE) {
    if(missing(R)) {
      stop("The constraint matrix R needs to be provided")
    } else if(!is.matrix(R)) {
      stop("The constraint matrix R has to be a matrix")
    } else if(nrow(as.matrix(R))!=p) {
      stop("The constraint matrix R has wrong dimension")
    } else if(ncol(as.matrix(R))>p) {
      stop("Why would you need a constraint matrix with more columns than rows? Please make sure your constraints make sense!")
    } else if(qr(as.matrix(R))$rank!=ncol(as.matrix(R))) {
      stop("The constraint matrix is not full column rank")
    }
  } else {
    if(missing(R)) {
      stop("a list of constraint matrices R needs to be provided")
    }
    if(!is.list(R) | length(R)!=M) {
      stop("The argument R should be a list of M constraint matrices R_{m} - one for each component model")
    }
    for(i1 in 1:M) {
      R0 = as.matrix(R[[i1]])
      q = ncol(R0)
      if(nrow(R0)!=p) {
        stop(paste("The constraint matrix R for regime", i1 ,"has wrong dimension"))
      } else if(q>p) {
        stop("Why would you need a constraint matrix with more columns than rows? Please make sure your constraints make sense!")
      } else if(qr(R0)$rank!=ncol(R0)) {
        stop(paste("The constraint matrix R for regime", i1 ,"is not full column rank"))
      }
    }
  }
}

#' @title Check p and M are correctly set
#'
#' @description \code{checkPM} checks that the arguments p and M are correctly set.
#' @inheritParams loglikelihood_int
#' @return Doesn't return anything, but throws and error if something is wrong.

checkPM <- function(p, M, GStMAR=FALSE) {
  if(GStMAR==TRUE) {
    if(length(M)!=2) {
      stop("For G-StMAR model the argument M should be a vector of length 2")
    }
    for(i1 in 1:2) {
      if(M[i1]<1 | M[i1]%%1!=0) {
        stop(paste0("Argument M[", i1, "] has to be positive integer"))
      }
    }
  } else if(M<1 | M%%1!=0) {
      stop("Argument M has to be positive integer")
  }
  if(p<1 | p%%1!=0) {
    stop("Argument p has to be positive integer")
  }
}

#' @title Calculate the number of parameters
#'
#' @description \code{nParams} calculates the number of parameters that should be in the parameter vector.
#' @inheritParams loglikelihood_int
#' @return returns the number of parameters.

nParams <- function(p, M, StMAR=FALSE, GStMAR=FALSE, restricted=FALSE, constraints=FALSE, R) {
  if(restricted==FALSE) {
    if(StMAR==TRUE) {
      if(constraints==FALSE) {
        d = M*(p+4)-1
      } else {
        d = 4*M-1+sum(vapply(1:M, function(i1) ncol(as.matrix(R[[i1]])), numeric(1)))
      }
    } else if(GStMAR==TRUE) {
      M1 = M[1]
      M2 = M[2]
      M = sum(M)
      if(constraints==FALSE) {
        d = M*(p+3)+M2-1
      } else {
        d = 3*M+M2-1+sum(vapply(1:M, function(i1) ncol(as.matrix(R[[i1]])), numeric(1)))
      }
    } else { # If GMAR==TRUE
      if(constraints==FALSE) {
        d = M*(p+3)-1
      } else {
        d = 3*M-1+sum(vapply(1:M, function(i1) ncol(as.matrix(R[[i1]])), numeric(1)))
      }
    }
  } else { # if restricted==TRUE
    if(StMAR==TRUE) {
      if(constraints==FALSE) {
        d = 4*M+p-1
      } else {
        d = 4*M+ncol(as.matrix(R))-1
      }
    } else if(GStMAR==TRUE) {
      M1 = M[1]
      M2 = M[2]
      M = sum(M)
      if(constraints==FALSE) {
        d = 3*M+M2+p-1
      } else {
        d = 3*M+M2+ncol(as.matrix(R))-1
      }
    } else { # if GMAR=TRUE
      if(constraints==FALSE) {
        d = 3*M+p-1
      } else {
        d = 3*M+ncol(as.matrix(R))-1
      }
    }
  }
  return(d)
}


#' @title Check that the logical arguments StMAR and GStMAR don't conflict
#'
#' @description \code{checkLogicals} checks that the logical arguments StMAR and GStMAR don't conflict.
#' @inheritParams loglikelihood_int
#' @return Doesn't return anything, but throws and error if something is wrong.
checkLogicals <- function(StMAR, GStMAR) {
  if(StMAR==TRUE & GStMAR==TRUE) {
    stop("Arguments StMAR and GStMAR are both set to be TRUE. You obviously can't consider both models at the same time!")
  }
}

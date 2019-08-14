#' @title Check the stationary and identification conditions of specified GMAR, StMAR or G-StMAR model.
#'
#' @description \code{isStationary_int} checks the stationary condition and \code{isIdentifiable} checks the identification conditions
#'  of the specified GMAR, StMAR or G-StMAR model.
#'
#' @inheritParams loglikelihood_int
#' @param params a real valued parameter vector specifying the model.
#'  \describe{
#'    \item{For \strong{non-restricted} models:}{
#'      \describe{
#'        \item{For \strong{GMAR} model:}{Size \eqn{(M(p+3)-1x1)} vector \strong{\eqn{\theta}}\eqn{=}(\strong{\eqn{\upsilon_{1}}},...,\strong{\eqn{\upsilon_{M}}},
#'          \eqn{\alpha_{1},...,\alpha_{M-1}}), where \strong{\eqn{\upsilon_{m}}}\eqn{=(\phi_{m,0},}\strong{\eqn{\phi_{m}}}\eqn{,
#'          \sigma_{m}^2)} and \strong{\eqn{\phi_{m}}}=\eqn{(\phi_{m,1},...,\phi_{m,p}), m=1,...,M}.}
#'        \item{For \strong{StMAR} model:}{Size \eqn{(M(p+4)-1x1)} vector (\strong{\eqn{\theta, \nu}})\eqn{=}(\strong{\eqn{\upsilon_{1}}},...,\strong{\eqn{\upsilon_{M}}},
#'          \eqn{\alpha_{1},...,\alpha_{M-1}, \nu_{1},...,\nu_{M}}).}
#'        \item{For \strong{G-StMAR} model:}{Size \eqn{(M(p+3)+M2-1x1)} vector (\strong{\eqn{\theta, \nu}})\eqn{=}(\strong{\eqn{\upsilon_{1}}},...,\strong{\eqn{\upsilon_{M}}},
#'          \eqn{\alpha_{1},...,\alpha_{M-1}, \nu_{M1+1},...,\nu_{M}}).}
#'      }
#'    }
#'    \item{For \strong{restricted} models:}{
#'      \describe{
#'        \item{For \strong{GMAR} model:}{Size \eqn{(3M+p-1x1)} vector \strong{\eqn{\theta}}\eqn{=(\phi_{1,0},...,\phi_{M,0},}\strong{\eqn{\phi}}\eqn{,
#'          \sigma_{1}^2,...,\sigma_{M}^2,\alpha_{1},...,\alpha_{M-1})}, where \strong{\eqn{\phi}}=\eqn{(\phi_{1},...,\phi_{M})}.}
#'        \item{For \strong{StMAR} model:}{Size \eqn{(4M+p-1x1)} vector (\strong{\eqn{\theta, \nu}})\eqn{=(\phi_{1,0},...,\phi_{M,0},}\strong{\eqn{\phi}}\eqn{,
#'          \sigma_{1}^2,...,\sigma_{M}^2,\alpha_{1},...,\alpha_{M-1}, \nu_{1},...,\nu_{M})}.}
#'        \item{For \strong{G-StMAR} model:}{Size \eqn{(3M+M2+p-1x1)} vector (\strong{\eqn{\theta, \nu}})\eqn{=(\phi_{1,0},...,\phi_{M,0},}\strong{\eqn{\phi}}\eqn{,
#'          \sigma_{1}^2,...,\sigma_{M}^2,\alpha_{1},...,\alpha_{M-1}, \nu_{M1+1},...,\nu_{M})}.}
#'      }
#'    }
#'  }
#'  Symbol \eqn{\phi} denotes an AR coefficient, \eqn{\sigma^2} a variance, \eqn{\alpha} a mixing weight and \eqn{\nu} a degrees of
#'  freedom parameter. In the \strong{G-StMAR} model the first \code{M1} components are \emph{GMAR-type} and the rest \code{M2} components
#'  are \emph{StMAR-type}.
#'  Note that in the case \strong{M=1} the parameter \eqn{\alpha} is dropped, and in the case of \strong{StMAR} or \strong{G-StMAR} model
#'  the degrees of freedom parameters \eqn{\nu_{m}} have to be larger than \eqn{2}.
#' @details These functions don't support models parametrized with general linear constraints.
#' @return Returns \code{TRUE} or \code{FALSE} accordingly.
#' @section Warning:
#'  These functions don't have any argument checks!
#' @references
#'  \itemize{
#'    \item Kalliovirta L., Meitz M. and Saikkonen P. 2015. Gaussian Mixture Autoregressive model for univariate time series.
#'            \emph{Journal of Time Series Analysis}, \strong{36}, 247-266.
#'    \item Meitz M., Preve D., Saikkonen P. 2018. A mixture autoregressive model based on Student's t-distribution.
#'            arXiv:1805.04010 \strong{[econ.EM]}.
#'    \item There are currently no published references for the G-StMAR model, but it's a straightforward generalization with
#'            theoretical properties similar to the GMAR and StMAR models.
#'  }

isStationary_int <- function(p, M, params, restricted=FALSE) {
  M <- sum(M)
  if(restricted == FALSE) {
    pars <- matrix(params[1:(M*(p + 2))], ncol=M)
    for(i1 in 1:M) {
      if(any(abs(polyroot(c(1, -pars[2:(p + 1), i1]))) <= 1 + 1e-8)) {
        return(FALSE)
      }
    }
  } else {
    absroots <- abs(polyroot(c(1, -params[(M + 1):(M + p)])))
    if(any(absroots <= 1 + 1e-8)) {
      return(FALSE)
    }
  }
  TRUE
}


#' @rdname isStationary_int
isIdentifiable <- function(p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL) {
  model <- match.arg(model)
  params <- removeAllConstraints(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints)
  M_orig <- M
  if(model == "G-StMAR") {
    M1 <- M[1]
    M2 <- M[2]
    M <- sum(M)
    if(M1 == 1 & M2 == 1) {
      return(TRUE)
    }
  }
  if(M == 1) {
    return(TRUE)
  }
  pars <- pick_pars(p=p, M=M_orig, params=params, model=model, restricted=FALSE, constraints=NULL)
  alphas <- pick_alphas(p=p, M=M_orig, params=params, model=model, restricted=FALSE, constraints=NULL)
  dfs <- pick_dfs(p=p, M=M_orig, params=params, model=model)

  if(model == "StMAR") {
    pars <- rbind(pars, dfs)
  } else if(model == "G-StMAR") {
    pars1 <- as.matrix(pars[,1:M1])
    pars2 <- as.matrix(pars[,(M1 + 1):M])
    pars2 <- rbind(pars2, dfs)
    alphas1 <- alphas[1:M1]
    alphas2 <- alphas[(M1 + 1):M]
  }

  alphas_unsorted <- function(alps) is.unsorted(rev(alps), strictly=TRUE)
  pars_dublicates <- function(prs0) anyDuplicated(t(prs0)) != 0

  if(model == "G-StMAR") { # Check GMAR parameters and StMAR-parameters separately
    if(M1 > 1) {
      if(alphas_unsorted(alphas1)) {
        return(FALSE)
      } else if(pars_dublicates(pars1)) {
        return(FALSE)
      }
    }
    if(M2 > 1) {
      if(alphas_unsorted(alphas2)) {
        return(FALSE)
      } else if(pars_dublicates(pars2)) {
        return(FALSE)
      }
    }
  } else { # if model != "G-StMAR
    if(alphas_unsorted(alphas)) {
      return(FALSE)
    } else if(pars_dublicates(pars)) {
      return(FALSE)
    }
  }
  TRUE
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
#' isStationary(1, 2, params12t, model="StMAR")
#'
#' # G-StMAR model
#' params12gs <- c(1, 0.1, 1, 2, 0.2, 2, 0.8, 20)
#' isStationary(1, c(1, 1), params12gs, model="G-StMAR")
#'
#' # Restricted GMAR model
#' params13r <- c(0.1, 0.2, 0.3, -0.99, 0.1, 0.2, 0.3, 0.5, 0.3)
#' isStationary(1, 3, params13r, restricted=TRUE)
#'
#' # Restricted StMAR model
#' params22tr <- c(-0.1, -0.2, 0.01, 0.99, 0.3, 0.4, 0.9, 3, 13)
#' isStationary(2, 2, params22tr, model="StMAR", restricted=TRUE)
#'
#' # Restricted G-StMAR model
#' params13gsr <- c(1, 2, 3, -0.99, 1, 2, 3, 0.5, 0.4, 20, 30)
#' isStationary(1, c(1, 2), params13gsr, model="G-StMAR", restricted=TRUE)
#'
#' # GMAR model as a mixture of AR(2) and AR(1) models
#' constraints <- list(diag(1, ncol=2, nrow=2), as.matrix(c(1, 0)))
#' params22c <- c(1.2, 0.8, 0.2, 0.3, 3.3, 0.7, 3, 0.8)
#' isStationary(2, 2, params22c, constraints=constraints)
#'
#' # Such StMAR(3,2) that the AR coefficients are restricted to be the
#' # same for both regimes and that the second AR coefficients are
#' # constrained to zero.
#' params32trc <- c(1, 2, 0.8, -0.3, 1, 2, 0.7, 11, 12)
#' isStationary(3, 2, params32trc, model="StMAR", restricted=TRUE,
#'              constraints=matrix(c(1, 0, 0, 0, 0, 1), ncol=2))
#' @export

isStationary <- function(p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL) {
  model <- match.arg(model)
  checkPM(p=p, M=M, model=model)
  check_params_length(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints)
  if(!is.null(constraints)) {
    checkConstraintMat(p=p, M=M, restricted=restricted, constraints=constraints)
    params <- reformConstrainedPars(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints)
  }
  isStationary_int(p=p, M=M, params=params, restricted=restricted)
}


#' @title Check the data is set correctly and correct if not
#'
#' @description \code{checkAndCorrectData} checks that the data is set correctly and corrects it if not.
#'  Throws an error if it can't convert the data to the correct form.
#'
#' @inheritParams loglikelihood_int
#' @return Returns a numeric column matrix containing the data.

checkAndCorrectData <- function(data, p) {
  if(anyNA(data)) {
    stop("The data contains NA values")
  } else if(length(data) < p+2) {
    stop("The data must contain at least p+2 values")
  }
  if(is.matrix(data)) {
    if(ncol(data) > 1) {
      stop("Only univariate time series are supported! For multivariate analysis, try the package 'gmvarkit'.")
    }
  }
  if(!is.ts(data)) {
    data <- as.ts(data)
  }
  data
}


#' @title Check the parameter vector is specified correctly
#'
#' @description \code{parameterChecks} checks dimension, restrictions and stationarity of the given parameters
#'   of GMAR, StMAR or G-StMAR model. Throws an error if any check fails. Does NOT consider the identifiability
#'   condition!
#'
#' @inheritParams loglikelihood_int
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
#' @return Throws an informative error if any check fails. Doesn't return anything.

parameterChecks <- function(p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL) {
  model <- match.arg(model)
  check_params_length(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints)
  params <- removeAllConstraints(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints)
  pars <- pick_pars(p=p, M=M, params=params, model=model, restricted=FALSE, constraints=NULL)
  alphas <- pick_alphas(p=p, M=M, params=params, model=model, restricted=FALSE, constraints=NULL)
  dfs <- pick_dfs(p=p, M=M, params=params, model=model)

  if(model == "StMAR" | model == "G-StMAR") {
    if(any(dfs <= 2)) {
      stop("The degrees of freedom parameters have to be larger than 2")
    } else if(any(dfs > 1e+5)) {
      stop("We have set an upper bound of 1e+5 for the degrees of freedom parameters
           in order to avoid numerical problems. This is not, however, restrictive
           since t distribution with dfs higher than this strongly resembles
           the Gaussian distribution by its shape.")
    }
  }

  if(sum(M) >= 2 & sum(alphas[-sum(M)]) >= 1) {
    stop("The mixing weights don't sum to one")
  } else if(any(alphas <= 0)) {
    stop("Mixing weight parameters have to be larger than zero")
  } else if(any(pars[p + 2,] <= 0)) {
    stop("Variance parameters have to be larger than zero")
  } else if(!isStationary_int(p=p, M=M, params=params, restricted=FALSE)) {
    stop("The model doesn't satisfy the stationary condition")
  }
}


#' @title Check the constraint matrices
#'
#' @description \code{checkConstraintMat} checks for some parts that the constraint matrices are correctly set.
#' @inheritParams loglikelihood_int
#' @return Doesn't return anything, but throws an informative error if finds out that something is wrong.

checkConstraintMat <- function(p, M, restricted=FALSE, constraints=NULL) {
  if(!is.null(constraints)) {
    M <- sum(M) # For G-StMAR
    if(restricted == TRUE) { # The constraints is a matrix
      if(!is.matrix(constraints)) {
        stop("The constraint matrix has to be a matrix (not in a list)")
      } else if(nrow(as.matrix(constraints)) != p) {
        stop("The constraint matrix has wrong dimension")
      } else if(ncol(as.matrix(constraints)) > p) {
        stop("The constraint matrix has more columns than rows?? Please make sure your constraints make sense!")
      } else if(qr(as.matrix(constraints))$rank != ncol(as.matrix(constraints))) {
        stop("The constraint matrix doesn't have full column rank")
      }
    } else { # The constraints is a list of matrices
      if(!is.list(constraints) | length(constraints) != M) {
        stop("The argument constraints should be a list of M constraint matrices - one for each mixture component")
      }
      for(i1 in 1:M) {
        C0 <- as.matrix(constraints[[i1]])
        q <- ncol(C0)
        if(nrow(C0) != p) {
          stop(paste("The constraint matrix for regime", i1 ,"has wrong dimension"))
        } else if(q > p) {
          stop(paste("The constraint matrix for regime", i1, "has more columns than rows?? Please make sure your constraints make sense!"))
        } else if(qr(C0)$rank != ncol(C0)) {
          stop(paste("The constraint matrix for regime", i1 ,"doesn't have full column rank"))
        }
      }
    }
  }
}


#' @title Check p and M are correctly set
#'
#' @description \code{checkPM} checks that the arguments p and M are correctly set.
#' @inheritParams loglikelihood_int
#' @return Doesn't return anything, but throws an informative error if something is wrong.

checkPM <- function(p, M, model=c("GMAR", "StMAR", "G-StMAR")) {
  model <- match.arg(model)
  if(model == "G-StMAR") {
    if(length(M) != 2 | !all_pos_ints(M)) {
      stop("For G-StMAR model the argument M should be length 2 a positive integer vector")
    }
  } else if(!all_pos_ints(M) | length(M) != 1) {
      stop("Argument M has to be positive integer")
  }
  if(!all_pos_ints(p)) {
    stop("Argument p has to be positive integer")
  }
}


#' @title Calculate the number of parameters
#'
#' @description \code{nParams} calculates the number of parameters that should be in the parameter vector.
#' @inheritParams loglikelihood_int
#' @return Returns the number of parameters.

nParams <- function(p, M, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL) {
  model <- match.arg(model)
  if(restricted == FALSE) {
    n_const_pars <- function(M, constraints) sum(vapply(1:M, function(i1) ncol(as.matrix(constraints[[i1]])), numeric(1)))
    if(model == "StMAR") {
      if(is.null(constraints)) {
        d <- M*(p + 4) - 1
      } else {
        d <- 4*M - 1 + n_const_pars(M=M, constraints=constraints)
      }
    } else if(model == "G-StMAR") {
      M1 <- M[1]
      M2 <- M[2]
      M <- sum(M)
      if(is.null(constraints)) {
        d <- M*(p + 3) + M2 - 1
      } else {
        d <- 3*M + M2 - 1 + n_const_pars(M=M, constraints=constraints)
      }
    } else { # If model == "GMAR"
      if(is.null(constraints)) {
        d <- M*(p + 3) - 1
      } else {
        d <- 3*M - 1 + n_const_pars(M=M, constraints=constraints)
      }
    }
  } else { # if restricted == TRUE
    if(model == "StMAR") {
      if(is.null(constraints)) {
        d <- 4*M + p - 1
      } else {
        d <- 4*M + ncol(as.matrix(constraints)) - 1
      }
    } else if(model == "G-StMAR") {
      M1 <- M[1]
      M2 <- M[2]
      M <- sum(M)
      if(is.null(constraints)) {
        d <- 3*M + M2 + p - 1
      } else {
        d <- 3*M + M2 + ncol(as.matrix(constraints)) - 1
      }
    } else { # if model == "GMAR"
      if(is.null(constraints)) {
        d <- 3*M + p - 1
      } else {
        d <- 3*M + ncol(as.matrix(constraints)) - 1
      }
    }
  }
  d
}



#' @title Check whether all arguments are positive scalar whole numbers
#'
#' @description \code{all_pos_ints} tells whether all the elements in a vector
#'   are strictly positive whole numbers.
#'
#' @param x a vector containing the elements to be tested.
#' @return Returns \code{TRUE} or \code{FALSE} accordingly.

all_pos_ints <- function(x) {
  all(vapply(x, function(x)  x %% 1 == 0 & length(x) == 1 & x >= 1, logical(1)))
}


#' @title Check that the argument model is correctly specified.
#'
#' @description \code{check_model} checks that the argument model is correctly specified.
#'
#' @inheritParams loglikelihood_int
#' @return Doesn't return anything, but throws and error if something is wrong.

check_model <- function(model) {
  if(!model %in% c("GMAR", "StMAR", "G-StMAR")) {
    stop("The argument 'model' has to be 'GMAR', 'StMAR' or 'G-StMAR'")
  }
}

#' @title Check that the parameter vector has the correct dimension
#'
#' @description \code{check_model} checks that the parameter vector has the correct dimension
#' @inheritParams loglikelihood_int
#' @inherit check_model return

check_params_length <- function(p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL) {
  if(length(params) != nParams(p=p, M=M, model=model, restricted=restricted, constraints=constraints)) {
    stop("The parameter vector has wrong dimension")
  }
}


#' @title Check that given object has class attribute 'gsmar'
#'
#' @description \code{check_gsmar} checks that that given object has class attribute 'gsmar'.
#'
#' @param object an object to be tested
#' @inherit check_model return

check_gsmar <- function(object) {
  if(!any(class(object) == "gsmar")) stop("The argument 'gsmar' has to be object of class 'gsmar' created with fitGSMAR or GSMAR.")
}


#' @title Check that given object contains data
#'
#' @description \code{check_data} checks that that given object contains data.
#'
#' @inheritParams check_gsmar
#' @inherit check_gsmar return

check_data <- function(object) {
  if(is.null(object$data)) stop("The model has to contain data! Data can be added with the function add_data.")
}



#' @title Warn about large degrees of freedom paramater values
#'
#' @description \code{warn_dfs} warns if the model contains large degrees of freedom paramater values
#'   possibly indicating unrealible numerical derivatives.
#'
#' @inheritParams check_gsmar
#' @inheritParams loglikelihood_int
#' @param warn_about warn about inaccurate derivatives or standard errors?
#' @details Either provide a class 'gsmar' object or specify the model by hand.
#' @return Doesn't return anything but throws a warning if any degrees of freedom parameters have value
#'   larger than 1000.

warn_dfs <- function(object, p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL, warn_about=c("derivs", "errors")) {

  if(!missing(object)) {
    p <- object$model$p
    M <- object$model$M
    params <- object$params
    model <- object$model$model
    restricted <- object$model$restricted
    constraints <- object$model$constraints
  }
  if(model %in% c("StMAR", "G-StMAR")) {
    pars <- removeAllConstraints(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints)
    dfs <- pick_dfs(p=p, M=M, params=pars, model=model)
    if(any(dfs > 100)) warning("The model contains overly large degrees of freedom parameter values. Consider switching to G-StMAR model by setting the corresponding regimes to GMAR type with the function 'stmar_to_gstmar'.")
  }
}

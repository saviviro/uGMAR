#' @title Check the stationarity and identification conditions of specified GMAR, StMAR, or G-StMAR model.
#'
#' @description \code{is_stationary_int} checks the stationarity condition and \code{is_identifiable} checks the identification condition
#'  of the specified GMAR, StMAR, or G-StMAR model.
#'
#' @inheritParams loglikelihood_int
#' @details \code{is_stationary_int} does not support models imposing linear constraints. In order to use it for a model imposing linear
#'  constraints, one needs to expand the constraints first to obtain a nonconstrained parameters vector.
#'
#'  Note that \code{is_stationary_int} returns \code{FALSE} for stationary parameter vectors if they are extremely close to the boundary
#'  of the stationarity region.
#'
#'  \code{is_identifiable} checks that the regimes are sorted according to the mixing weight parameters and that there are no dublicate
#'  regimes.
#' @return Returns \code{TRUE} or \code{FALSE} accordingly.
#' @section Warning:
#'  These functions don't have any argument checks!
#' @references
#'  \itemize{
#'    \item Kalliovirta L., Meitz M. and Saikkonen P. 2015. Gaussian Mixture Autoregressive model for univariate time series.
#'            \emph{Journal of Time Series Analysis}, \strong{36}, 247-266.
#'    \item Meitz M., Preve D., Saikkonen P. 2018. A mixture autoregressive model based on Student's t-distribution.
#'            arXiv:1805.04010 \strong{[econ.EM]}.
#'    \item Virolainen S. 2020. A mixture autoregressive model based on Gaussian and Student's t-distribution.	arXiv:2003.05221 [econ.EM].
#'  }

is_stationary_int <- function(p, M, params, restricted=FALSE) {
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


#' @rdname is_stationary_int
is_identifiable <- function(p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL) {
  model <- match.arg(model)
  params <- remove_all_constraints(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints)
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

  if(model == "G-StMAR") { # Check GMAR parameters and StMAR parameters separately
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


#' @title Check the stationary condition of specified GMAR, StMAR, or G-StMAR model.
#'
#' @description \code{is_stationary} checks the stationarity condition of the specified GMAR, StMAR, or G-StMAR model.
#'
#' @inheritParams loglikelihood
#' @details This function falsely returns \code{FALSE} for stationary models when the parameter is extremely close
#'  to the boundary of the stationarity region.
#' @inherit is_stationary_int return references
#' @examples
#' # GMAR model
#' params22 <- c(0.4, 0.39, 0.6, 0.3, 0.4, 0.1, 0.6, 0.3, 0.8)
#' is_stationary(p=2, M=2, params=params22)
#'
#' # StMAR model
#' params12t <- c(-0.3, 1, 0.9, 0.1, 0.8, 0.6, 0.7, 10, 12)
#' is_stationary(p=1, M=2, params=params12t, model="StMAR")
#'
#' # G-StMAR model
#' params12gs <- c(1, 0.1, 1, 2, 0.2, 2, 0.8, 20)
#' is_stationary(p=1, M=c(1, 1), params=params12gs, model="G-StMAR")
#'
#' # Restricted GMAR model
#' params13r <- c(0.1, 0.2, 0.3, -0.99, 0.1, 0.2, 0.3, 0.5, 0.3)
#' is_stationary(p=1, M=3, params=params13r, restricted=TRUE)
#'
#' # Such StMAR(3, 2) that the AR coefficients are restricted to be the
#' # same for both regimes and that the second AR coefficients are
#' # constrained to zero.
#' params32trc <- c(1, 2, 0.8, -0.3, 1, 2, 0.7, 11, 12)
#' is_stationary(p=3, M=2, params=params32trc, model="StMAR", restricted=TRUE,
#'               constraints=matrix(c(1, 0, 0, 0, 0, 1), ncol=2))
#' @export

is_stationary <- function(p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL) {
  model <- match.arg(model)
  check_pM(p=p, M=M, model=model)
  check_params_length(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints)
  if(!is.null(constraints)) {
    check_constraint_mat(p=p, M=M, restricted=restricted, constraints=constraints)
    params <- reform_constrained_pars(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints)
  }
  is_stationary_int(p=p, M=M, params=params, restricted=restricted)
}


#' @title Check that the data is set correctly and correct if not
#'
#' @description \code{check_and_correct_data} checks that the data is set correctly and
#'  throws an error if there is something wrong with the data.
#'
#' @inheritParams loglikelihood_int
#' @return Returns the data as a class 'ts' object.

check_and_correct_data <- function(data, p) {
  if(anyNA(data)) {
    stop("The data contains NA values")
  } else if(length(data) < p+2) {
    stop("The data must contain at least p+2 values")
  } else if(!is.numeric(data)) {
    stop("The data should be numeric")
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
#' @description \code{parameter_checks} checks dimension, restrictions, and stationarity of the given parameter
#'   of a GMAR, StMAR, or G-StMAR model. Throws an error if any check fails. Does NOT consider the identifiability
#'   condition!
#'
#' @inheritParams loglikelihood_int
#' @return Throws an informative error if any check fails. Does not return anything.

parameter_checks <- function(p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL) {
  model <- match.arg(model)
  check_params_length(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints)
  params <- remove_all_constraints(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints)
  pars <- pick_pars(p=p, M=M, params=params, model=model, restricted=FALSE, constraints=NULL)
  alphas <- pick_alphas(p=p, M=M, params=params, model=model, restricted=FALSE, constraints=NULL)
  dfs <- pick_dfs(p=p, M=M, params=params, model=model)

  if(model == "StMAR" | model == "G-StMAR") {
    if(any(dfs <= 2)) {
      stop("The degrees of freedom parameters have to be larger than 2")
    } else if(any(dfs > 1e+5)) {
      stop("We have set an upper bound of 1e+5 for the degrees of freedom parameters
           in order to avoid numerical problems. This is not, however, restrictive
           since t-distribution with such a large degrees of freedom parameter value
           already strongly resembles the Gaussian distribution by its shape. Moreover,
           the conditional variance of the corresponding regime is close to the constant
           parameter 'sigma^2' if the degrees of freedom parameter is very large.")
    }
  }

  if(sum(M) >= 2 & sum(alphas[-sum(M)]) >= 1) {
    stop("The mixing weight parameters don't sum to one")
  } else if(any(alphas <= 0)) {
    stop("The mixing weight parameters have to be strictly larger than zero")
  } else if(any(pars[p + 2,] <= 0)) {
    stop("The variance parameters have to be strictly larger than zero")
  } else if(!is_stationary_int(p=p, M=M, params=params, restricted=FALSE)) {
    stop("The model doesn't satisfy the stationarity condition")
  }
}


#' @title Check the constraint matrices
#'
#' @description \code{check_constraint_mat} checks for some parts that the constraint matrices are correctly set.
#' @inheritParams loglikelihood_int
#' @return Doesn't return anything but throws an informative error if finds out that something is wrong.

check_constraint_mat <- function(p, M, restricted=FALSE, constraints=NULL) {
  if(!is.null(constraints)) {
    M <- sum(M)
    if(restricted == TRUE) { # 'constraints' is a single matrix
      if(!is.matrix(constraints)) {
        stop("The constraint matrix has to be a matrix (not in a list)")
      } else if(nrow(as.matrix(constraints)) != p) {
        stop("The constraint matrix has wrong dimension")
      } else if(ncol(as.matrix(constraints)) > p) {
        stop("The constraint matrix has more columns than rows!? Please make sure your constraints make sense!")
      } else if(qr(as.matrix(constraints))$rank != ncol(as.matrix(constraints))) {
        stop("The constraint matrix does not have full column rank")
      }
    } else { # restricted == FALSE, so 'constraints' is a list of matrices
      if(!is.list(constraints) | length(constraints) != M) {
        stop("The argument constraints should be a list of M constraint matrices - one for each mixture component")
      }
      for(i1 in 1:M) {
        C0 <- as.matrix(constraints[[i1]])
        q <- ncol(C0)
        if(nrow(C0) != p) {
          stop(paste("The constraint matrix for regime", i1 ,"has wrong dimension"))
        } else if(q > p) {
          stop(paste("The constraint matrix for regime", i1, "has more columns than rows!? Please make sure your constraints make sense!"))
        } else if(qr(C0)$rank != ncol(C0)) {
          stop(paste("The constraint matrix for regime", i1 ,"does not have full column rank"))
        }
      }
    }
  }
}


#' @title Check that p and M are correctly set
#'
#' @description \code{check_pM} checks that the arguments p and M are correctly set.
#' @inheritParams loglikelihood_int
#' @return Doesn't return anything but throws an informative error if something is wrong.

check_pM <- function(p, M, model=c("GMAR", "StMAR", "G-StMAR")) {
  model <- match.arg(model)
  if(model == "G-StMAR") {
    if(length(M) != 2 | !all_pos_ints(M)) {
      stop("For a G-StMAR model the argument M should be a length 2 positive integer vector")
    }
  } else if(!all_pos_ints(M) | length(M) != 1) {
      stop("Argument M has to be positive integer")
  }
  if(!all_pos_ints(p) | length(p) != 1) {
    stop("Argument p has to be positive integer")
  }
}


#' @title Calculate the number of parameters
#'
#' @description \code{n_params} calculates the number of parameters that should be in the parameter vector.
#' @inheritParams loglikelihood_int
#' @return Returns the number of parameters.

n_params <- function(p, M, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL) {
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



#' @title Check whether all arguments are stricly positive natural numbers
#'
#' @description \code{all_pos_ints} tells whether all the elements in a vector
#'   are strictly positive natural numbers.
#'
#' @param x a vector containing the elements to be tested.
#' @return Returns \code{TRUE} or \code{FALSE} accordingly.

all_pos_ints <- function(x) {
  all(vapply(x, function(x)  x %% 1 == 0 & length(x) == 1 & x >= 1, logical(1)))
}


#' @title Check that the argument 'model' is correctly specified.
#'
#' @description \code{check_model} checks that the argument 'model' is correctly specified.
#'
#' @inheritParams loglikelihood_int
#' @return Doesn't return anything but throws and error if something is wrong.

check_model <- function(model) {
  if(!model %in% c("GMAR", "StMAR", "G-StMAR")) {
    stop("The argument 'model' has to be 'GMAR', 'StMAR', or 'G-StMAR'")
  }
}

#' @title Check that the parameter vector has the correct dimension
#'
#' @description \code{check_model} checks that the parameter vector has the correct dimension.
#' @inheritParams loglikelihood_int
#' @inherit check_model return

check_params_length <- function(p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL) {
  if(length(params) != n_params(p=p, M=M, model=model, restricted=restricted, constraints=constraints)) {
    stop("The parameter vector has wrong dimension")
  }
}


#' @title Check that given object has class attribute 'gsmar'
#'
#' @description \code{check_gsmar} checks that the given object has class attribute 'gsmar'.
#'
#' @param object an object to be tested
#' @param object_name the name of the tested object
#' @inherit check_model return

check_gsmar <- function(object, object_name="gsmar") {
  if(!any(class(object) == "gsmar")) stop(paste("The argument", object_name, "has to be an object of class 'gsmar', often created with function 'fitGSMAR' or 'GSMAR'."))
}


#' @title Check that given object contains data
#'
#' @description \code{check_data} checks that a given object contains data.
#'
#' @inheritParams check_gsmar
#' @inherit check_gsmar return

check_data <- function(object) {
  if(is.null(object$data)) stop("The model has to contain data! Data can be added with the function add_data.")
}



#' @title Warn about large degrees of freedom parameter values
#'
#' @description \code{warn_dfs} warns if the model contains large degrees of freedom parameter values.
#'
#' @inheritParams check_gsmar
#' @inheritParams loglikelihood_int
#' @param warn_about warn about inaccurate derivatives or standard errors?
#' @details Either provide a class 'gsmar' object or specify the model by hand.
#' @return Doesn't return anything but throws a warning if any degrees of freedom parameters have value
#'   larger than 100.

warn_dfs <- function(object, p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), warn_about=c("derivs", "errors")) {

  if(!missing(object)) {
    M <- object$model$M
    params <- object$params
    model <- object$model$model
  }
  if(model %in% c("StMAR", "G-StMAR")) {
    M2 <- ifelse(model == "G-StMAR", M[2], M)
    d <- length(params)
    dfs <- params[(d - M2 + 1):d]
    if(any(dfs > 100)) {
      warning("The model contains overly large degrees of freedom parameter values. Consider switching to a G-StMAR model by setting the corresponding regimes to be GMAR type with the function 'stmar_to_gstmar'.")
    }
  }
}


#' @title Warn about near-unit-roots in some regimes
#'
#' @description \code{warn_ar_roots} warns if the model contains near-unit-roots in some regimes
#'
#' @inheritParams simulateGSMAR
#' @param tol if some root is smaller that \code{1 + tol}, a warning is thrown
#' @details Warns if some moduli of the autoregressive polynomial's roots are close to one.
#' @return Doesn't return anything.

warn_ar_roots <- function(gsmar, tol=0.005) {
  ar_roots <- get_ar_roots(gsmar)
  near_nonstat <- vapply(1:sum(gsmar$model$M), function(i1) any(abs(ar_roots[[i1]]) < 1 + tol), logical(1))
  if(any(near_nonstat)) {
    my_string <- ifelse(sum(near_nonstat) == 1,
                        paste("Regime", which(near_nonstat),"has near-unit-roots!"),
                        paste("Regimes", paste(which(near_nonstat), collapse=" and ") ,"have near-unit-roots!"))
    warning(paste(my_string, "Consider building a model from the next-largest local maximum with the function 'alt_gsmar' by adjusting its argument 'which_largest'."))
  }
}

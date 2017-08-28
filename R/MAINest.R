#' @import stats
#'
#' @title Estimate Gaussian or Student's t Mixture Autoregressive model
#'
#' @description \code{fitGMAR} estimates GMAR or StMAR model in two phases. It uses genetic algorithm to find parameter values close to the maximum point
#'  of the log-likelihood function and then uses them as starting values for quasi-Newton method to find the maximum point.
#' @inheritParams GAfit
#' @param nCalls an (optional) positive integer specifying how many rounds of estimation should be done.
#'  The estimation results may vary from round to round because of multimodality of the log-likelihood function
#'  and randomness associated with the genetic algorithm. Default is \code{round(10 + 9*log(M))}.
#' @param multicore an (optional) logical argument defining whether parallel computing should be used in the estimation process.
#'  Highly recommended and default is \code{TRUE}.
#' @param ncores an (optional) positive integer specifying the number of cores to be used in the estimation process.
#'  Default is that the number of available cores is detected with \code{parallel::detectCores()} and all them are used. Ignored if \code{multicore==FALSE}.
#' @param printRes an (optional) logical argument defining whether results should be printed or not. Default is \code{TRUE}.
#' @param runTests an (optional) logical argument defining whether quantile residual tests for the estimated model should be performed or not. Default is \code{FALSE}.
#' @details The user should consider adjusting \code{ar0scale} and/or \code{sigmascale} accordingly to the best knowledge about the process.
#'
#'   Note that \code{fitGMAR} can't verify whether the found estimates denote the global or just a local maximum point.
#'   For more reliable results one should increase the number of estimation rounds (\code{nCalls}) to be performed.
#' @return Returns a list with...
#' \describe{
#'   \item{\code{$estimates}}{The estimated parameter vector...
#'  \describe{
#'    \item{For \strong{non-restricted} models:}{
#'      \describe{
#'        \item{For \strong{GMAR} model:}{Size \eqn{(M(p+3)-1x1)} vector \strong{\eqn{\theta}}\eqn{=}(\strong{\eqn{\upsilon_{1}}},...,\strong{\eqn{\upsilon_{M}}},
#'          \eqn{\alpha_{1},...,\alpha_{M-1}}), where \strong{\eqn{\upsilon_{m}}}\eqn{=(\phi_{m,0},}\strong{\eqn{\phi_{m}}}\eqn{,
#'          \sigma_{m}^2)} and \strong{\eqn{\phi_{m}}}=\eqn{(\phi_{m,1},...,\phi_{m,p}), m=1,...,M}.}
#'        \item{For \strong{StMAR} model:}{Size \eqn{(M(p+4)-1x1)} vector (\strong{\eqn{\theta, \nu}})\eqn{=}(\strong{\eqn{\upsilon_{1}}},...,\strong{\eqn{\upsilon_{M}}},
#'          \eqn{\alpha_{1},...,\alpha_{M-1}, \nu_{1},...,\nu_{M}}).}
#'        \item{With \strong{linear constraints}:}{Parameter vector as descripted above, but vectors \strong{\eqn{\phi_{m}}} replaced with vectors \strong{\eqn{\psi_{m}}}
#'          that satisfy \strong{\eqn{\phi_{m}}}\eqn{=}\strong{\eqn{R_{m}\psi_{m}}} for all \eqn{m=1,...,M}, where
#'          \strong{\eqn{\psi_{m}}}\eqn{=(\psi_{m,1},...,\psi_{m,q_{m}})}.}
#'      }
#'    }
#'    \item{For \strong{restricted} models:}{
#'      \describe{
#'        \item{For \strong{GMAR} model:}{Size \eqn{(3M+p-1x1)} vector \strong{\eqn{\theta}}\eqn{=(\phi_{1,0},...,\phi_{M,0},}\strong{\eqn{\phi}}\eqn{,
#'          \sigma_{1}^2,...,\sigma_{M}^2,\alpha_{1},...,\alpha_{M-1})}, where \strong{\eqn{\phi}}=\eqn{(\phi_{1},...,\phi_{M})}.}
#'        \item{For \strong{StMAR} model:}{Size \eqn{(4M+p-1x1)} vector (\strong{\eqn{\theta, \nu}})\eqn{=(\phi_{1,0},...,\phi_{M,0},}\strong{\eqn{\phi}}\eqn{,
#'          \sigma_{1}^2,...,\sigma_{M}^2,\alpha_{1},...,\alpha_{M-1}, \nu_{1},...,\nu_{M})}.}
#'        \item{With \strong{linear constraints}:}{Parameter vector as descripted above, but vector \strong{\eqn{\phi}} replaced with vector \strong{\eqn{\psi}}
#'          that satisfies \strong{\eqn{\phi}}\eqn{=}\strong{\eqn{R\psi}}, where
#'          \strong{\eqn{\psi}}\eqn{=(\psi_{1},...,\psi_{q})}.}
#'      }
#'    }
#'  }
#'  Symbol \eqn{\phi} denotes an AR coefficient, \eqn{\sigma^2} a variance, \eqn{\alpha} a mixing weight and \eqn{\nu} a degrees of
#'  freedom parameter.}
#'   \item{\code{$stdErrors}}{Approximate standard errors of the estimates. \code{NA} values may sometimes occur because the observed information matrix is numerically estimated.}
#'   \item{\code{$loglikelihood}}{Log-likelihood value of the estimated model.}
#'   \item{\code{$IC}}{A data frame containing information criteria scores of the estimated model: \code{$AIC}, \code{$BIC}, \code{$HQIC}.}
#'   \item{\code{$quantileResiduals}}{A numeric vector containing the quantile residuals of the estimated model.}
#'   \item{\code{$mixingWeights}}{A numeric matrix containing the mixing weights of the estimated model (i:th column for i:th regime).}
#'   \item{\code{$allEstimates}}{A list of estimated parameter vectors from all of the estimation rounds.}
#'   \item{\code{$allLoglikelihoods}}{A numeric vector containing the log-likelihood values from all of the estimation rounds. Corresponds to \code{$allEstimates}.}
#'   \item{\code{$converged}}{A logical vector containing information whether the quasi-Newton algorithm converged successfully or not. Corresponds to \code{$allEstimates}.}
#'   \item{\code{$normality}}{A data frame containing results from the normality test. Returned only if \code{runTests==TRUE}.}
#'   \item{\code{$autocorrelation}}{A data frame containing results from the autocorrelation tests. Returned only if \code{runTests==TRUE}.}
#'   \item{\code{$cond.heteroscedasticity}}{A data frame containing results from the conditional heteroscedasticity tests. Returned only if \code{runTests==TRUE}.}
#'   \item{\code{$unconstrainedEstimates}}{A numeric parameter vector denoting the estimates without any constraints (if given any). That is instead of
#'     vectors \strong{\eqn{\psi_{m}}} the estimates are parametrized with vectors \strong{\eqn{\phi_{m}}} calculated from
#'     \strong{\eqn{\phi_{m}}}\eqn{=}\strong{\eqn{R_{m}\psi_{m}}}, or in the case of restricted models
#'     \strong{\eqn{\phi}}\eqn{=}\strong{\eqn{R\psi}}. Returned only if \code{constraints==TRUE}.}
#' }
#' @section Printed results:
#'   The results printed out regarding the genetic algorithm and quasi-Newton estimations are the log-likelihood values
#'   the algorithms ended up with. The lowest value, mean value and largest value are printed to give perspective.
#'
#'   If quantile residual tests are run, the results from the tests are printed so that the letter "N" means normality test, "A" autocorrelation test
#'   and "H" conditional heteroscedasticity test. The numbers right next to "A" and "H" indicate the number of lags used
#'   in each test. The statistics following them are the corresponding test statistics and p-values.
#'   \code{NA} values mean that it was not (numerically) possible for the code to calculate all the necessary estimates for the tests.
#' @section Suggested packages:
#'   Install the suggested package "pbapply" if you wish to see a progress bar during parallel computing.
#'
#'   For faster evaluation of the quantile residuals of StMAR model install the suggested package "gsl".
#'   Note that for large StMAR models with large data the evaluations for the quantile residual tests may take
#'   significantly long time without the package "gsl".
#' @section The optimization algorithms:
#'   The genetic algorithm is mostly based on the description by \emph{Dorsey R. E. ja Mayer W. J. (1995)}. It uses individually
#'   adaptive crossover and mutation rates described by \emph{Patnaik L.M. and Srinivas M. (1994)}, with slight modifications.
#'
#'   The quasi-Newton method is implemented with function \code{optim} from the package \code{stats}.
#' @references
#'  \itemize{
#'    \item Kalliovirta L., Meitz M. and Saikkonen P. (2015) Gaussian Mixture Autoregressive model for univariate time series.
#'          \emph{Journal of Time Series Analysis}, \strong{36}, 247-266.
#'    \item Kalliovirta L. (2012) Misspecification tests based on quantile residuals.
#'          \emph{The Econometrics Journal}, \strong{15}, 358-393.
#'    \item Dorsey R. E. ja Mayer W. J. (1995) Genetic algorithms for estimation problems with multiple optima, nondifferentiability, and other irregular features.
#'          \emph{Journal of Business & Economic Statistics}, \strong{13}, 53-66.
#'    \item Patnaik L.M. and Srinivas M. (1994) Adaptive Probabilities of Crossover and Mutation in Genetic Algorithms.
#'          \emph{Transactions on Systems, Man and Cybernetics} \strong{24}, 656-667.
#'    \item Lutkepohl H. New Introduction to Multiple Time Series Analysis,
#'          \emph{Springer}, 2005.
#'    \item Galbraith, R., Galbraith, J., (1974). On the inverses of some patterned matrices arising
#'            in the theory of stationary time series. \emph{Journal of Applied Probability} \strong{11}, 63-71.
#'    \item References regarding the StMAR model and general linear constraints will be updated after they are published.
#'  }
#' @examples
#' \donttest{
#' # GMAR model
#' fit12 <- fitGMAR(VIX, 1, 2, ar0scale=c(3, 2), runTests=TRUE)
#'
#' # Restricted GMAR model
#' fit12r <- fitGMAR(VIX, 1, 2, restricted=TRUE, nCalls=10,
#'                   runTests=TRUE)
#'
#' # StMAR model
#' fit12t <- fitGMAR(VIX, 1, 2, StMAR=TRUE, ar0scale=c(3, 2))
#'
#' # Non-mixture version of StMAR model
#' fit11t <- fitGMAR(VIX, 1, 1, StMAR=TRUE)
#'
#' # Fit GMAR model that is a mixture of AR(1) and such AR(3) model that the
#' # second AR coeffiecient is constrained to zero.
#' R <- list(matrix(c(1, 0, 0, 0, 0, 1), ncol=2), as.matrix(c(1, 0, 0)))
#' fit32c <- fitGMAR(VIX, 3, 2, constraints=TRUE, R=R, ar0scale=c(3, 2))
#'
#' # Fit such constrained StMAR(3, 1) model that the second order AR coefficient
#' # is constrained to zero.
#' R0 <- matrix(c(1, 0, 0, 0, 0, 1), ncol=2)
#' fit31tc <- fitGMAR(VIX, 3, 1, StMAR=TRUE, constraints=TRUE, R=list(R0))
#'
#' # Fit such StMAR(3,2) that the AR coefficients are restricted to be
#' # the same for both regimes and that the second AR coefficients are
#' # constrained to zero.
#' fit32trc <- fitGMAR(VIX, 3, 2, StMAR=TRUE, restricted=TRUE, constraints=TRUE,
#'                     R=matrix(c(1, 0, 0, 0, 0, 1), ncol=2))
#' }
#' @export

fitGMAR <- function(data, p, M, StMAR=FALSE, restricted=FALSE, constraints=FALSE, R, conditional=TRUE, nCalls, multicore=TRUE, ncores, initpop=FALSE, ngen, popsize, smartMu, ar0scale, sigmascale, printRes=TRUE, runTests=FALSE) {

  checkPM(p, M)
  data = checkAndCorrectData(data, p)
  epsilon = round(log(.Machine$double.xmin)+10)
  if(constraints==TRUE) {
    checkConstraintMat(p, M, R, restricted=restricted)
  }
  if(missing(R)) {
    R = NULL
  }
  d = nParams(p=p, M=M, StMAR=StMAR, restricted=restricted, constraints=constraints, R=R)
  if(missing(popsize)) {
    popsize = 10*d
  }
  if(missing(ngen)) {
    ngen = max(round(0.1*length(data)), 200)
  }
  if(missing(smartMu)) {
    smartMu = min(100, round(0.5*ngen))
  }
  if(missing(nCalls)) {
    nCalls = round(10 + 9*log(M))
  } else if(nCalls<1 | nCalls%%1!=0) {
    stop("nCalls has to be positive integer")
  }
  if(missing(ar0scale)) {
    avg = mean(data); T1 = length(data)
    c0 = (t(data-avg)%*%(data-avg))/T1
    c1 = (t((data-avg)[1:(T1-1)])%*%(data-avg)[2:T1])/(T1-1)
    ar0scale = c(1.5*avg*(1-c1/c0), max(c0, 4))
  } else if(length(ar0scale)!=2) {
    stop("ar0scale is wrong dimension")
  } else if(ar0scale[2]<=0) {
    stop("the second element of ar0scale should be larger than zero")
  }
  if(missing(sigmascale)) {
    sigmascale = 1+sd(data)
  } else if(length(sigmascale)!=1) {
    stop("sigmascale is wrong dimension")
  } else if(sigmascale<=0) {
    stop("sigmascale should be larger than zero")
  }
  if(missing(ncores)) {
    if(multicore==TRUE) {
      ncores = parallel::detectCores()
    } else {
      ncores = 1
    }
  } else if(multicore==FALSE) {
    ncores = 1
  } else if(ncores<1 | ncores%%1!=0) {
    stop("ncores has to be positive integer")
  }
  if(nCalls<ncores) {
    ncores = nCalls
  }
  print(paste("Using", ncores, "cores for", nCalls, "estimation rounds..."))

  #############
  ## FITTING ##
  #############

  ### Genetic algorithm ###

  if(multicore==FALSE) {
    print("Optimizing with genetic algorithm...")
    GAfit_tmp <- function(round) {
      ret = GAfit(data, p, M, StMAR=StMAR, restricted=restricted, constraints=constraints, R=R, conditional=conditional,
                  ngen=ngen, popsize=popsize, smartMu=smartMu, ar0scale=ar0scale, sigmascale=sigmascale, epsilon=epsilon)
      print(paste0(round, "/", nCalls))
      return(ret)
    }
    GAresults = lapply(1:nCalls, function(round) GAfit_tmp(round))
  } else {
    cl = parallel::makeCluster(ncores)
    parallel::clusterExport(cl, c("data", "p", "M", "ngen", "popsize", "smartMu", "ar0scale", "sigmascale", "conditional",
                                  "restricted", "constraints", "R", "initpop", "GAfit", "loglikelihood_int", "isStationary_int",
                                  "isIdentifiable", "StMAR", "sortComponents", "smartIndividual_int", "randomIndividual_int",
                                  "checkAndCorrectData", "reformParameters", "parameterChecks", "reformConstrainedPars",
                                  "checkConstraintMat", "mixingWeights_int", "extractRegime", "changeRegime", "epsilon",
                                  "nParams"), envir = environment())
    parallel::clusterEvalQ(cl, c(library(Brobdingnag)))

    print("Optimizing with genetic algorithm...")
    if(requireNamespace("pbapply", quietly = TRUE)) {
      parallel::clusterEvalQ(cl, library(pbapply))
      GAresults = pbapply::pblapply(1:nCalls, function(x) GAfit(data=data, p=p, M=M, StMAR=StMAR, restricted=restricted,
                                                                constraints=constraints, R=R, conditional=conditional,
                                                                ngen=ngen, popsize=popsize, smartMu=smartMu, ar0scale=ar0scale,
                                                                sigmascale=sigmascale, initpop=initpop, epsilon=epsilon), cl=cl)
    } else {
      GAresults = parallel::parLapply(cl, 1:nCalls, function(x) GAfit(data=data, p=p, M=M, StMAR=StMAR, restricted=restricted,
                                                                      constraints=constraints, R=R, conditional=conditional,
                                                                      ngen=ngen, popsize=popsize, smartMu=smartMu, ar0scale=ar0scale,
                                                                      sigmascale=sigmascale, initpop=initpop, epsilon=epsilon))
    }
    parallel::stopCluster(cl=cl)
  }
  loks = sapply(1:nCalls, function(i1) loglikelihood_int(data=data, p=p, M=M, params=GAresults[[i1]], StMAR=StMAR, restricted=restricted,
                                                        constraints=constraints, R=R, conditional=conditional, boundaries=TRUE,
                                                        checks=FALSE, returnTerms=FALSE, epsilon=epsilon))
  if(printRes==TRUE) {
    print("Results from genetic algorithm:")
    print(paste("lowest value:", round(min(loks), 3)))
    print(paste("mean value:", round(mean(loks), 3)))
    print(paste("largest value:", round(max(loks), 3)))
  }

  ### Quasi-Newton ###

  # Logarithmize dfs to get them to the same range as other parameters
  manipulateDFS <- function(M, params, logDFS=TRUE) {
    dfs = params[(d-M+1):d]
    if(logDFS==TRUE) {
      params[(d-M+1):d] = log(dfs) # log dfs
    } else {
      params[(d-M+1):d] = exp(dfs) # exp dfs (from log to normal)
    }
    return(params)
  }
  if(StMAR==TRUE) {
    GAresults = lapply(1:nCalls, function(i1) manipulateDFS(M=M, params=GAresults[[i1]], logDFS=TRUE))
  }

  # Function to maximize loglikelihood
  f = function(params) {
    if(StMAR==TRUE) {
      params = manipulateDFS(M=M, params=params, logDFS=FALSE) # Unlogarithmize dfs
    }
    loglikelihood_int(data=data, p=p, M=M, params=params, StMAR=StMAR, restricted=restricted, constraints=constraints, R=R,
                  boundaries=TRUE, conditional=conditional, checks=FALSE, returnTerms=FALSE, epsilon=epsilon)
  }

  # Calculate gradient of the log-likelihood function using central finite difference
  h = 1e-8
  I = diag(rep(1, d))
  gr <- function(params) {
    grad = numeric(d)
    for(i1 in 1:d) {
      grad[i1] = (f(params+I[i1,]*h)-f(params-I[i1,]*h))/(2*h)
    }
    return(grad)
  }

  if(multicore==TRUE) {
    cl = parallel::makeCluster(ncores)
    parallel::clusterExport(cl, c("GAresults", "f", "manipulateDFS", "data", "p", "M", "conditional", "restricted", "constraints",
                                  "R", "loglikelihood_int", "isStationary_int", "isIdentifiable", "sortComponents",
                                  "checkAndCorrectData", "StMAR", "parameterChecks", "reformParameters",
                                  "reformConstrainedPars", "checkConstraintMat", "epsilon", "gr"), envir = environment())
    parallel::clusterEvalQ(cl, c(library(Brobdingnag)))

    print("Optimizing with quasi-Newton method...")
    if(requireNamespace("pbapply", quietly = TRUE)) {
      parallel::clusterEvalQ(cl, library(pbapply))
      NEWTONresults = pbapply::pblapply(1:nCalls, function(i1) optim(par=GAresults[[i1]], fn=f, gr=gr, method=c("BFGS"),
                                                                     control=list(fnscale=-1)), cl=cl)
    } else {
      NEWTONresults = parallel::parLapply(cl, 1:nCalls, function(i1) optim(par=GAresults[[i1]], fn=f, gr=gr, method=c("BFGS"),
                                                                           control=list(fnscale=-1)))
    }
    parallel::stopCluster(cl=cl)
  } else {
    print("optimizing with quasi-Newton method...")
    NEWTONresults = lapply(1:nCalls, function(i1) optim(par=GAresults[[i1]], fn=f, gr=gr, method=c("BFGS"), control=list(fnscale=-1)))
  }

  loks = vapply(1:nCalls, function(i1) NEWTONresults[[i1]]$value, numeric(1))
  converged = vapply(1:nCalls, function(i1) NEWTONresults[[i1]]$convergence==0, logical(1))
  if(constraints==FALSE) {
    newtonEstimates = lapply(1:nCalls, function(i1) sortComponents(p, M, NEWTONresults[[i1]]$par, StMAR=StMAR, restricted=restricted))
  } else { # To maintain the order of matrices R, models with constraints won't be sorted
    newtonEstimates = lapply(1:nCalls, function(i1) NEWTONresults[[i1]]$par)
  }

  # Unlogarithmize dfs
  if(StMAR==TRUE) {
    newtonEstimates = lapply(1:nCalls, function(i1) manipulateDFS(M=M, params=newtonEstimates[[i1]], logDFS=FALSE))
  }

  if(printRes==TRUE) {
    print("Results from quasi-Newton:")
    print(paste("lowest value:", round(min(loks), 3)))
    print(paste("mean value:", round(mean(loks), 3)))
    print(paste("largest value:", round(max(loks), 3)))
  }

  # Obtain the estimates
  bestind = which(loks==max(loks))[1]
  bestfit = NEWTONresults[[bestind]]
  params = newtonEstimates[[bestind]]
  if(constraints==FALSE) {
    params = sortComponents(p, M, params, StMAR=StMAR, restricted=restricted) # Sort the parameter vector by alphas
  }
  mw = mixingWeights_int(data, p, M, params, StMAR=StMAR, restricted=restricted, constraints=constraints, R=R, epsilon=epsilon)

  # Warnings and notifications
  if(any(vapply(1:M, function(i1) sum(mw[,i1]>0.05)<0.01*length(data), logical(1)))) {
    print("NOTE: At least one the regimes in the estimated model seems to be wasted! Consider re-estimating the model for steadier results or re-specify the model.")
  }
  if(bestfit$convergence==1) {
    print("NOTE: Iteration limit was reached when estimating the best fitting individual!")
  }

  #######################
  ## TESTS AND WRAP-UP ##
  #######################

  # Quantile residual tests
  if(runTests==TRUE) {
    print("Performing quantile residual tests...")
    qrTests = quantileResidualTests(data, p, M, params, StMAR=StMAR, restricted=restricted, constraints=constraints, R=R,
                                    lagsAC=c(1, 2, 4, 6, 8, 10), lagsCH=c(1, 2, 4, 6, 8, 10), nsimu=2000, printRes=printRes)
  }

  # Information criteria
  loglik = bestfit$value
  T0 = length(data) - p
  AIC = -2*loglik + 2*d
  BIC = -2*loglik + d*log(T0)
  HQIC = -2*loglik + 2*d*log(log(T0))
  IC = data.frame(AIC, BIC, HQIC)
  if(printRes==TRUE) {
    print(round(IC))
  }

  # Standard errors of the estimates
  stdErrors = standardErrors(data=data, p=p, M=M, params=params, StMAR=StMAR, restricted=restricted,
                             constraints=constraints, R=R, conditional=conditional, epsilon=epsilon)
  estimates = data.frame(params, stdErrors)
  colnames(estimates) = c("estimate", "stdError")
  if(printRes==TRUE) {
    print(round(estimates, 3))
  }

  # Collect the stuff to be returned
  results = list(params, stdErrors, loglik, IC,
                 quantileResiduals_int(data, p, M, params, StMAR=StMAR, restricted=restricted, constraints=constraints, R=R, epsilon=epsilon),
                 mw, newtonEstimates, loks, converged)
  names(results) = c("estimates", "stdErrors", "loglikelihood", "IC", "quantileResiduals", "mixingWeights",
                     "allEstimates", "allLoglikelihoods", "converged")
  if(runTests==TRUE) {
    results[[length(results)+1]] = qrTests$normality
    names(results)[length(results)] = "normality"
    results[[length(results)+1]] = qrTests$autocorrelation
    names(results)[length(results)] = "autocorrelation"
    results[[length(results)+1]] = qrTests$cond.heteroscedasticity
    names(results)[length(results)] = "cond.heteroscedasticity"
  }
  if(constraints==TRUE) {
    results[[length(results)+1]] = reformConstrainedPars(p=p, M=M, params=params, StMAR=StMAR, restricted=restricted, R=R)
    names(results)[length(results)] = "unconstrainedEstimates"
  }
  return(results)
}


#' @title Calculate standard errors for estimates of GMAR or StMAR model
#'
#' @description \code{standardErrors} numerically approximates standard errors for the given estimates of GMAR or StMAR model
#'
#' @inheritParams loglikelihood_int

standardErrors <- function(data, p, M, params, StMAR=FALSE, restricted=FALSE, constraints=FALSE, R, conditional=TRUE, epsilon) {

  # Function to derivate
  fn <- function(params) {
    loglikelihood_int(data=data, p=p, M=M, params=params, StMAR=StMAR, restricted=restricted, constraints=constraints, R=R,
                  boundaries=TRUE, conditional=conditional, checks=FALSE, returnTerms=FALSE, epsilon=epsilon)
  }
  h0 = c(6e-06, 1e-06, 0.001) # Difference
  d = length(params)
  I = diag(1, ncol=d, nrow=d) # Indicates which parameter is derivated

  for(j1 in 1:length(h0)) {
    h = h0[j1]
    Hess = matrix(ncol=d, nrow=d)

    # Calculate the second derivatives
    for(i1 in 1:d) {
      for(i2 in i1:d) {
        dr1 = (fn(params + h*I[i1,] + h*I[i2,]) - fn(params - h*I[i1,] + h*I[i2,]))/(2*h)
        dr2 = (fn(params + h*I[i1,] - h*I[i2,]) - fn(params - h*I[i1,] - h*I[i2,]))/(2*h)
        Hess[i1, i2] = (dr1 - dr2)/(2*h)
        Hess[i2, i1] = Hess[i1, i2] # Take use of symmetry
      }
    }

    # Inverse of the observed information matrix
    invObsInf <- tryCatch(solve(-Hess), error=function(cond) return(matrix(0, ncol=d, nrow=d)))

    # Calculate the standard errors if possible: break loop if all calculated and change the difference if not
    if(all(diag(invObsInf)>0) | j1==length(h0)) {
      stdErrors = sqrt(diag(invObsInf))
      if(all(stdErrors==0)) {
        stdErrors = rep(NA, d)
      }
      break;
    }
  }
  return(stdErrors)
}


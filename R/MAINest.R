#' @import stats
#'
#' @title Estimate Gaussian or Student's t Mixture Autoregressive model
#'
#' @description \code{fitGSMAR} estimates GMAR, StMAR or G-StMAR model in two phases: in the first phase it uses genetic algorithm
#'   to find starting values for gradient based variable metric algorithm (also known as quasi-Newton method), which it then uses
#'   to finalize the estimation in the second phase. Parallel computing is used to perform multiple rounds of estimations in parallel.
#' @inheritParams GAfit
#' @param nCalls a positive integer specifying how many rounds of estimation should be done.
#'  The estimation results may vary from round to round because of multimodality of the log-likelihood function
#'  and randomness associated with the genetic algorithm.
#' @param nCores a positive integer specifying the number of cores to be used in the estimation process.
#'  Default is that the number of available cores is detected with \code{parallel::detectCores()} and all of them are used.
#' @param maxit maximum number of iterations in the variable metric algorithm.
#' @param printRes should the estimation results be printed?
#' @param runTests should quantile residuals tests be performed after the estimation?
#' @param ... additional settings passed to the function \code{GAfit()} employing the genetic algorithm.
#' @details
#'  Because of complexity and multimodality of the log-likelihood function, it's \strong{not guaranteed} that the estimation
#'  algorithms will end up in the global maximum point. It's expected that most of the estimation rounds will end up in some local maximum
#'  point instead, and therefore a number of estimation rounds is required for reliable results. Because of the nature of the models,
#'  the estimation may fail particularly in the cases where the number of mixture components is chosen too large.
#'
#'  If the iteration limit in the variable metric algorithm (\code{maxit}) is reached, one can continue estimation by iterating
#'  more with the function \code{iterate_more()}.
#'
#'  The genetic algorithm is mostly based on the description by \emph{Dorsey and Mayer (1995)}.
#'  It uses (slightly modified) individually adaptive crossover and mutation rates described by \emph{Patnaik and Srinivas (1994)}
#'  and employs (50\%) fitness inheritance discussed by \emph{Smith, Dike and Stegmann (1995)}.
#'
#'  The variable metric algorithm (or quasi-Newton method) used in the second phase is implemented with function
#'  \code{optim} from the package \code{stats}.
#'
#'  Some mixture components of StMAR model may sometimes get very large degrees of freedom parameter estimates. Such estimates may, for
#'  example, cause computing the quantile residual tests to fail. However, such mixture components are very much similar to the components
#'  of GMAR model. It's hence advisable to further estimate a G-StMAR model by allowing the mixture components with large degrees of freedom
#'  estimates to be GMAR type.
#' @return Returns an object of class \code{'gsmar'} defining the estimated GMAR, StMAR or G-StMAR model. The returned object contains
#'   empirical mixing weights, quantile residuals and quantile residual test results if the tests were performed. In addition, the returned
#'   object contains the estimates and log-likelihood values from all the estimation rounds. The estimated parameter vector can be obtained
#'   at \code{gsmar$params} (and the corresponding approximate standard errors at \code{gsmar$std_errors}) and it is...
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
#'        \item{With \strong{linear constraints}:}{Replace the vectors \strong{\eqn{\phi_{m}}} with vectors \strong{\eqn{\psi_{m}}} and provide a  list of constraint
#'          matrices \strong{C} that satisfy \strong{\eqn{\phi_{m}}}\eqn{=}\strong{\eqn{R_{m}\psi_{m}}} for all \eqn{m=1,...,M}, where
#'          \strong{\eqn{\psi_{m}}}\eqn{=(\psi_{m,1},...,\psi_{m,q_{m}})}.}
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
#'        \item{With \strong{linear constraints}:}{Replace the vector \strong{\eqn{\phi}} with vector \strong{\eqn{\psi}} and provide a constraint matrix
#'          \strong{\eqn{C}} that satisfies \strong{\eqn{\phi}}\eqn{=}\strong{\eqn{R\psi}}, where
#'          \strong{\eqn{\psi}}\eqn{=(\psi_{1},...,\psi_{q})}.}
#'      }
#'    }
#'  }
#'  Symbol \eqn{\phi} denotes an AR coefficient, \eqn{\sigma^2} a variance, \eqn{\alpha} a mixing weight and \eqn{\nu} a degrees of
#'  freedom parameter. If \code{parametrization=="mean"} just replace each intercept term \eqn{\phi_{m,0}} with regimewise mean
#'  \eqn{\mu_m = \phi_{m,0}/(1-\sum\phi_{i,m})}. In the \strong{G-StMAR} model the first \code{M1} components are \emph{GMAR-type}
#'  and the rest \code{M2} components are \emph{StMAR-type}.
#'  Note that in the case \strong{M=1} the parameter \eqn{\alpha} is dropped, and in the case of \strong{StMAR} or \strong{G-StMAR} model
#'  the degrees of freedom parameters \eqn{\nu_{m}} have to be larger than \eqn{2}.
#' @section S3 methods:
#'  The following S3 methods are supported for class \code{'gsmar'} objects: \code{print}, \code{summary}, \code{plot},
#'  \code{logLik}, \code{residuals}.
#' @section Suggested packages:
#'  For faster evaluation of the quantile residuals of StMAR and G-StMAR models install the suggested package "gsl".
#'  Note that for large StMAR and G-StMAR models with large data the evaluations for the quantile residual tests may take
#'  significantly long time without the package "gsl".
#' @seealso \code{\link{GSMAR}}, \code{\link{iterate_more}}, \code{\link{add_data}}, \code{\link{swap_parametrization}},
#'  \code{\link{get_gradient}}, \code{\link{simulateGSMAR}}, \code{\link{predict.gsmar}}, \code{\link{diagnosticPlot}},
#'  , \code{\link{quantileResidualTests}}
#' @references
#'  \itemize{
#'    \item Dorsey R. E. and Mayer W. J. 1995. Genetic algorithms for estimation problems with multiple optima,
#'          nondifferentiability, and other irregular features. \emph{Journal of Business & Economic Statistics},
#'          \strong{13}, 53-66.
#'    \item Kalliovirta L., Meitz M. and Saikkonen P. 2015. Gaussian Mixture Autoregressive model for univariate time series.
#'          \emph{Journal of Time Series Analysis}, \strong{36}, 247-266.
#'    \item Meitz M., Preve D., Saikkonen P. 2018. A mixture autoregressive model based on Student's t-distribution.
#'          arXiv:1805.04010 \strong{[econ.EM]}.
#'    \item Patnaik L.M. and Srinivas M. 1994. Adaptive Probabilities of Crossover and Mutation in Genetic Algorithms.
#'          \emph{Transactions on Systems, Man and Cybernetics} \strong{24}, 656-667.
#'    \item Smith R.E., Dike B.A., Stegmann S.A. 1995. Fitness inheritance in genetic algorithms.
#'          \emph{Proceedings of the 1995 ACM Symposium on Applied Computing}, 345-350.
#'    \item There are currently no published references for G-StMAR model, but it's a straightforward generalization with
#'          theoretical properties similar to GMAR and StMAR models.
#'  }
#' @examples
#' \donttest{
#' # These are long running examples and use parallel computing
#'
#' # GMAR model
#' fit12 <- fitGSMAR(VIX, 1, 2, runTests=TRUE)
#' fit12
#' summary(fit12)
#' plot(fit12)
#'
#' # Restricted GMAR model
#' fit12r <- fitGSMAR(VIX, 1, 2, restricted=TRUE,
#'  parametrization="mean", nCalls=10)
#' fit12r
#' summary(fit12r)
#'
#' # Non-mixture version of StMAR model
#' fit11t <- fitGSMAR(VIX, 1, 1, model="StMAR", nCores=1, nCalls=1)
#' fit11t
#'
#' # StMAR model, 100 estimations rounds
#' fit12t <- fitGSMAR(VIX, 1, 2, model="StMAR", nCalls=100)
#' fit12t
#'
#' # Restricted StMAR model (implied by the StMAR(1,2) model)
#' fit12tr <- fitGSMAR(VIX, 1, 2, model="StMAR", restricted=TRUE)
#' fit12tr
#'
#' # G-StMAR model (implied by the StMAR(1,2) models), 100 estimation rounds
#' fit12gs <- fitGSMAR(VIX, 1, c(1, 1), model="G-StMAR", nCalls=100)
#' fit12gs
#'
#' # Restricted G-StMAR model (implied by the previous StMAR and G-StMAR models)
#' fit12gsr <- fitGSMAR(VIX, 1, c(1, 1), model="G-StMAR", restricted=TRUE)
#' fit12gsr
#'
#' # Fit GMAR model that is a mixture of AR(1) and such AR(3) model that the
#' # second AR coeffiecient is constrained to zero.
#' constraints <- list(matrix(c(1, 0, 0, 0, 0, 1), ncol=2), as.matrix(c(1, 0, 0)))
#' fit32c <- fitGSMAR(VIX, 3, 2, constraints=constraints)
#' fit32c
#'
#' # Fit such constrained StMAR(3, 1) model that the second order AR coefficient
#' # is constrained to zero.
#' constraints <- list(matrix(c(1, 0, 0, 0, 0, 1), ncol=2))
#' fit31tc <- fitGSMAR(VIX, 3, 1, model="StMAR", constraints=constraints)
#' fit31tc
#'
#' # Fit such StMAR(3,2) model that the AR coefficients are restricted to be
#' # the same for both regimes and that the second AR coefficients are
#' # constrained to zero.
#' fit32trc <- fitGSMAR(VIX, 3, 2, model="StMAR", restricted=TRUE,
#'                      constraints=matrix(c(1, 0, 0, 0, 0, 1), ncol=2))
#' fit32trc
#' }
#' @export

fitGSMAR <- function(data, p, M, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL, conditional=TRUE,
                     parametrization=c("intercept", "mean"), nCalls=round(10 + 9*log(sum(M))), nCores=min(nCalls, parallel::detectCores()),
                     maxit=100, printRes=TRUE, runTests=FALSE, ...) {
  on.exit(closeAllConnections())
  model <- match.arg(model)
  check_model(model)
  parametrization <- match.arg(parametrization)
  stopifnot(parametrization %in% c("intercept", "mean"))
  if(!all_pos_ints(c(nCalls, nCores, maxit))) stop("Arguments nCalls, nCores and maxit have to be positive integers")
  checkPM(p, M, model=model)
  data <- checkAndCorrectData(data, p)
  if(!is.null(constraints)) checkConstraintMat(p, M, restricted=restricted, constraints=constraints)
  d <- nParams(p=p, M=M, model=model, restricted=restricted, constraints=constraints)
  dot_params <- list(...)
  minval <- ifelse(is.null(dot_params$minval), -(10^(ceiling(log10(length(data))) + 1) - 1), dot_params$minval)
  red_criteria <- ifelse(rep(is.null(dot_params$red_criteria), 2), c(0.05, 0.01), dot_params$red_criteria)

  if(nCalls < nCores) {
    nCores <- nCalls
    message("nCores was set to be larger than the number of estimation rounds: using nCores = nCalls")
  }
  cat(paste("Using", nCores, "cores for", nCalls, "estimation rounds..."), "\n")


  ### Genetic algorithm optimization ###

  cl <- parallel::makeCluster(nCores)
  parallel::clusterExport(cl, ls(environment(fitGSMAR)), envir = environment(fitGSMAR)) # assign all variables from package:uGMAR
  parallel::clusterEvalQ(cl, c(library(Brobdingnag), library(pbapply)))

  cat("Optimizing with genetic algorithm...\n")
  GAresults <- pbapply::pblapply(1:nCalls, function(x) GAfit(data=data, p=p, M=M, model=model, restricted=restricted,
                                                             constraints=constraints, conditional=conditional,
                                                             parametrization=parametrization, ...), cl=cl)
  parallel::stopCluster(cl=cl)


  loks <- vapply(1:nCalls, function(i1) loglikelihood_int(data=data, p=p, M=M, params=GAresults[[i1]], model=model,
                                                          restricted=restricted, constraints=constraints, conditional=conditional,
                                                          parametrization=parametrization, boundaries=TRUE, checks=FALSE,
                                                          minval=minval), numeric(1))

  if(printRes == TRUE) {
    cat("Results from genetic algorithm:\n")
    cat(paste("lowest value: ", round(min(loks), 3)), "\n")
    cat(paste("mean value:   ", round(mean(loks), 3)), "\n")
    cat(paste("largest value:", round(max(loks), 3)), "\n")
  }

  ### Variable metric algorithm optimization ###

  # Logarithmize dfs to get them to the same range as other parameters:
  # this allows the df-parameters to "explode" more sensitively to implicate GMAR-type components.
  manipulateDFS <- function(M, params, model, FUN) {
    FUN <- match.fun(FUN)
    M2 <- ifelse(model == "StMAR", M, M[2])
    params[(d - M2 + 1):d] <- FUN(params[(d - M2 + 1):d])
    params
  }
  if(model == "StMAR" | model == "G-StMAR") {
    GAresults <- lapply(1:nCalls, function(i1) manipulateDFS(M=M, params=GAresults[[i1]], model=model, FUN=log))
  }

  # Function to maximize loglikelihood
  f <- function(params) {
     if(model == "StMAR" | model == "G-StMAR") {
       params <- manipulateDFS(M=M, params=params, model=model, FUN=exp) # Unlogarithmize dfs
     }
     tryCatch(loglikelihood_int(data=data, p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints,
                                conditional=conditional, parametrization=parametrization, boundaries=TRUE, checks=FALSE, minval=minval),
              error=function(e) minval)
  }

  # Calculate gradient of the log-likelihood function using central finite difference
  h <- 6e-6
  I <- diag(rep(1, d))
  gr <- function(params) {
    vapply(1:d, function(i1) (f(params + I[i1,]*h) - f(params - I[i1,]*h))/(2*h), numeric(1))
  }

  cl <- parallel::makeCluster(nCores)
  parallel::clusterExport(cl, ls(environment(fitGSMAR)), envir = environment(fitGSMAR)) # assign all variables from package:uGMAR
  parallel::clusterEvalQ(cl, c(library(Brobdingnag), library(pbapply)))
  cat("Optimizing with variable metric algorithm...\n")
  NEWTONresults <- pbapply::pblapply(1:nCalls, function(i1) optim(par=GAresults[[i1]], fn=f, gr=gr, method=c("BFGS"),
                                                                  control=list(fnscale=-1, maxit=maxit)), cl=cl)
  parallel::stopCluster(cl=cl)

  loks <- vapply(1:nCalls, function(i1) NEWTONresults[[i1]]$value, numeric(1))
  converged <- vapply(1:nCalls, function(i1) NEWTONresults[[i1]]$convergence == 0, logical(1))

  if(is.null(constraints)) {
    newtonEstimates <- lapply(1:nCalls, function(i1) sortComponents(p=p, M=M, params=NEWTONresults[[i1]]$par, model=model, restricted=restricted))
  } else {
    newtonEstimates <- lapply(1:nCalls, function(i1) NEWTONresults[[i1]]$par)
  }

  # Unlogarithmize dfs
  if(model == "StMAR" | model == "G-StMAR") {
    newtonEstimates <- lapply(1:nCalls, function(i1) manipulateDFS(M=M, params=newtonEstimates[[i1]], model=model, FUN=exp))
  }

  if(printRes == TRUE) {
    cat("Results from variable metric algorithm:\n")
    cat(paste("lowest value: ", round(min(loks), 3)), "\n")
    cat(paste("mean value:   ", round(mean(loks), 3)), "\n")
    cat(paste("largest value:", round(max(loks), 3)), "\n")
  }

  # Obtain the estimates
  bestind <- which(loks==max(loks))[1]
  bestfit <- NEWTONresults[[bestind]]
  params <- newtonEstimates[[bestind]]
  mw <- mixingWeights_int(data, p, M, params, model=model, restricted=restricted, constraints=constraints,
                          parametrization=parametrization, to_return="mw")

  # Warnings and notifications
  if(any(vapply(1:sum(M), function(i1) sum(mw[,i1] > red_criteria[1]) < red_criteria[2]*length(data), logical(1)))) {
    message("At least one of the mixture components in the estimated model seems to be wasted!")
  }
  if(bestfit$convergence == 1) message("Iteration limit was reached when estimating the best fitting individual!")


  ### Tests, estimates, standard errors, IC ###
  tmp_gsmar <- GSMAR(data, p, M, params=params, model=model, restricted=restricted, constraints=constraints,
                     conditional=conditional, parametrization=parametrization, calc_std_errors=FALSE)

  # Quantile residual tests
  if(runTests == TRUE) {
    cat("Performing quantile residual tests...\n")
    qr_tests <- quantileResidualTests(tmp_gsmar, lagsAC=c(1, 2, 5, 10), lagsCH=c(1, 2, 5, 10), nsimu=2000, printRes=printRes)
    if(printRes == TRUE) cat("\n")
  } else {
    qr_tests <- NULL
  }

  # Information criteria
  obs <- ifelse(conditional, length(data) - p, length(data))
  loglik <- bestfit$value
  IC <- get_IC(loglik=loglik, npars=d, obs=obs)

  # Standard errors of the estimates
  std_errors <- standardErrors(data=data, p=p, M=M, params=params, model=model, restricted=restricted,
                               constraints=constraints, conditional=conditional, parametrization=parametrization,
                               minval=minval)
  if(all(is.na(std_errors))) message("Unable to calculate approximate standard errors")


  ### Wrap up ###
  qresiduals <- quantileResiduals_int(data, p, M, params, model=model, restricted=restricted, constraints=constraints,
                                      parametrization=parametrization)

  ret <- structure(list(data=data,
                        model=list(p=p,
                                   M=M,
                                   model=model,
                                   restricted=restricted,
                                   constraints=constraints,
                                   conditional=conditional,
                                   parametrization=parametrization),
                        params=params,
                        std_errors=std_errors,
                        mixing_weights=mw,
                        quantile_residuals=qresiduals,
                        loglik=structure(loglik,
                                         class="logLik",
                                         df=d),
                        IC=IC,
                        all_estimates=newtonEstimates,
                        all_logliks=loks,
                        which_converged=converged,
                        qr_tests=qr_tests),
                   class="gsmar")
  cat("Finished!\n")
  ret
}


#' @title Maximum likelihood estimation of GMAR, StMAR or G-StMAR model with preliminary estimates
#'
#' @description \code{iterate_more} uses variable metric algorithm to finalize maximum likelihood
#'  estimation of GMAR, StMAR or G-StMAR model (object of class \code{'gsmarar'}) which already has preliminary estimates.
#'
#' @inheritParams simulateGSMAR
#' @inheritParams fitGSMAR
#' @details The main purpose of \code{iterate_more()} is to provide a simple and convenient tool to finalize
#'   the estimation when the maximum number of iterations is reached when estimating a model with the
#'   main estimation function \code{fitGSMAR()}. It's just a simple wrapper around function \code{optim()}
#'   from the package \code{stats} and \code{GSMAR()} from the package \code{uGMAR}.
#' @return Returns an object of class \code{'gsmar'} defining the estimated model. Can be used
#'   to work with other functions provided in \code{uGMAR}.
#' @seealso \code{fitGSMAR()}, \code{GSMAR()}, \code{optim()}
#' @inherit GSMAR references
#' @examples
#' \donttest{
#' # Estimate GMAR model with only 50 generations of genetic algorithm and
#' # only 1 iteration in variable metric algorithm
#' fit12 <- fitGSMAR(VIX, 1, 2, maxit=1, ngen=50)
#' fit12
#'
#' # Iterate more since iteration limit was reached
#' fit12 <- iterate_more(fit12)
#' fit12
#' }
#' @export

iterate_more <- function(gsmar, maxit=100) {
  check_gsmar(gsmar)
  stopifnot(maxit %% 1 == 0 & maxit >= 1)
  minval <- -(10^(ceiling(log10(length(gsmar$data))) + 1) - 1)

  fn <- function(params) {
    tryCatch(loglikelihood_int(data=gsmar$data, p=gsmar$model$p, M=gsmar$model$M, params=params,
                               model=gsmar$model$model, restricted=gsmar$model$restricted, constraints=gsmar$model$constraints,
                               conditional=gsmar$model$conditional, parametrization=gsmar$model$parametrization,
                               boundaries=TRUE, checks=FALSE, to_return="loglik", minval=minval), error=function(e) minval)
  }
  gr <- function(params) {
    calc_gradient(x=params, fn=fn)
  }

  res <- optim(par=gsmar$params, fn=fn, gr=gr, method=c("BFGS"), control=list(fnscale=-1, maxit=maxit))
  if(res$convergence == 1) message("The maximum number of iterations was reached! Consired iterating more.")

  GSMAR(data=gsmar$data, p=gsmar$model$p, M=gsmar$model$M, params=res$par,
        restricted=gsmar$model$restricted, constraints=gsmar$model$constraints,
        conditional=gsmar$model$conditional, parametrization=gsmar$model$parametrization,
        calc_std_errors=TRUE)
}

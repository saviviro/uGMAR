#' @import stats
#'
#' @title Estimate Gaussian or Student's t Mixture Autoregressive model
#'
#' @description \code{fitGSMAR} estimates GMAR, StMAR, or G-StMAR model in two phases. In the first phase, a genetic algorithm is employed
#'   to find starting values for a gradient based method. In the second phase, the gradient based variable metric algorithm is utilized to
#'   accurately converge to a local maximum or a saddle point near each starting value. Parallel computing is used to conduct multiple
#'   rounds of estimations in parallel.
#' @inheritParams GAfit
#' @param ncalls a positive integer specifying how many rounds of estimation should be conducted.
#'  The estimation results may vary from round to round because of multimodality of the log-likelihood function
#'  and the randomness associated with the genetic algorithm.
#' @param ncores the number of CPU cores to be used in the estimation process.
#' @param maxit the maximum number of iterations for the variable metric algorithm.
#' @param seeds a length \code{ncalls} vector containing the random number generator seed for each call to the genetic algorithm,
#'   or \code{NULL} for not initializing the seed. Exists for the purpose of creating reproducible results.
#' @param printRes should the estimation results be printed?
#' @param runTests should quantile residuals tests be performed after the estimation?
#' @param ... additional settings passed to the function \code{GAfit} employing the genetic algorithm.
#' @details
#'  Because of complexity and multimodality of the log-likelihood function, it's \strong{not guaranteed} that the estimation
#'  algorithm will end up in the global maximum point. It's often expected that most of the estimation rounds will end up in
#'  some local maximum point instead, and therefore a number of estimation rounds is required for reliable results. Because
#'  of the nature of the models, the estimation may fail particularly in the cases where the number of mixture components is
#'  chosen too large. Note that the genetic algorithm is designed to avoid solutions with mixing weights of some regimes
#'  too close to zero at almost all times ('redudant regimes') but the settings can, however, be adjusted (see ?GAfit).
#'
#'  If the iteration limit for the variable metric algorithm (\code{maxit}) is reached, one can continue the estimation by
#'  iterating more with the function \code{iterate_more}.
#'
#'  The core of the genetic algorithm is mostly based on the description by \emph{Dorsey and Mayer (1995)}. It utilizes
#'  a slightly modified version the individually adaptive crossover and mutation rates described by \emph{Patnaik and Srinivas (1994)}
#'  and employs (50\%) fitness inheritance discussed by \emph{Smith, Dike and Stegmann (1995)}. Large (in absolute value) but stationary
#'  AR parameter values are generated with the algorithm proposed by Monahan (1984).
#'
#'  The variable metric algorithm (or quasi-Newton method, Nash (1990, algorithm 21)) used in the second phase is implemented
#'  with function the \code{optim} from the package \code{stats}.
#'
#'  Some mixture components of the StMAR model may sometimes get very large estimates for the degrees of freedom parameters. Such parameters
#'  are weakly identified and induce various numerical problems. However, mixture components with large degree of freedom parameters are
#'  similar to the mixture components of the GMAR model. It's hence advisable to further estimate a G-StMAR model by allowing the mixture
#'  components with large degrees of freedom parameter estimates to be GMAR type.
#' @return Returns an object of class \code{'gsmar'} defining the estimated GMAR, StMAR or G-StMAR model. The returned object contains
#'   estimated mixing weights, some conditional and unconditional moments, quantile residuals, and quantile residual test results
#'   if the tests were performed. Note that the first p observations are taken as the initial values so the mixing weights, conditional
#'   moments, and quantile residuals start from the p+1:th observation (interpreted as t=1).In addition, the returned object contains
#'   the estimates and log-likelihood values from all of the estimation rounds.
#'   The estimated parameter vector can be obtained as \code{gsmar$params} (and the corresponding approximate standard errors as
#'   \code{gsmar$std_errors}) and it's...
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
#'          matrices \strong{C} that satisfy \strong{\eqn{\phi_{m}}}\eqn{=}\strong{\eqn{C_{m}\psi_{m}}} for all \eqn{m=1,...,M}, where
#'          \strong{\eqn{\psi_{m}}}\eqn{=(\psi_{m,1},...,\psi_{m,q_{m}})}.}
#'      }
#'    }
#'    \item{For \strong{restricted} models:}{
#'      \describe{
#'        \item{For \strong{GMAR} model:}{Size \eqn{(3M+p-1x1)} vector \strong{\eqn{\theta}}\eqn{=(\phi_{1,0},...,\phi_{M,0},}\strong{\eqn{\phi}}\eqn{,
#'          \sigma_{1}^2,...,\sigma_{M}^2,\alpha_{1},...,\alpha_{M-1})}, where \strong{\eqn{\phi}}=\eqn{(\phi_{1},...,\phi_{p})}.}
#'        \item{For \strong{StMAR} model:}{Size \eqn{(4M+p-1x1)} vector (\strong{\eqn{\theta, \nu}})\eqn{=(\phi_{1,0},...,\phi_{M,0},}\strong{\eqn{\phi}}\eqn{,
#'          \sigma_{1}^2,...,\sigma_{M}^2,\alpha_{1},...,\alpha_{M-1}, \nu_{1},...,\nu_{M})}.}
#'        \item{For \strong{G-StMAR} model:}{Size \eqn{(3M+M2+p-1x1)} vector (\strong{\eqn{\theta, \nu}})\eqn{=(\phi_{1,0},...,\phi_{M,0},}\strong{\eqn{\phi}}\eqn{,
#'          \sigma_{1}^2,...,\sigma_{M}^2,\alpha_{1},...,\alpha_{M-1}, \nu_{M1+1},...,\nu_{M})}.}
#'        \item{With \strong{linear constraints}:}{Replace the vector \strong{\eqn{\phi}} with vector \strong{\eqn{\psi}} and provide a constraint matrix
#'          \strong{\eqn{C}} that satisfies \strong{\eqn{\phi}}\eqn{=}\strong{\eqn{C\psi}}, where
#'          \strong{\eqn{\psi}}\eqn{=(\psi_{1},...,\psi_{q})}.}
#'      }
#'    }
#'  }
#'  Symbol \eqn{\phi} denotes an AR coefficient, \eqn{\sigma^2} a variance, \eqn{\alpha} a mixing weight, and \eqn{\nu} a degrees of
#'  freedom parameter. If \code{parametrization=="mean"} just replace each intercept term \eqn{\phi_{m,0}} with regimewise mean
#'  \eqn{\mu_m = \phi_{m,0}/(1-\sum\phi_{i,m})}. In the \strong{G-StMAR} model, the first \code{M1} components are \emph{GMAR type}
#'  and the rest \code{M2} components are \emph{StMAR type}.
#'  Note that in the case \strong{M=1} the parameter \eqn{\alpha} is dropped, and in the case of \strong{StMAR} or \strong{G-StMAR} model
#'  the degrees of freedom parameters \eqn{\nu_{m}} have to be larger than \eqn{2}.
#' @section S3 methods:
#'  The following S3 methods are supported for class \code{'gsmar'} objects: \code{print}, \code{summary}, \code{plot},
#'  \code{logLik}, \code{residuals}.
#' @section Suggested packages:
#'  For faster evaluation of the quantile residuals for StMAR and G-StMAR models, install the suggested package "gsl".
#'  Note that for large StMAR and G-StMAR models with large data, performing the quantile residual tests may take
#'  significantly long time without the package "gsl".
#' @seealso \code{\link{GSMAR}}, \code{\link{iterate_more}}, , \code{\link{stmar_to_gstmar}}, \code{\link{add_data}},
#'  \code{\link{profile_logliks}}, \code{\link{swap_parametrization}}, \code{\link{get_gradient}}, \code{\link{simulateGSMAR}}, \code{\link{predict.gsmar}},
#'   \code{\link{diagnostic_plot}}, \code{\link{quantile_residual_tests}}, \code{\link{cond_moments}}, \code{\link{uncond_moments}}, \code{\link{LR_test}}, \code{\link{Wald_test}}
#' @references
#'  \itemize{
#'    \item Dorsey R. E. and Mayer W. J. 1995. Genetic algorithms for estimation problems with multiple optima,
#'          nondifferentiability, and other irregular features. \emph{Journal of Business & Economic Statistics},
#'          \strong{13}, 53-66.
#'    \item Kalliovirta L., Meitz M. and Saikkonen P. 2015. Gaussian Mixture Autoregressive model for univariate time series.
#'          \emph{Journal of Time Series Analysis}, \strong{36}, 247-266.
#'    \item Meitz M., Preve D., Saikkonen P. 2018. A mixture autoregressive model based on Student's t-distribution.
#'          arXiv:1805.04010 \strong{[econ.EM]}.
#'    \item Monahan J.F. 1984. A Note on Enforcing Stationarity in Autoregressive-Moving Average Models.
#'          \emph{Biometrica} \strong{71}, 403-404.
#'    \item Nash J. 1990. Compact Numerical Methods for Computers. Linear algebra and Function Minimization.
#'          \emph{Adam Hilger}.
#'    \item Patnaik L.M. and Srinivas M. 1994. Adaptive Probabilities of Crossover and Mutation in Genetic Algorithms.
#'          \emph{Transactions on Systems, Man and Cybernetics} \strong{24}, 656-667.
#'    \item Smith R.E., Dike B.A., Stegmann S.A. 1995. Fitness inheritance in genetic algorithms.
#'          \emph{Proceedings of the 1995 ACM Symposium on Applied Computing}, 345-350.
#'    \item Virolainen S. 2020. A mixture autoregressive model based on Gaussian and Student's t-distribution.	arXiv:2003.05221 [econ.EM].
#'  }
#' @examples
#' \donttest{
#' # These are long running examples that use parallel computing
#'
#' # GMAR model
#' fit12 <- fitGSMAR(simudata, p=1, M=2, model="GMAR")
#' summary(fit12)
#' plot(fit12)
#' profile_logliks(fit12)
#'
#' # StMAR model
#' fit42 <- fitGSMAR(data=T10Y1Y, p=4, M=2, model="StMAR",
#'                   ncalls=16, ncores=4, seeds=1:16)
#' fit42
#' summary(fit42)
#' plot(fit42)
#'
#' # Restricted StMAR model: plot also the individual statistics with
#' # their approximate critical bounds using the given data
#' fit42r <- fitGSMAR(T10Y1Y, p=4, M=2, model="StMAR", restricted=TRUE,
#'                    ncores=4)
#' fit42r
#' plot(fit42)
#'
#' # Non-mixture version of StMAR model
#' fit101t <- fitGSMAR(T10Y1Y, p=10, M=1, model="StMAR", ncores=1, ncalls=1)
#' diagnostic_plot(fit101t)
#'
#' # G-StMAR model with one GMAR type and one StMAR type regime
#' fit42g <- fitGSMAR(T10Y1Y, p=4, M=c(1, 1), model="G-StMAR")
#' diagnostic_plot(fit42g)
#'
#' # GMAR model; seeds for reproducibility
#' fit43gm <- fitGSMAR(T10Y1Y, p=4, M=3, model="GMAR", ncalls=16,
#'   seeds=1:16)
#' fit43gm
#'
#' # Restricted GMAR model
#' fit43gmr <- fitGSMAR(T10Y1Y, p=4, M=3, model="GMAR", ncalls=12,
#'   restricted=TRUE, seeds=1:12)
#' fit43gmr
#'
#'
#' # The following three examples demonstrate how to apply linear constraints
#' # to the AR parameters.
#'
#' # Two-regime GMAR p=2 model with the second AR coeffiecient of
#' # of the second regime contrained to zero.
#' constraints <- list(diag(1, ncol=2, nrow=2), as.matrix(c(1, 0)))
#' fit22c <- fitGSMAR(T10Y1Y, p=2, M=2, constraints=constraints)
#' fit22c
#'
#' # Such constrained StMAR(3, 1) model that the second order AR coefficient
#' # is constrained to zero.
#' constraints <- list(matrix(c(1, 0, 0, 0, 0, 1), ncol=2))
#' fit31tc <- fitGSMAR(T10Y1Y, p=3, M=1, model="StMAR", constraints=constraints)
#' fit31tc
#'
#' # Such StMAR(3,2) that the AR coefficients are restricted to be
#' # the same for both regimes and that the second AR coefficients are
#' # constrained to zero.
#' fit32rc <- fitGSMAR(T10Y1Y, p=3, M=2, model="StMAR", restricted=TRUE,
#'  constraints=matrix(c(1, 0, 0, 0, 0, 1), ncol=2))
#' fit32rc
#' }
#' @export

fitGSMAR <- function(data, p, M, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL, conditional=TRUE,
                     parametrization=c("intercept", "mean"), ncalls=round(10 + 9*log(sum(M))), ncores=min(2, ncalls, parallel::detectCores()),
                     maxit=300, seeds=NULL, printRes=TRUE, runTests=FALSE, ...) {
  on.exit(closeAllConnections())
  if(!all_pos_ints(c(ncalls, ncores, maxit))) stop("Arguments ncalls, ncores and maxit have to be positive integers")
  if(!is.null(seeds) && length(seeds) != ncalls) stop("The argument 'seeds' needs be NULL or a vector of length 'ncalls'")
  model <- match.arg(model)
  check_model(model)
  parametrization <- match.arg(parametrization)
  check_pM(p, M, model=model)
  data <- check_and_correct_data(data, p)
  check_constraint_mat(p, M, restricted=restricted, constraints=constraints)
  d <- n_params(p=p, M=M, model=model, restricted=restricted, constraints=constraints)
  dot_params <- list(...)

  minval <- ifelse(is.null(dot_params$minval), get_minval(data), dot_params$minval)
  red_criteria <- ifelse(rep(is.null(dot_params$red_criteria), 2), c(0.05, 0.01), dot_params$red_criteria)

  if(ncores > parallel::detectCores()) {
    ncores <- parallel::detectCores()
    message("ncores was set to be larger than the number of detected: using ncores = parallel::detectCores()")
  }
  if(ncalls < ncores) {
    ncores <- ncalls
    message("ncores was set to be larger than the number of estimation rounds: using ncores = ncalls")
  }
  cat(paste("Using", ncores, "cores for", ncalls, "estimation rounds..."), "\n")


  ### Optimization with the genetic algorithm ###

  cl <- parallel::makeCluster(ncores)
  parallel::clusterExport(cl, ls(environment(fitGSMAR)), envir=environment(fitGSMAR)) # assign all variables from package:uGMAR
  parallel::clusterEvalQ(cl, c(library(Brobdingnag), library(pbapply)))

  cat("Optimizing with a genetic algorithm...\n")
  GAresults <- pbapply::pblapply(1:ncalls, function(i1) GAfit(data=data, p=p, M=M, model=model, restricted=restricted,
                                                             constraints=constraints, conditional=conditional,
                                                             parametrization=parametrization, seed=seeds[i1], ...), cl=cl)

  loks <- vapply(1:ncalls, function(i1) loglikelihood_int(data=data, p=p, M=M, params=GAresults[[i1]], model=model,
                                                          restricted=restricted, constraints=constraints, conditional=conditional,
                                                          parametrization=parametrization, boundaries=TRUE, checks=FALSE,
                                                          minval=minval), numeric(1))

  if(printRes) {
    printloks <- function() {
        printfun <- function(txt, FUN) cat(paste(txt, round(FUN(loks), 3)), "\n")
        printfun("The lowest loglik: ", min)
        printfun("The mean loglik:   ", mean)
        printfun("The largest loglik:", max)
    }
    cat("Results from the genetic algorithm:\n")
    printloks()
  }

  ### Optimization with the variable metric algorithm ###

  # Logarithmize dfs to get overly large dfs values to the same range as other parameters.
  # This adjusts the difference 'h' larger for larger dfs parameters in non-log scale to
  # avoid numerical problems associated with overly large degrees of freedom values.
  manipulateDFS <- function(M, params, model, FUN) {
    FUN <- match.fun(FUN)
    M2 <- ifelse(model == "StMAR", M, M[2])
    params[(d - M2 + 1):d] <- FUN(params[(d - M2 + 1):d])
    params
  }
  if(model == "StMAR" | model == "G-StMAR") {
    GAresults <- lapply(1:ncalls, function(i1) manipulateDFS(M=M, params=GAresults[[i1]], model=model, FUN=log))
  }

  # Function to maximize loglikelihood
  f <- function(params) {
     if(model == "StMAR" | model == "G-StMAR") {
       params <- manipulateDFS(M=M, params=params, model=model, FUN=exp) # Unlogarithmize dfs for calculating log-likelihood
     }
     tryCatch(loglikelihood_int(data=data, p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints,
                                conditional=conditional, parametrization=parametrization, boundaries=TRUE, checks=FALSE, minval=minval),
              error=function(e) minval)
  }

  # Calculate gradient of the log-likelihood function using central finite difference approximation
  I <- diag(rep(1, d))
  h <- 6e-6
  gr <- function(params) {
    vapply(1:d, function(i1) (f(params + I[i1,]*h) - f(params - I[i1,]*h))/(2*h), numeric(1))
  }

  cat("Optimizing with a variable metric algorithm...\n")
  NEWTONresults <- pbapply::pblapply(1:ncalls, function(i1) optim(par=GAresults[[i1]], fn=f, gr=gr, method="BFGS",
                                                                  control=list(fnscale=-1, maxit=maxit)), cl=cl)
  parallel::stopCluster(cl=cl)

  loks <- vapply(1:ncalls, function(i1) NEWTONresults[[i1]]$value, numeric(1))
  converged <- vapply(1:ncalls, function(i1) NEWTONresults[[i1]]$convergence == 0, logical(1))

  if(is.null(constraints)) {
    newtonEstimates <- lapply(1:ncalls, function(i1) sort_components(p=p, M=M, params=NEWTONresults[[i1]]$par, model=model, restricted=restricted))
  } else {
    newtonEstimates <- lapply(1:ncalls, function(i1) NEWTONresults[[i1]]$par)
  }

  # Unlogarithmize dfs
  if(model == "StMAR" | model == "G-StMAR") {
    newtonEstimates <- lapply(1:ncalls, function(i1) manipulateDFS(M=M, params=newtonEstimates[[i1]], model=model, FUN=exp))
  }

  if(printRes) {
    cat("Results from the variable metric algorithm:\n")
    printloks()
  }

  # Obtain the estimates
  bestind <- which(loks == max(loks))[1]
  bestfit <- NEWTONresults[[bestind]]
  params <- newtonEstimates[[bestind]]
  mw <- mixing_weights_int(data, p, M, params, model=model, restricted=restricted, constraints=constraints,
                          parametrization=parametrization, to_return="mw")

  # Warnings and notifications
  if(any(vapply(1:sum(M), function(i1) sum(mw[,i1] > red_criteria[1]) < red_criteria[2]*length(data), logical(1)))) {
    message("At least one of the mixture components in the estimated model seems to be wasted!")
  }
  if(bestfit$convergence == 1) message("Iteration limit was reached when estimating the best fitting individual! Iterate more with the function 'iterate_more'.")

  # Quantile residual tests
  if(runTests) {
    cat("Performing quantile residual tests...\n")
    tmp_gsmar <- GSMAR(data, p, M, params=params, model=model, restricted=restricted, constraints=constraints,
                       conditional=conditional, parametrization=parametrization, calc_std_errors=FALSE)
    qr_tests <- quantile_residual_tests(tmp_gsmar, lagsAC=c(1, 2, 5, 10), lagsCH=c(1, 2, 5, 10), nsimu=2000, printRes=printRes)
    if(printRes) cat("\n")
  } else {
    qr_tests <- NULL
  }

  ### Wrap up ###
  ret <- GSMAR(data=data, p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints,
               conditional=conditional, parametrization=parametrization, calc_qresiduals=TRUE,
               calc_cond_moments=TRUE, calc_std_errors=TRUE)
  ret$all_estimates <- newtonEstimates
  ret$all_logliks <- loks
  ret$which_converged <- converged
  ret$qrtests <- qr_tests

  cat("Finished!\n")
  ret
}


#' @title Maximum likelihood estimation of GMAR, StMAR, or G-StMAR model with preliminary estimates
#'
#' @description \code{iterate_more} uses a variable metric algorithm to finalize maximum likelihood
#'  estimation of a GMAR, StMAR or G-StMAR model (object of class \code{'gsmar'}) which already has
#'  preliminary estimates.
#'
#' @inheritParams simulateGSMAR
#' @inheritParams fitGSMAR
#' @inheritParams GSMAR
#' @details The main purpose of \code{iterate_more} is to provide a simple and convenient tool to finalize
#'   the estimation when the maximum number of iterations is reached when estimating a model with the
#'   main estimation function \code{fitGSMAR}. \code{iterate_more} is essentially a wrapper for the functions
#'   \code{optim} from the package \code{stats} and \code{GSMAR} from the package \code{uGMAR}.
#' @return Returns an object of class \code{'gsmar'} defining the estimated model.
#' @seealso \code{\link{fitGSMAR}}, \code{\link{GSMAR}}, \code{\link{stmar_to_gstmar}},
#'   \code{\link{profile_logliks}}, \code{\link{optim}}
#' @inherit GSMAR references
#' @examples
#' \donttest{
#' # Estimate GMAR model with only 50 generations of genetic algorithm and
#' # only 1 iteration in variable metric algorithm
#' fit13 <- fitGSMAR(T10Y1Y, 1, 3, maxit=1, ngen=50, ncalls=1, seeds=1)
#' fit13
#'
#' # Iterate more since iteration limit was reached
#' fit13 <- iterate_more(fit13)
#' fit13
#' }
#' @export

iterate_more <- function(gsmar, maxit=100, custom_h=NULL, calc_std_errors=TRUE) {
  check_gsmar(gsmar)
  stopifnot(maxit %% 1 == 0 & maxit >= 1)
  if(!is.null(custom_h)) stopifnot(length(custom_h) == length(gsmar$params))
  minval <- get_minval(gsmar$data)

  fn <- function(params) {
    tryCatch(loglikelihood_int(data=gsmar$data, p=gsmar$model$p, M=gsmar$model$M, params=params,
                               model=gsmar$model$model, restricted=gsmar$model$restricted, constraints=gsmar$model$constraints,
                               conditional=gsmar$model$conditional, parametrization=gsmar$model$parametrization,
                               boundaries=TRUE, checks=FALSE, to_return="loglik", minval=minval), error=function(e) minval)
  }
  gr <- function(params) {
    if(is.null(custom_h)) {
      varying_h <- get_varying_h(p=gsmar$model$p, M=gsmar$model$M, params=params, model=gsmar$model$model)
    } else {
      varying_h <- custom_h
    }
    calc_gradient(x=params, fn=fn, varying_h=varying_h)
  }

  res <- optim(par=gsmar$params, fn=fn, gr=gr, method=c("BFGS"), control=list(fnscale=-1, maxit=maxit))
  if(res$convergence == 1) message("The maximum number of iterations was reached! Consider iterating more.")

  GSMAR(data=gsmar$data, p=gsmar$model$p, M=gsmar$model$M, params=res$par, model=gsmar$model$model,
        restricted=gsmar$model$restricted, constraints=gsmar$model$constraints,
        conditional=gsmar$model$conditional, parametrization=gsmar$model$parametrization,
        calc_qresiduals=TRUE, calc_cond_moments=TRUE, calc_std_errors=calc_std_errors, custom_h=custom_h)
}



#' @title Returns the default smallest allowed log-likelihood for given data.
#'
#' @description \code{get_minval} returns the default smallest allowed log-likelihood for given data.
#'
#' @inheritParams GAfit
#' @details This function exists simply to avoid dublication inside the package.
#' @return Returns \code{-(10^(ceiling(log10(length(data))) + 1) - 1)}
#' @seealso \code{\link{fitGSMAR}}, \code{\link{GAfit}}

get_minval <- function(data) {
  -(10^(ceiling(log10(length(data))) + 1) - 1)
}

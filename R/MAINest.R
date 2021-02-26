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
#' @param print_res should the estimation results be printed?
#' @param ... additional settings passed to the function \code{GAfit} employing the genetic algorithm.
#' @details
#'  Because of complexity and multimodality of the log-likelihood function, it's \strong{not guaranteed} that the estimation
#'  algorithm will end up in the global maximum point. It's often expected that most of the estimation rounds will end up in
#'  some local maximum point instead, and therefore a number of estimation rounds is required for reliable results. Because
#'  of the nature of the models, the estimation may fail particularly in the cases where the number of mixture components is
#'  chosen too large. Note that the genetic algorithm is designed to avoid solutions with mixing weights of some regimes
#'  too close to zero at almost all times ("redundant regimes") but the settings can, however, be adjusted (see ?GAfit).
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
#'  \strong{Addiotional Notes about the estimates:}
#'
#'  Sometimes the found MLE is very close to the boundary of the stationarity region some regime, the related variance parameter
#'  is very small, and the associated mixing weights are "spiky". This kind of estimates often maximize the log-likelihood function
#'  for a technical reason that induces by the endogenously determined mixing weights. In such cases, it might be more appropriate
#'  to consider the next-best local maximum point of the log-likelihood function that is well inside the parameter space. Models based
#'  local-only maximum points can be built with the function \code{alt_gsmar} by adjusting the argument \code{which_largest}
#'  accordingly.
#'
#'  Some mixture components of the StMAR model may sometimes get very large estimates for the degrees of freedom parameters. Such parameters
#'  are weakly identified and induce various numerical problems. However, mixture components with large degree of freedom parameters are
#'  similar to the mixture components of the GMAR model. It's hence advisable to further estimate a G-StMAR model by allowing the mixture
#'  components with large degrees of freedom parameter estimates to be GMAR type with the function \code{stmar_to_gstmar}.
#' @return Returns an object of class \code{'gsmar'} defining the estimated GMAR, StMAR or G-StMAR model. The returned object contains
#'   estimated mixing weights, some conditional and unconditional moments, and quantile residuals. Note that the first \code{p}
#'   observations are taken as the initial values, so the mixing weights, conditional moments, and quantile residuals start from
#'   the \code{p+1}:th observation (interpreted as t=1). In addition, the returned object contains the estimates and log-likelihoods
#'   from all of the estimation rounds. See \code{?GSMAR} for the form of the parameter vector, if needed.
#' @section S3 methods:
#'  The following S3 methods are supported for class \code{'gsmar'} objects: \code{print}, \code{summary}, \code{plot},
#'  \code{logLik}, \code{residuals}.
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
#' ## These are long running examples that use parallel computing.
#' ## The below examples take approximately 90 seconds to run.
#'
#' ## Note that the number of estimation rounds (ncalls) is relatively small
#' ## in the below examples. For reliable results, a large number of estimation
#' ## rounds is recommended!
#'
#' # GMAR model
#' fit12 <- fitGSMAR(data=simudata, p=1, M=2, model="GMAR", ncalls=4, seeds=1:4)
#' summary(fit12)
#' plot(fit12)
#' profile_logliks(fit12)
#' diagnostic_plot(fit12)
#'
#' # StMAR model (boundary estimate + large degrees of freedom)
#' fit42t <- fitGSMAR(data=M10Y1Y, p=4, M=2, model="StMAR", ncalls=6, seeds=1:6)
#' summary(fit42t, digits=4) # Four almost-unit roots in the 2nd regime!
#' plot(fit42t) # Spiking mixing weights!
#' fit42t_alt <- alt_gsmar(fit42t, which_largest=2) # The second largest local max
#' summary(fit42t_alt) # Overly large 2nd regime degrees of freedom estimate!
#' fit42gs <- stmar_to_gstmar(fit42t_alt) # Switch to G-StMAR model
#' summary(fit42gs) # Finally, an appropriate model!
#' plot(fit42gs)
#'
#' # Restricted StMAR model
#' fit42r <- fitGSMAR(M10Y1Y, p=4, M=2, model="StMAR", restricted=TRUE,
#'                    ncalls=2, seeds=1:2)
#' fit42r
#'
#' # G-StMAR model with one GMAR type and one StMAR type regime
#' fit42gs <- fitGSMAR(M10Y1Y, p=4, M=c(1, 1), model="G-StMAR",
#'                     ncalls=1, seeds=4)
#' fit42gs
#'
#' # The following three examples demonstrate how to apply linear constraints
#' # to the autoregressive (AR) parameters.
#'
#' # Two-regime GMAR p=2 model with the second AR coeffiecient of
#' # of the second regime contrained to zero.
#' C22 <- list(diag(1, ncol=2, nrow=2), as.matrix(c(1, 0)))
#' fit22c <- fitGSMAR(M10Y1Y, p=2, M=2, constraints=C22, ncalls=1, seeds=6)
#' fit22c
#'
#' # StMAR(3, 1) model with the second order AR coefficient constrained to zero.
#' C31 <- list(matrix(c(1, 0, 0, 0, 0, 1), ncol=2))
#' fit31tc <- fitGSMAR(M10Y1Y, p=3, M=1, model="StMAR", constraints=C31,
#'                     ncalls=1, seeds=1)
#' fit31tc
#'
#' # Such StMAR(3, 2) model that the AR coefficients are restricted to be
#' # the same for both regimes and the second AR coefficients are
#' # constrained to zero.
#' fit32rc <- fitGSMAR(M10Y1Y, p=3, M=2, model="StMAR", restricted=TRUE,
#'                     constraints=matrix(c(1, 0, 0, 0, 0, 1), ncol=2),
#'                     ncalls=1, seeds=1)
#' fit32rc
#' }
#' @export

fitGSMAR <- function(data, p, M, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL, conditional=TRUE,
                     parametrization=c("intercept", "mean"), ncalls=round(10 + 9*log(sum(M))), ncores=min(2, ncalls, parallel::detectCores()),
                     maxit=300, seeds=NULL, print_res=TRUE, ...) {
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

  # We check and warn about deprecated arguments found in dot arguments
  if(!is.null(dot_params$printRes)) {
    warning("The argument 'printRes' is deprecated! Use 'print_res' instead!")
    print_res <- dot_params$printRes
  }

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

  if(print_res) {
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

  if(print_res) {
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

  ### Wrap up ###
  ret <- GSMAR(data=data, p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints,
               conditional=conditional, parametrization=parametrization, calc_qresiduals=TRUE,
               calc_cond_moments=TRUE, calc_std_errors=TRUE)
  ret$all_estimates <- newtonEstimates
  ret$all_logliks <- loks
  ret$which_converged <- converged
  ret$which_round <- bestind # Which estimation round induced the largest log-likelihood?
  warn_ar_roots(ret)
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
#' # Estimate GMAR model with on only 1 iteration in variable metric algorithm
#' fit12 <- fitGSMAR(simudata, p=1, M=2, maxit=1, ncalls=1, seeds=1)
#' fit12
#'
#' # Iterate more since iteration limit was reached
#' fit12 <- iterate_more(fit12)
#' fit12
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

  ret <- GSMAR(data=gsmar$data, p=gsmar$model$p, M=gsmar$model$M, params=res$par, model=gsmar$model$model,
               restricted=gsmar$model$restricted, constraints=gsmar$model$constraints,
               conditional=gsmar$model$conditional, parametrization=gsmar$model$parametrization,
               calc_qresiduals=TRUE, calc_cond_moments=TRUE, calc_std_errors=calc_std_errors, custom_h=custom_h)

  ret$all_estimates <- gsmar$all_estimates
  ret$all_logliks <- gsmar$all_logliks
  ret$which_converged <- gsmar$which_converged
  if(!is.null(gsmar$which_round)) {
    ret$which_round <- gsmar$which_round
    ret$all_estimates[[gsmar$which_round]] <- ret$params
    ret$all_logliks[gsmar$which_round] <- ret$loglik
    ret$which_converged[gsmar$which_round] <- res$convergence == 0
  }
  warn_ar_roots(ret)
  ret
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

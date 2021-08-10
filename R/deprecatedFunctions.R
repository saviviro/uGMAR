#' @title DEPRECATED, USE \code{is_stationary} INSTEAD! Check the stationary condition of specified GMAR, StMAR, or G-StMAR model.
#'
#' @description \code{isStationary} checks the stationarity condition of the specified GMAR, StMAR, or G-StMAR model.
#'  DEPRECATED, USE \code{is_stationary} INSTEAD!
#'
#' @inheritParams loglikelihood
#' @details DEPRECATED, USE \code{is_stationary} INSTEAD!
#' This function falsely returns \code{FALSE} for stationary models when the parameter is extremely close
#'  to the boundary of the stationarity region.
#' @inherit is_stationary_int return references
#' @export

isStationary <- function(p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL) {
  .Deprecated("is_stationary")
  is_stationary(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints)
}


#' @title DEPRECATED, USE \code{diagnostic_plot} INSTEAD! Quantile residual based diagnostic plots for GMAR, StMAR, and G-StMAR models
#'
#' @description \code{diagnosticPlot} plots quantile residual time series, normal QQ-plot, autocorrelation function,
#'  and squared quantile residual autocorrelation function. There is an option to also plot the individual statistics
#'  associated with the quantile residual tests (for autocorrelation and conditional heteroskedasticity) divided by
#'  their approximate standard errors with their approximate 95\% critical bounds (see Kalliovirta 2012, Section 3).
#'  DEPRECATED, USE \code{diagnostic_plot} INSTEAD!
#'
#' @inheritParams diagnostic_plot
#' @details DEPRECATED, USE \code{diagnostic_plot} INSTEAD!
#'
#'  Sometimes the individual statistics are not plotted because it was not (numerically) possible
#'  to calculate all the required statistics. This may suggest that the model is misspecified.
#'
#'  The dashed lines plotted with autocorrelation functions (for quantile residuals and their squares) are
#'  plus-minus \eqn{1.96*T^{-1/2}} where \eqn{T} is the sample size (minus the \eqn{p} initial values for
#'  conditional models).
#' @inherit diagnostic_plot return references
#' @section Suggested packages:
#'   Install the suggested package "gsl" for faster evaluations in the cases of StMAR and G-StMAR models.
#'   For large StMAR and G-StMAR models with large data the calculations to obtain the individual statistics
#'   may take a significantly long time without the package "gsl".
#' @seealso \code{\link{profile_logliks}}, \code{\link{get_foc}}, \code{\link{fitGSMAR}}, \code{\link{cond_moment_plot}}, \code{\link{quantile_residual_tests}},
#'  \code{\link{quantile_residual_plot}}, \code{\link{simulate.gsmar}}, \code{\link{LR_test}}, \code{\link{Wald_test}}
#' @export

diagnosticPlot <- function(gsmar, nlags=20, nsimu=1, plot_indstats=FALSE) {
  .Deprecated("diagnostic_plot")
  diagnostic_plot(gsmar=gsmar, nlags=nlags, nsimu=nsimu, plot_indstats=plot_indstats)
}


#' @title DEPRECATED, USE \code{quantile_residual_plot} INSTEAD! Plot quantile residual time series and histogram
#'
#' @description \code{quantileResidualsPlot} plots quantile residual time series and histogram.
#'  DEPRECATED, USE \code{quantile_residual_plot} INSTEAD!
#'
#' @inheritParams add_data
#' @details DEPRECATED, USE \code{quantile_residual_plot} INSTEAD!
#' @inherit quantile_residual_plot return references
#' @seealso  \code{\link{profile_logliks}}, \code{\link{diagnostic_plot}}, \code{\link{fitGSMAR}}, \code{\link{GSMAR}},
#'  \code{\link{quantile_residual_tests}}, \code{\link{simulate.gsmar}}
#' @export

quantileResidualPlot <- function(gsmar) {
  .Deprecated("quantile_residual_plot")
  quantile_residual_plot(gsmar)
}


#' @title DEPRECATED, USE \code{random_ind} OR \code{smart_ind} INSTEAD!
#'  Create random GMAR, StMAR, or G-StMAR model compatible parameter vector
#'
#' @description \code{randomIndividual} creates a random GMAR, StMAR, or G-StMAR model compatible mean-parametrized parameter vector.
#'    DEPRECATED, USE \code{random_ind} INSTEAD!
#'
#' \code{smartIndividual} creates a random GMAR, StMAR, or G-StMAR model compatible parameter vector close to argument \code{params}.
#'   Sometimes returns exactly the given parameter vector.  DEPRECATED, USE \code{smart_ind} INSTEAD!
#'
#' @inheritParams random_ind_int
#' @param whichRandom deprecated! Use \code{which_random} instead!
#' @param meanscale deprecated! Use \code{mu_scale} instead!
#' @param sigmascale deprecated! Use \code{sigma_scale} instead!
#' @inherit random_ind_int return references
#' @details DEPRECATED, USE \code{random_ind} OR \code{smart_ind} INSTEAD!
#'
#'   These functions can be used, for example, to create initial populations for the genetic algorithm. Mean-parametrization
#'   (instead of intercept terms \eqn{\phi_{m,0}}) is assumed.
#' @export

randomIndividual <- function(p, M, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL, mu_scale, sigma_scale,
                             forcestat=FALSE, meanscale=NULL, sigmascale=NULL) {
  .Deprecated("random_ind")
  if(!is.null(meanscale)) {
    warning("The argument 'meanscale' is deprecated! Use 'mu_scale' instead!")
    mu_scale <- meanscale
  }
  if(!is.null(sigmascale)) {
    warning("The argument 'sigmascale' is deprecated! Use 'sigma_scale' instead!")
    sigma_scale <- sigmascale
  }
  random_ind(p=p, M=M, model=model, restricted=restricted, constraints=constraints,
             mu_scale=mu_scale, sigma_scale=sigma_scale, forcestat=forcestat)
}


#' @rdname randomIndividual
#' @inheritParams smart_ind_int
#' @export

smartIndividual <- function(p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL,
                      mu_scale, sigma_scale, accuracy, which_random=numeric(0), forcestat=FALSE, whichRandom=NULL,
                      meanscale=NULL, sigmascale=NULL) {
  .Deprecated("smart_ind")
  if(!is.null(whichRandom)) {
    warning("The argument 'whichRandom' is deprecated! Use 'whichRandom' instead!")
    which_random <- whichRandom
  }
  if(!is.null(meanscale)) {
    warning("The argument 'meanscale' is deprecated! Use 'mu_scale' instead!")
    mu_scale <- meanscale
  }
  if(!is.null(sigmascale)) {
    warning("The argument 'sigmascale' is deprecated! Use 'sigma_scale' instead!")
    sigma_scale <- sigmascale
  }
  smart_ind(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints,
            mu_scale=mu_scale, sigma_scale=sigma_scale, accuracy=accuracy, which_random=which_random,
            forcestat=forcestat)
}


#' @title DEPRECATED, USE \code{mixing_weights} INSTEAD! Calculate mixing weights of GMAR, StMAR or G-StMAR model
#'
#' @description \code{mixingWeights} calculates the mixing weights of the specified GMAR, StMAR or G-StMAR model
#'  and returns them as a matrix. DEPRECATED, USE \code{mixing_weights} INSTEAD!
#'
#' @inheritParams mixing_weights_int
#' @details DEPRECATED, USE \code{mixing_weights} INSTEAD!
#' @inherit mixing_weights_int return references
#' @export

mixingWeights <- function(data, p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL,
                           parametrization=c("intercept", "mean")) {
  .Deprecated("mixing_weights")
  mixing_weights(data=data, p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints,
                 parametrization=parametrization)
}


#' @title DEPRECATED, USE \code{cond_moments} INSTEAD! Calculate conditional moments of GMAR, StMAR, or G-StMAR model
#'
#' @description \code{condMoments} calculates the regime specific conditional means and variances and total
#'  conditional means and variances of the specified GMAR, StMAR or G-StMAR model.
#'  DEPRECATED, USE \code{cond_moments} INSTEAD!
#'
#' @inheritParams cond_moments
#' @inherit cond_moments return references
#' @export

condMoments <- function(data, p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL,
                         parametrization=c("intercept", "mean"), to_return=c("regime_cmeans", "regime_cvars", "total_cmeans", "total_cvars")) {
  .Deprecated("cond_moments")
  cond_moments(data=data, p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints,
               parametrization=parametrization, to_return=to_return)
}


#' @title DEPRECATED, USE \code{cond_moment_plot} INSTEAD! Conditional mean or variance plot for GMAR, StMAR, and G-StMAR models
#'
#' @description \code{condmomentPlot} plots the one-step in-sample conditional means/variances of the model along with
#'  the time series contained in the model (e.g. the time series the model was fitted to). Also plots
#'  the regimewise conditional means/variances multiplied with the mixing weights.
#'  DEPRECATED, USE \code{cond_moment_plot} INSTEAD!
#'
#' @inheritParams cond_moment_plot
#' @details DEPRECATED, USE \code{cond_moment_plot} INSTEAD!
#' @inherit cond_moment_plot return references
#' @seealso \code{\link{profile_logliks}}, \code{\link{diagnostic_plot}}, \code{\link{fitGSMAR}}, \code{\link{GSMAR}}, \code{\link{quantile_residual_tests}},
#'  \code{\link{quantileResidualPlot}}
#' @export

condmomentPlot <- function(gsmar, which_moment=c("mean", "variance")) {
  .Deprecated("cond_moment_plot")
  cond_moment_plot(gsmar=gsmar, which_moment=which_moment)
}


#' @title DEPRECATED, USE \code{quantile_residuals} INSTEAD! Compute quantile residuals of GMAR, StMAR, or G-StMAR model
#'
#' @description \code{quantileResiduals} computes the quantile residuals of the specified GMAR, StMAR, or G-StMAR model.
#'  DEPRECATED, USE \code{quantile_residuals} INSTEAD!
#'
#' @inheritParams quantile_residuals
#' @details DEPRECATED, USE \code{quantile_residuals} INSTEAD!
#' @inherit quantile_residuals return references
#' @section Suggested packages:
#'   Install the suggested package "gsl" for faster evaluation of the quantile residuals of StMAR and G-StMAR models.
#' @export

quantileResiduals <- function(data, p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE,
                               constraints=NULL, parametrization=c("intercept", "mean")) {
  .Deprecated("quantile_residuals")
  quantile_residuals(data=data, p=p, M=M, params=params, model=model, restricted=restricted,
                     constraints=constraints, parametrization=parametrization)
}


#' @title DEPRECATED, USE \code{quantile_residual_tests} INSTEAD! Quantile residual tests for GMAR, StMAR , and G-StMAR models
#'
#' @description \code{quantileResidualTests} performs quantile residual tests for GMAR, StMAR,
#'  and G-StMAR models, testing normality, autocorrelation, and conditional heteroscedasticity
#'  of the quantile residuals. DEPRECATED, USE \code{quantile_residual_tests} INSTEAD!
#'
#' @inheritParams quantile_residual_tests
#' @param lagsAC deprecated! Use \code{lags_ac} instead.
#' @param lagsCH deprecated! Use \code{lags_ch} instead.
#' @param printRes deprecated! Use \code{print_res} instead.
#' @details DEPRECATED! USE \code{quantile_residual_tests} INSTEAD!
#'
#'   For a correctly specified GSMAR model employing the maximum likelihood estimator, the quantile residuals
#'   are asymptotically independent with standard normal distribution. They can hence be used in a similar
#'   manner to conventional Pearson's residuals. For more details about quantile residual based diagnostics,
#'   and in particular, about the quantile residual tests, see the cited article by \emph{Kalliovirta (2012)}.
#' @section Suggested packages:
#'   Install the suggested package "gsl" for faster evaluations in the cases of StMAR and G-StMAR models.
#'   For large StMAR and G-StMAR models with large data, the evaluations may take significantly long time
#'   without the package "gsl".
#' @inherit quantile_residual_tests return references
#' @seealso \code{\link{profile_logliks}}, \code{\link{fitGSMAR}}, \code{\link{GSMAR}}, \code{\link{diagnostic_plot}},
#'  \code{\link{predict.gsmar}}, \code{\link{get_test_Omega}},
#' @export

quantileResidualTests <- function(gsmar, lags_ac=c(1, 3, 6, 12), lags_ch=lags_ac, nsimu=1, print_res=TRUE,
                                  lagsAC=NULL, lagsCH=NULL, printRes=NULL) {
  if(!is.null(lagsAC)) {
    print("The argument 'lagsAC' is deprecated! Use 'lags_ac' instead!")
    lags_ac <- lagsAC
  }
  if(!is.null(lagsCH)) {
    print("The argument 'lagsCH' is deprecated! Use 'lags_ch' instead!")
    lags_ch <- lagsCH
  }
  if(!is.null(printRes)) {
    print("The argument 'printRes' is deprecated! Use 'print_res' instead!")
    print_res <- printRes
  }

  .Deprecated("quantile_residual_tests")
  quantile_residual_tests(gsmar, lags_ac=lags_ac, lags_ch=lags_ch, nsimu=nsimu, print_res=print_res)
}


#' @import stats
#'
#' @title DEPRECATED, USE \code{simulate.gsmar} INSTEAD! Simulate observations from GMAR, StMAR, and G-StMAR processes
#'
#' @description \code{simulateGSMAR} simulates observations from the specified GMAR, StMAR, or G-StMAR process.
#'  Can be utilized for forecasting future values of the process. DEPRECATED, USE \code{simulate.gsmar} INSTEAD!
#'
#' @inheritParams simulate.gsmar
#' @inheritParams add_data
#' @param nsimu a positive integer specifying how many values (ahead from \code{init_values}) will be simulated.
#' @details DEPRECATED, USE \code{simulate.gsmar} INSTEAD!
#'
#'  The argument \code{ntimes} is intended for forecasting: a GSMAR process can be forecasted by simulating its
#'  possible future values. One can perform a large number of sets of simulations and calculate the sample quantiles from
#'  the simulated values to obtain prediction intervals. See the forecasting example below for a hand-on demonstration.
#' @return If \code{drop==TRUE} and \code{ntimes==1} (default): \code{$sample} and \code{$component} are vectors
#'  and \code{$mixing_weights} is a (\code{nsim}\eqn{xM}) matrix. Otherwise, returns a list with...
#'   \describe{
#'     \item{\code{$sample}}{a size (\code{nsim}\eqn{x}\code{ntimes}) matrix containing the simulated values.}
#'     \item{\code{$component}}{a size (\code{nsim}\eqn{x}\code{ntimes}) matrix containing the information from which
#'      mixture component each value was generated from.}
#'     \item{\code{$mixing_weights}}{a size (\code{nsim}\eqn{xMx}\code{ntimes}) array containing the mixing weights corresponding to the
#'      sample: the dimension \code{[i, , ]} is the time index, the dimension \code{[, i, ]} indicates the regime, and the dimension
#'      \code{[, , i]} indicates the i:th set of simulations.}
#'   }
#' @seealso \code{\link{fitGSMAR}}, \code{\link{GSMAR}}, \code{\link{predict.gsmar}},
#'  \code{\link{add_data}}, \code{\link{cond_moments}}, \code{\link{mixing_weights}}
#' @inherit loglikelihood references
#' @export

simulateGSMAR <- function(object, nsim, init_values=NULL, ntimes=1, drop=TRUE, gsmar=NULL, nsimu=NULL) {
  if(!is.null(gsmar)) {
    print("The argument 'gsmar' is deprecated! Use 'gsmar' instead!")
    object <- gsmar
  }
  if(!is.null(nsimu)) {
    print("The argument 'nsimu' is deprecated! Use 'nsim' instead!")
    nsim <- nsimu
  }
  .Deprecated("simulate.gsmar")

  simulate.gsmar(object=object, nsim=nsim, seed=NULL, init_values=init_values, ntimes=ntimes, drop=drop)
}

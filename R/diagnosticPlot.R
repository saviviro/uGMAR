#' @import stats
#' @import graphics
#' @importFrom grDevices rgb
#'
#' @title Quantile residual based diagnostic plots for GMAR, StMAR, and G-StMAR models
#'
#' @description \code{diagnosticPlot} plots quantile residual time series, normal QQ-plot, autocorrelation function,
#'  and squared quantile residual autocorrelation function. There is an option to also plot the individual statistics
#'  associated with the quantile residual tests (for autocorrelation and conditional heteroskedasticity) divided by
#'  their approximate standard errors with their approximate 95\% critical bounds (see Kalliovirta 2012, Section 3).
#'
#' @inheritParams simulateGSMAR
#' @param nlags a positive integer specifying how many lags should be calculated for the autocorrelation and
#'  conditional heteroscedasticity statistics.
#' @param nsimu a positive integer specifying to how many simulated values from the process the covariance
#'  matrix "Omega" (used to compute the tests) should be based on. Larger number of simulations may result
#'  more reliable tests but takes longer to compute. If smaller than data size, then "Omega" will be based
#'  on the given data. Ignored if \code{plot_indstats==FALSE}.
#' @param plot_indstats set \code{TRUE} if the individual statistics discussed in Kalliovirta (2012) should be
#'  plotted with their approximate 95\% critical bounds (this may take some time).
#' @details Sometimes the individual statistics are not plotted because it was not (numerically) possible
#'  to calculate all the required statistics. This may suggest that the model is misspecified.
#'
#'  The dashed lines plotted with autocorrelation functions (for quantile residuals and their squares) are
#'  plus-minus \eqn{1.96*T^{-1/2}} where \eqn{T} is the sample size (minus the \eqn{p} initial values for
#'  conditional models).
#' @return \code{diagnosticPlot} only plots to a graphical device and does not return anything. Use the
#'  function \code{quantileResidualTests} in order to obtain the individual statistics.
#' @inherit quantileResidualTests references
#' @section Suggested packages:
#'   Install the suggested package "gsl" for faster evaluations in the cases of StMAR and G-StMAR models.
#'   For large StMAR and G-StMAR models with large data the calculations to obtain the individual statistics
#'   may take a significantly long time without the package "gsl".
#' @seealso \code{\link{profile_logliks}}, \code{\link{get_foc}}, \code{\link{fitGSMAR}}, \code{\link{GSMAR}}, \code{\link{quantileResidualTests}},
#'  \code{\link{quantileResidualPlot}}, \code{\link{simulateGSMAR}}
#' @examples
#' \donttest{
#' # GMAR model
#' fit12 <- fitGSMAR(data=logVIX, p=1, M=2, model="GMAR")
#' diagnosticPlot(fit12)
#'
#' # Restricted GMAR model: plot also the individual statistics with
#' # their approximate critical bounds using the given data
#' fit12r <- fitGSMAR(logVIX, 1, 2, model="GMAR", restricted=TRUE)
#' diagnosticPlot(fit12r, nlags=10, nsimu=1, plot_indstats=TRUE)
#'
#' # Non-mixture version of StMAR model
#' fit11t <- fitGSMAR(logVIX, 1, 1, model="StMAR", ncores=1, ncalls=1)
#' diagnosticPlot(fit11t)
#'
#' # G-StMAR model
#' fit12gs <- fitGSMAR(logVIX, 1, M=c(1, 1), model="G-StMAR")
#' diagnosticPlot(fit12gs)
#'
#' # Restricted G-StMAR-model
#' fit12gsr <- fitGSMAR(logVIX, 1, M=c(1, 1), model="G-StMAR",
#'  restricted=TRUE)
#' diagnosticPlot(fit12gsr)
#'
#' # GMAR model as a mixture of AR(2) and AR(1) models
#' constraints <- list(diag(1, ncol=2, nrow=2), as.matrix(c(1, 0)))
#' fit22c <- fitGSMAR(logVIX, 2, 2, constraints=constraints)
#' diagnosticPlot(fit22c)
#'
#' # Such StMAR(3,2) that the AR coefficients are restricted to be
#' # the same for both regimes and that the second AR coefficients are
#' # constrained to zero.
#' fit32rc <- fitGSMAR(logVIX, 3, 2, model="StMAR", restricted=TRUE,
#'  constraints=matrix(c(1, 0, 0, 0, 0, 1), ncol=2))
#' diagnosticPlot(fit32rc)
#' }
#' @export

diagnosticPlot <- function(gsmar, nlags=20, nsimu=1, plot_indstats=FALSE) {
  if(!all_pos_ints(c(nlags, nsimu))) stop("The arguments 'nlags' and 'nsimu' have to be a strictly positive integers")
  check_gsmar(gsmar)
  check_data(gsmar)
  nsimu <- max(nsimu, length(data))
  data <- gsmar$data
  n_obs <- ifelse(gsmar$model$conditional, length(data) - gsmar$model$p, length(data))

  if(is.null(gsmar$quantile_residuals)) {
    qresiduals <- quantileResiduals_int(data=data, p=gsmar$model$p, M=gsmar$model$M, params=gsmar$params,
                                        model=gsmar$model$mode, restricted=gsmar$model$restricted,
                                        constraints=gsmar$model$constraints, parametrization=gsmar$model$parametrization)
  } else {
    qresiduals <- gsmar$quantile_residuals
  }
  if(anyNA(qresiduals)) {
    n_na <- sum(is.na(qresiduals))
    qresiduals <- qresiduals[!is.na(qresiduals)]
    warning(paste(n_na, "missing values removed from quantile residuals. Check the parameter estimates for possible problems (border of the prm space, large dfs, etc)?"))
  }
  old_par <- par(no.readonly=TRUE) # Save old settings
  on.exit(par(old_par)) # Restore the settings before quitting
  if(plot_indstats) {
    par(mfrow=c(3, 2), mar=c(2.1, 2.1, 2.1, 0.8))
  } else {
    par(mfrow=c(2, 2), mar=c(2.1, 2.1, 2.1, 0.8))
  }

  # Plot quantile residuals time series and qq-plot
  yaxt1 <- round(min(qresiduals))
  yaxt2 <- round(max(qresiduals))
  plot(qresiduals, yaxt="n", type="l", col=rgb(0, 0, 0, 1), ylab="", xlab="", main="Quantile residuals")
  axis(2, at=yaxt1:yaxt2, labels=yaxt1:yaxt2)
  abline(h=0, col=rgb(1, 0, 0, 0.3), lwd=2)
  qqnorm(qresiduals, yaxt="n", ylab="", xlab="", main="Normal QQ-plot")
  axis(2, at=yaxt1:yaxt2, labels=FALSE)
  qqline(qresiduals, col=rgb(1, 0, 0, 0.8))

  # Plot autocorrelation function of quantile residuals and their squares
  qr_acf <- as.matrix(acf(qresiduals, lag.max=nlags, type="correlation", plot=FALSE)$acf[2:(nlags + 1)])
  qrsquare_acf <- as.matrix(acf(qresiduals^2, lag.max=nlags, type="correlation", plot=FALSE)$acf[2:(nlags + 1)])
  yaxt0 <- max(0.3, round(max(c(qr_acf, qrsquare_acf) + 0.05, abs(c(qr_acf, qrsquare_acf)) + 0.05), 1))
  ticks1 <- round(seq(-yaxt0, yaxt0, by=0.1), 1)
  ticks2 = seq(-yaxt0, yaxt0, by=0.05)

  plot_qr_acf <- function(vals_to_plot, main) {
    plot(0, 0, type="n", yaxt="n", xlim=c(0, nlags + 0.1), ylim=c(-yaxt0, yaxt0), xlab="", ylab="", main=main)
    axis(2, at=ticks1, labels=ticks1)
    abline(h=ticks2, col=rgb(0.1, 0.1, 0.1, 0.2))
    abline(h=0)
    abline(v=seq(0, nlags, by=5), col=rgb(0.1, 0.1, 0.1, 0.2))
    segments(x0=1:nlags, y0=rep(0, nlags), x1=1:nlags, y1=vals_to_plot)
    points(1:nlags, vals_to_plot, pch=20, col="blue")
    abline(h=c(-1.96*n_obs^{-1/2}, 1.96*n_obs^{-1/2}), col=rgb(0, 0, 1, 0.8), lty=2)
  }

  plot_qr_acf(qr_acf, main="Qres ACF")
  plot_qr_acf(qrsquare_acf, main="Qres^2 ACF")

  if(plot_indstats) {
    # Obtain tests statistics
    qrtest <- quantileResidualTests(gsmar, lagsAC=1:nlags, lagsCH=1:nlags, nsimu=nsimu, printRes=FALSE)

    plot_inds <- function(which_ones) { # ac_res or ch_res
      res <- qrtest[[which(names(qrtest) == which_ones)]]
      inds_normalized <- res$indStat/res$stdError
      plot(0, 0, type="n", yaxt="n", xlim=c(0, nlags + 0.1), ylim=c(-3, 3), xlab="", ylab="",
           main=ifelse(which_ones == "ac_res", "Autocovariance statistics", "Cond. h.sked. statistics"))
      abline(h=c(-3:3), col=rgb(0.1, 0.1, 0.1, 0.2))
      abline(h=0, lty=1)
      abline(v=seq(0, nlags, by=5), col=rgb(0.1, 0.1, 0.1, 0.2))
      abline(h=1.96, col=rgb(0, 0, 1, 0.8), lty=2)
      abline(h=-1.96, col=rgb(0, 0, 1, 0.8), lty=2)
      segments(x0=1:nlags, y0=rep(0, nlags), x1=1:nlags, y1=inds_normalized)
      points(1:nlags, inds_normalized, pch=20, col="blue")
      axis(2, at=-3:3, labels=-3:3)
    }

    plot_inds("ac_res")
    plot_inds("ch_res")
  }
}


#' @import stats
#' @import graphics
#' @importFrom grDevices rgb
#'
#' @title Plot quantile residual time series and histogram
#'
#' @description \code{quantileResidualsPlot} plots quantile residual time series and histogram.
#'
#' @inheritParams simulateGSMAR
#' @return  Only plots to a graphical device and doesn't return anything.
#' @inherit quantileResiduals references
#' @seealso  \code{\link{profile_logliks}}, \code{\link{diagnosticPlot}}, \code{\link{fitGSMAR}}, \code{\link{GSMAR}},
#'  \code{\link{quantileResidualTests}}, \code{\link{simulateGSMAR}}
#' @examples
#' \donttest{
#' # GMAR model
#' fit12 <- fitGSMAR(data=logVIX, p=1, M=2, model="GMAR")
#' quantileResidualPlot(fit12)
#'
#' # Non-mixture version of StMAR model
#' fit11t <- fitGSMAR(logVIX, 1, 1, model="StMAR", ncores=1, ncalls=1)
#' quantileResidualPlot(fit11t)
#'
#' # Restricted G-StMAR-model
#' fit12gsr <- fitGSMAR(logVIX, 1, M=c(1, 1), model="G-StMAR",
#'  restricted=TRUE)
#' quantileResidualPlot(fit12gsr)
#'
#' # Such StMAR(3,2) that the AR coefficients are restricted to be
#' # the same for both regimes and that the second AR coefficients are
#' # constrained to zero.
#' fit32rc <- fitGSMAR(logVIX, 3, 2, model="StMAR", restricted=TRUE,
#'  constraints=matrix(c(1, 0, 0, 0, 0, 1), ncol=2))
#' quantileResidualPlot(fit32rc)
#' }
#' @export

quantileResidualPlot <- function(gsmar) {
  check_gsmar(gsmar)
  check_data(gsmar)
  if(is.null(gsmar$quantile_residuals)) {
    qresiduals <- quantileResiduals_int(data=gsmar$data, p=gsmar$model$p, M=gsmar$model$M, params=gsmar$params,
                                        model=gsmar$model$mode, restricted=gsmar$model$restricted,
                                        constraints=gsmar$model$constraints, parametrization=gsmar$model$parametrization)
  } else {
    qresiduals <- gsmar$quantile_residuals
  }
  if(anyNA(qresiduals)) {
    n_na <- sum(is.na(qresiduals))
    qresiduals <- qresiduals[!is.na(qresiduals)]
    warning(paste(n_na, "missing values removed from quantile residuals. Check the parameter estimates for possible problems (border of the prm space, large dfs, etc)?"))
  }
  old_par <- par(no.readonly = TRUE) # Save old settings
  on.exit(par(old_par)) # Restore the settings before quitting
  par(mfrow=c(2, 1), mar=c(2.6, 2.6, 2.1, 1.6))
  plot(qresiduals, type="l", ylab="", xlab="", main="Quantile residuals")
  abline(h=c(-1.96, 0, 1.96), lty=2, col="red")
  hs <- hist(qresiduals, breaks="Scott", probability=TRUE, col="skyblue", plot=TRUE, ylim=c(0, 0.5))
  xc <- seq(from=min(hs$breaks), to=max(hs$breaks), length.out=500)
  lines(x=xc, y=dnorm(xc), lty=2, col="red")
}


#' @import stats
#' @import graphics
#' @importFrom grDevices rgb
#'
#' @title Plot profile log-likehoods around the estimates
#'
#' @description \code{profile_logliks} plots profile log-likelihoods around the estimates.
#'
#' @inheritParams simulateGSMAR
#' @param scale a numeric scalar specifying the interval plotted for each estimate: the estimate plus-minus \code{abs(scale*estimate)}.
#' @param nrows how many rows should be in the plot-matrix? The default is \code{max(ceiling(log2(nparams) - 1), 1)}.
#' @param ncols how many columns should be in the plot-matrix? The default is \code{ceiling(nparams/nrows)}.
#'   Note that \code{nrows*ncols} should not be smaller than the number of parameters.
#' @param precission at how many points should each profile log-likelihood be evaluated at?
#' @details The red vertical line points the estimate.
#'
#'   Be aware that the profile log-likelihood function is subject to a numerical error due to limited float-point
#'   precission when considering extremely large parameter values, say, overly large degrees freedom estimates (see the
#'   related example below).
#' @return  Only plots to a graphical device and doesn't return anything.
#' @inherit loglikelihood references
#' @seealso  \code{\link{quantileResidualPlot}}, \code{\link{diagnosticPlot}}, \code{\link{fitGSMAR}}, \code{\link{GSMAR}},
#'  \code{\link{quantileResidualTests}}, \code{\link{simulateGSMAR}}
#' @examples
#' \donttest{
#' # GMAR model
#' fit12 <- fitGSMAR(data=logVIX, p=1, M=2, model="GMAR")
#' profile_logliks(fit12)
#'
#' # Non-mixture version of StMAR model
#' fit11t <- fitGSMAR(logVIX, 1, 1, model="StMAR", ncores=1, ncalls=1)
#' profile_logliks(fit11t)
#'
#' # Restricted G-StMAR-model
#' fit12gsr <- fitGSMAR(logVIX, 1, M=c(1, 1), model="G-StMAR",
#'  restricted=TRUE)
#' profile_logliks(fit12gsr)
#'
#' # Extremely large degrees of freedom numerical error demonstration
#' fit12t <- fitGSMAR(logVIX, 1, 2, model="StMAR", ncores=1,
#'  ncalls=1, seeds=1)
#' profile_logliks(fit12t, scale=0.00001)
#' # See the last figure? Surface of the profile log-likelihood function
#' # should be flat around that large degrees of freedom!
#' }
#' @export

profile_logliks <- function(gsmar, scale=0.02, nrows, ncols, precission=200) {
  check_gsmar(gsmar)
  check_data(gsmar)
  p <- gsmar$model$p
  M_orig <- gsmar$model$M
  M <- sum(M_orig)
  params <- gsmar$params
  npars <- length(params)
  constraints <- gsmar$model$constraints
  restricted <- gsmar$model$restricted
  parametrization <- gsmar$model$parametrization
  if(missing(nrows)) nrows <- max(ceiling(log2(npars) - 1), 1)
  if(missing(ncols)) ncols <- ceiling(npars/nrows)
  stopifnot(all_pos_ints(c(nrows, ncols)) && nrows*ncols >= npars)

  old_par <- par(no.readonly = TRUE) # Save old settings
  on.exit(par(old_par)) # Restore the settings before quitting
  par(mar=c(2.1, 2.1, 1.6, 1.1), mfrow=c(nrows, ncols))

  # In order to get the labels right, we first determine which indeces in params
  # correspond to which parameters: different procedure for restricted model.
  if(restricted == FALSE) {
    if(!is.null(constraints)) {
      all_q <- vapply(1:length(constraints), function(m) ncol(constraints[[m]]), numeric(1)) # q = number of estimated AR parameters in a regime
    } else {
      all_q <- rep(p, M) # Consider non-constrained models as special cases of constrainted models
    }
    cum_q <- c(0, cumsum(all_q + 2)) # Indeces in parameter vector after which the regime changes (before alphas)
  } else { # restricted == TRUE
    cum_q <- rep(0, times=M + 1) # For techinal reasons so that "first arguments" below are always false
    q <- ifelse(is.null(constraints), p, ncol(constraints)) # Constraints is a singe matrix
  }

  for(i1 in seq_len(npars)) { # Go though the parameters
    pars <- params
    range <- abs(scale*pars[i1])
    vals <- seq(from=pars[i1] - range, to=pars[i1] + range, length.out=precission) # Loglik to be evaluated at these values of the parameter considered
    logliks <- vapply(vals, function(val) {
      new_pars <- pars
      new_pars[i1] <- val # Change the single parameter value
      loglikelihood_int(data=gsmar$data, p=p, M=M_orig, params=new_pars, model=gsmar$model$model, restricted=restricted,
                        constraints=constraints, conditional=gsmar$model$conditional, parametrization=parametrization,
                        boundaries=TRUE, checks=FALSE, minval=NA)
    }, numeric(1))

    if(sum(is.na(logliks)) == precission) stop("Profile log-likelihood function is too peaky - increase precission (also estimates might be peculiar)")

    # In order to get the labels right, we first determine which parameter is in question.
    # For readability of the code, we do the cases of restricted and unrestricted models
    # complitely separately at the expense of some dublicate code.
    if(restricted == FALSE) {
      if(i1 <= max(cum_q)) { # phi and sigma^2 parameters first
        m <- sum(i1 > cum_q) # Which regime are we considering
        if(i1 == cum_q[m + 1]) { # sigma^2
          main <- substitute(sigma[foo]^2, list(foo=m))
        } else if(i1 <= max(cum_q)) { # phi_{m,0},...,phi_{m,p}
          if(i1 == cum_q[m] + 1) { # phi_{m,0}
            mylist <- list(foo=paste0(m, ",", 0))
            if(parametrization == "intercept") {
              main <- substitute(phi[foo], mylist)
            } else { # Different label for mean parametrization
              main <- substitute(mu[foo], mylist)
            }
          } else {  # phi_{m,1},...,phi_{m,p}
            p0 <- i1 - cum_q[m] - 1 # p0 = 1,...,p; minus 1 from phi_0
            mylist  <- list(foo=paste0(m, ",", p0))
            if(!is.null(constraints)) {
              # The AR parameters are not generally the same as phi-parameters with linear constraints
              mylist$phi <- "AR"
            }
            main <- substitute(phi[foo], mylist)
          }
        }
      } else { # alphas and degrees of freedom if any
        if(M == 1) {
          # No alphas, so if we are here there must be dfs parameters
          m <- 1
          main <- substitute(nu[foo], list(foo=m))
        } else {
          if(i1 <= max(cum_q) + M - 1) { # The alphas first
            m <- i1 - max(cum_q)
            main <- substitute(alpha[foo], list(foo=m))
          } else { # Finally, degrees of freedom
            m <- i1 - max(cum_q) - (M - 1)
            main <- substitute(nu[foo], list(foo=m))
          }
        }
      }
    } else { ## restricted == TRUE, taken care separately for readability of the code
        if(i1 <= M + q) { # phi_{m,0},...,phi_{m,p}
          if(i1 <= M) { # phi_{m,0}
            m <- i1
            mylist <- list(foo=paste0(m, ",", 0))
            if(parametrization == "intercept") {
              main <- substitute(phi[foo], mylist)
            } else { # Different label for mean parametrization
              main <- substitute(mu[foo], mylist)
            }
          } else { # phi_{m,1},...,phi_{m,p}
            m <- "m"
            p0 <- i1 - M
            if(is.null(constraints)) {
              mylist  <- list(foo=paste0(m, ",", p0))
            } else { # The AR parameters are not generally the same as phi-parameters with linear constraints
              mylist  <- list(phi="AR", foo=paste0(m, ",", p0))
            }
            main <- substitute(phi[foo], mylist)
          }
        } else if(i1 <= 2*M + q) { # sigma^2
          m <- i1 - M - q
          main <- substitute(sigma[foo]^2, list(foo=m))
        } else { # alphas and degrees of freedom
          if(M == 1) { # No alphas if this is true
            m <- 1
            main <- substitute(nu[foo]^2, list(foo=m))
          } else {
            if(i1 <= 3*M + q - 1) { # The alphas first
              m <- i1 - 2*M - q
              main <- substitute(alpha[foo], list(foo=m))
            } else { # Finally, degrees of freedom
              m <- i1 - (3*M + q - 1)
              main <- substitute(nu[foo], list(foo=m))
            }
          }
        }
    }
    plot(x=vals, y=logliks, type="l", main=main)
    abline(v=pars[i1], col="red") # Points the estimate
  }
}

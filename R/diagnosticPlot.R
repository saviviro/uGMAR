#' @import stats
#' @import graphics
#' @importFrom grDevices rgb
#'
#' @title Quantile residual based diagnostic plots for GMAR, StMAR and G-StMAR models
#'
#' @description \code{diagnosticPlot} plots quantile residual time series, normal QQ-plot, autocorrelation function
#'  and squared quantile residual autocorrelation function. There is an option to also plot the individual statistics
#'  associated with the quantile residual tests (for autocorrelation and conditional heteroskedasticity) divided by
#'  their approximate standard errors with their approximate 95\% critical bounds.
#'
#' @inheritParams simulateGSMAR
#' @param nlags a positive integer specifying how many lags should be calculated for the autocorrelation and
#'  conditional heteroscedasticity statistics.
#' @param nsimu a positive integer specifying to how many simulated values from the process the covariance
#'  matrix "Omega" (used to compute the tests) should be based on. Larger number of simulations may result
#'  more reliable tests. If smaller than data size, then it will be based on the given data.
#'  Ignored if \code{plot_indstats==FALSE}.
#' @param plot_indstats set \code{TRUE} if the individual statistics discussed in Kalliovirta (2012) should be
#'  plotted with their approximate 95\% critical bounds (this may take some time).
#' @details Sometimes the individual statistics are not plotted because it's not (numerically) possible for
#'  to calculate all the necessary estimates required. This may suggest that the model is misspecified.
#'
#'  The dashed lines plotted with autocorrelation functions (for quantile residuals and their squares) are
#'  plus-minus \eqn{1.96*T^{-1/2}}.
#' @return \code{diagnosticPlot} only plots to a graphical device and doesn't return anything. Use the
#'  function \code{quantileResidualTests} in order to obtain the individual statistics.
#' @inherit quantileResidualTests references
#' @section Suggested packages:
#'   Install the suggested package "gsl" for faster evaluations in the cases of StMAR and G-StMAR models.
#'   For large StMAR and G-StMAR models with large data the calculations to obtain the individual statistics
#'   may take a significantly long time without the package "gsl".
#' @seealso \code{\link{fitGSMAR}}, \code{\link{GSMAR}}, \code{\link{quantileResidualTests}}, \code{\link{quantileResidualPlot}}, \code{\link{simulateGSMAR}}
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

diagnosticPlot <- function(gsmar, nlags=20, nsimu=2000, plot_indstats=FALSE) {
  if(!all_pos_ints(c(nlags, nsimu))) stop("The arguments nlags and nsimu have to be a strictly positive integers")
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
  old_par <- par(no.readonly = TRUE) # Save old settings
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
#' @title Ploy quantile residual time series and kernel density
#'
#' @description \code{quantileResidualsPlot} plots quantile residual time series and histogram.
#'
#' @inheritParams simulateGSMAR
#' @return  Only plots to a graphical device and doesn't return anything.
#' @inherit quantileResiduals
#' @seealso \code{\link{diagnosticPlot}}, \code{\link{fitGSMAR}}, \code{\link{GSMAR}}, \code{\link{quantileResidualTests}}, \code{\link{simulateGSMAR}}
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

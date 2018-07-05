#' @import stats
#' @import graphics
#' @importFrom grDevices rgb
#'
#' @title Quantile residual based diagnostic plots for GMAR, StMAR and G-StMAR models
#'
#' @description \code{plotGMAR} plots quantile residual time series, normal QQ-plot, autocorrelation function and squared
#'   quantile residual autocorrelation function. There is an option to also plot the individual statistics associated with the quantile residual
#'   tests (for autocorrelation and conditional heteroskedasticity) divided by their approximate standard errors,
#'   with their approximate 95\% critical bounds.
#'
#' @inheritParams loglikelihood
#' @param nlags an (optional) positive integer specifying how many lags should be calculated for the autocorrelation and conditional heteroscedasticity statistics. Default is 20.
#' @param nsimu an (optional) positive integer specifying to how many process's simulated values the covariance matrix "Omega" should be based on.
#'  If smaller than data size, then it will be based on the given data. Ignored if \code{approxBounds==FALSE}. Default is 2000.
#' @param approxBounds set \code{TRUE} if the Kalliovirta's (2012) individual statistics and their approximate 95\% critical bounds should be plotted (this may take some time). Default is FALSE.
#' @details Sometimes the individual statistics are not plotted because it was not (numerically) possible for the code to
#'   calculate all the necessary estimates required. This may suggest that the model is misspecified.
#'
#'   The dashed lines plotted with autocorrelation functions (for quantile residuals and their squares) are plus-minus \eqn{1.96*T^{-1/2}}.
#' @return \code{plotGMAR} only plots to a graphical device and doesn't return anything. Use the function \code{quantileResidualTests} if you wish to obtain the statistics.
#' @inherit quantileResidualTests references
#' @section Suggested packages:
#'   Install the suggested package "gsl" for faster evaluations in the cases of StMAR and G-StMAR models.
#'   For large StMAR and G-StMAR models with large data the evaluations for approximating the critical bounds may take
#'   significantly long time without the package "gsl".
#' @examples
#' \donttest{
#' # GMAR model
#' params13 <- c(1.4, 0.88, 0.26, 2.46, 0.82, 0.74, 5.0, 0.68, 5.2, 0.72, 0.2)
#' plotGMAR(VIX, 1, 3, params13)
#'
#' # Restricted GMAR model: plot also the individual statistics with their approximate
#' # critical bounds using the given data
#' params12r <- c(1.4, 1.8, 0.88, 0.29, 3.18, 0.84)
#' plotGMAR(VIX, 1, 2, params12r, restricted=TRUE, nsimu=1, approxBounds=TRUE)
#'
#' # StMAR model
#' params12t <- c(1.38, 0.88, 0.27, 3.8, 0.74, 3.15, 0.8, 100, 3.6)
#' plotGMAR(VIX, 1, 2, params12t, StMAR=TRUE)
#'
#' # G-StMAR model (similar to the StMAR model above)
#' params12gs <- c(1.38, 0.88, 0.27, 3.8, 0.74, 3.15, 0.8, 3.6)
#' plotGMAR(VIX, 1, c(1,1), params12gs, GStMAR=TRUE)
#'
#' # Restricted G-StMAR-model
#' params13gsr <- c(1.3, 1, 1.4, 0.8, 0.4, 2, 0.2, 0.25, 0.15, 20)
#' plotGMAR(VIX, 1, c(2, 1), params13gsr, GStMAR=TRUE, restricted=TRUE)
#'
#' # GMAR model as a mixture of AR(2) and AR(1) models
#' R <- list(diag(1, ncol=2, nrow=2), as.matrix(c(1, 0)))
#' params22c <- c(1.2, 0.85, 0.04, 0.3, 3.3, 0.77, 2.8, 0.77)
#' plotGMAR(VIX, 2, 2, params22c, constraints=TRUE, R=R)
#'
#' # Such StMAR(3,2) that the AR coefficients are restricted to be
#' # the same for both regimes and that the second AR coefficients are
#' # constrained to zero.
#' params32trc <- c(2.2, 1.8, 0.88, -0.03, 2.4, 0.27, 0.40, 3.9, 1000)
#' plotGMAR(VIX, 3, 2, params32trc, StMAR=TRUE, restricted=TRUE,
#'          constraints=TRUE, R=matrix(c(1, 0, 0, 0, 0, 1), ncol=2))
#' }
#' @export

plotGMAR <- function(data, p, M, params, StMAR=FALSE, GStMAR=FALSE, restricted=FALSE, constraints=FALSE, R, nlags=20, nsimu=2000, approxBounds=FALSE) {

  checkLogicals(StMAR=StMAR, GStMAR=GStMAR)
  checkPM(p, M, GStMAR=GStMAR)
  M_orig = M
  if(GStMAR==TRUE) {
    M1 = M[1]
    M2 = M[2]
    M = sum(M)
  }
  if(length(params)!=nParams(p=p, M=M_orig, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted, constraints=constraints, R=R)) {
    stop("The parameter vector has wrong dimension")
  }
  if(nlags<1) {
    stop("The number of lags has to be equal or greater than one")
  }
  data = checkAndCorrectData(data, p)
  nsimu = max(nsimu, length(data))
  qresiduals = quantileResiduals_int(data, p, M_orig, params, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted, constraints=constraints, R=R)
  old_par = par(no.readonly=T) # Save old settings
  on.exit(par(old_par)) # Restore the settings before quitting
  if(approxBounds==TRUE) {
    par(mfrow=c(3, 2), mar=c(2.1, 2.1, 2.1, 0.8)) # Set new temporary settings
  } else {
    par(mfrow=c(2, 2), mar=c(2.1, 2.1, 2.1, 0.8))
  }

  # Plot quantile residuals time series and qq-plot
  yaxt1 = round(min(qresiduals)); yaxt2=round(max(qresiduals))
  plot(qresiduals, yaxt="n", type="l", col=rgb(0, 0, 0, 1), ylab="", xlab="", main="Quantile residuals")
  axis(2, at=yaxt1:yaxt2, labels=yaxt1:yaxt2)
  abline(h=0, col=rgb(1, 0, 0, 0.3), lwd=2)
  qqnorm(qresiduals, yaxt="n", ylab="", xlab="", main="Normal QQ-plot")
  axis(2, at=yaxt1:yaxt2, labels=FALSE)
  qqline(qresiduals, col=rgb(1, 0, 0, 0.8))

  # Plot autocorrelation function of quantile residuals
  qr_acf = as.matrix(acf(qresiduals, lag.max=nlags, type="correlation", plot=FALSE)$acf[2:(nlags+1)])
  qrsquare_acf = as.matrix(acf(qresiduals^2, lag.max=nlags, type="correlation", plot=FALSE)$acf[2:(nlags+1)])
  yaxt0 = max(0.3, round(max(c(qr_acf, qrsquare_acf)+0.05, abs(c(qr_acf, qrsquare_acf))+0.05), 1))
  ticks1 = round(seq(-yaxt0, yaxt0, by=0.1), 1); ticks2 = seq(-yaxt0, yaxt0, by=0.05)
  plot(0, 0, type="n", yaxt="n", xlim=c(0, nlags+1), ylim=c(-yaxt0, yaxt0), xlab="", ylab="", main="Quantile residual acf")
  axis(2, at=ticks1, labels=ticks1)
  abline(h=ticks2, col=rgb(0.1, 0.1, 0.1, 0.2)); abline(h=0); abline(v=seq(0, nlags, by=5), col=rgb(0.1, 0.1, 0.1, 0.2))
  segments(x0=1:nlags, y0=rep(0, nlags), x1=1:nlags, y1=qr_acf)
  points(1:nlags, qr_acf, pch=20, col="blue")
  abline(h=c(-1.96*length(data)^{-1/2}, 1.96*length(data)^{-1/2}), col=rgb(0, 0, 1, 0.8), lty=2)

  # Plot autocorrelation function of squared quantile residuals
  plot(0, 0, type="n", yaxt="n", xlim=c(0, nlags+1), ylim=c(-yaxt0, yaxt0), xlab="", ylab="", main="Squared quant. res. acf")
  axis(2, at=ticks1, labels=FALSE)
  abline(h=ticks2, col=rgb(0.1, 0.1, 0.1, 0.2)); abline(h=0); abline(v=seq(0, nlags, by=5), col=rgb(0.1, 0.1, 0.1, 0.2))
  segments(x0=1:nlags, y0=rep(0, nlags), x1=1:nlags, y1=qrsquare_acf)
  points(1:nlags, qrsquare_acf, pch=20, col="blue")
  abline(h=c(-1.96*length(data)^{-1/2}, 1.96*length(data)^{-1/2}), col=rgb(0, 0, 1, 0.8), lty=2)

  if(approxBounds==TRUE) {
    # Obtain tests statistics
    testResults = quantileResidualTests(data, p, M_orig, params, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted, constraints=constraints,
                                        R=R, lagsAC=1:nlags, lagsCH=1:nlags, nsimu=nsimu, printRes=FALSE)

    # Plot autocorrelation statistics
    ck_normalized = testResults$autocorrelation$indStat/testResults$autocorrelation$stdError
    plot(0, 0, type="n", yaxt="n", xlim=c(0, nlags+1), ylim=c(-3, 3), xlab="", ylab="", main="Autocovariance statistics")
    abline(h=0, lty=1); abline(h=1.96, col=rgb(0, 0, 1, 0.8), lty=2); abline(h=-1.96, col=rgb(0, 0, 1, 0.8), lty=2)
    segments(x0=1:nlags, y0=rep(0,nlags), x1=1:nlags, y1=ck_normalized)
    points(1:nlags, ck_normalized, pch=20, col="blue")
    axis(2, at=-3:3, labels=-3:3)

    # Plot conditional heteroscedasticity statistics
    dk_normalized = testResults$cond.heteroscedasticity$indStat/testResults$cond.heteroscedasticity$stdError
    plot(0, 0, type="n", yaxt="n", xlim=c(0, nlags+1), ylim=c(-3, 3), xlab="", ylab="", main="Cond. h.sked. statistics")
    abline(h=0, lty=1); abline(h=1.96, col=rgb(0, 0, 1, 0.8), lty=2); abline(h=-1.96, col=rgb(0, 0, 1, 0.8), lty=2)
    segments(x0=1:nlags, y0=rep(0,nlags), x1=1:nlags, y1=dk_normalized)
    points(1:nlags, dk_normalized, pch=20, col="blue")
    axis(2, at=-3:3, labels=FALSE)
  }
}

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
#' @seealso \code{\link{fitGSMAR}}, \code{\link{GSMAR}}, \code{\link{quantileResidualTests}}, \code{\link{predict.gsmar}}
#' @examples
#' \donttest{
#' # GMAR model
#' params13 <- c(1.4, 0.88, 0.26, 2.46, 0.82, 0.74, 5.0, 0.68, 5.2, 0.72, 0.2)
#' gmar13 <- GSMAR(data=VIX, p=1, M=3, params=params13, model="GMAR")
#' diagnosticPlot(gmar13)
#'
#' # Restricted GMAR model: plot also the individual statistics with
#' # their approximate critical bounds using the given data
#' params12r <- c(1.4, 1.8, 0.88, 0.29, 3.18, 0.84)
#' gmar12r <- GSMAR(data=VIX, p=1, M=2, params=params12r, model="GMAR",
#'  restricted=TRUE)
#' diagnosticPlot(gmar12r, nlags=10, nsimu=1, plot_indstats=TRUE)
#'
#' # StMAR model
#' params12t <- c(1.38, 0.88, 0.27, 3.8, 0.74, 3.15, 0.8, 100, 3.6)
#' stmar12 <- GSMAR(data=VIX, p=1, M=2, params=params12t, model="StMAR")
#' diagnosticPlot(stmar12)
#'
#' # G-StMAR model (similar to the StMAR model above)
#' params12gs <- c(1.38, 0.88, 0.27, 3.8, 0.74, 3.15, 0.8, 3.6)
#' gstmar12 <- GSMAR(data=VIX, p=1, M=c(1, 1), params=params12gs,
#'  model="G-StMAR")
#' diagnosticPlot(gstmar12)
#'
#' # Restricted G-StMAR-model
#' params13gsr <- c(1.3, 1, 1.4, 0.8, 0.4, 2, 0.2, 0.25, 0.15, 20)
#' gstmar13r <- GSMAR(data=VIX, p=1, M=c(2, 1), params=params13gsr,
#'  model="G-StMAR", restricted=TRUE)
#' diagnosticPlot(gstmar13r)
#'
#' # GMAR model as a mixture of AR(2) and AR(1) models
#' constraints <- list(diag(1, ncol=2, nrow=2), as.matrix(c(1, 0)))
#' params22c <- c(1.2, 0.85, 0.04, 0.3, 3.3, 0.77, 2.8, 0.77)
#' gmar22c <- GSMAR(data=VIX, p=2, M=2, params=params22c,
#'  model="GMAR", constraints=constraints)
#' diagnosticPlot(gmar22c)
#'
#' # Such StMAR(3,2) that the AR coefficients are restricted to be
#' # the same for both regimes and that the second AR coefficients are
#' # constrained to zero.
#' params32trc <- c(2.2, 1.8, 0.88, -0.03, 2.4, 0.27, 0.40, 3.9, 1000)
#' stmar32rc <- GSMAR(data=VIX, p=3, M=2, params=params32trc, model="StMAR",
#'  restricted=TRUE, constraints=matrix(c(1, 0, 0, 0, 0, 1), ncol=2))
#' diagnosticPlot(stmar32rc)
#' }
#' @export

diagnosticPlot <- function(gsmar, nlags=20, nsimu=2000, plot_indstats=FALSE) {
  if(!all_pos_ints(c(nlags, nsimu))) stop("The arguments nlags and nsimu have to be a strictly positive integers")
  check_gsmar(gsmar)
  check_data(gsmar)
  data <- gsmar$data
  p <- gsmar$model$p
  M <- gsmar$model$M
  params <- gsmar$params
  model <- gsmar$model$model
  restricted <- gsmar$model$restricted
  constraints <- gsmar$model$constraints
  parametrization <- gsmar$model$parametrization
  nsimu <- max(nsimu, length(data))
  qresiduals <- quantileResiduals_int(data=data, p=p, M=M, params=params, model=model, restricted=restricted,
                                      constraints=constraints, parametrization=parametrization)
  old_par <- par(no.readonly = TRUE) # Save old settings
  on.exit(par(old_par)) # Restore the settings before quitting
  if(plot_indstats == TRUE) {
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

  plot_qr_acf <- function(vals_to_plot) {
    plot(0, 0, type="n", yaxt="n", xlim=c(0, nlags + 0.1), ylim=c(-yaxt0, yaxt0), xlab="", ylab="", main="Quantile residual acf")
    axis(2, at=ticks1, labels=ticks1)
    abline(h=ticks2, col=rgb(0.1, 0.1, 0.1, 0.2))
    abline(h=0)
    abline(v=seq(0, nlags, by=5), col=rgb(0.1, 0.1, 0.1, 0.2))
    segments(x0=1:nlags, y0=rep(0, nlags), x1=1:nlags, y1=vals_to_plot)
    points(1:nlags, vals_to_plot, pch=20, col="blue")
    abline(h=c(-1.96*length(data)^{-1/2}, 1.96*length(data)^{-1/2}), col=rgb(0, 0, 1, 0.8), lty=2)
  }

  plot_qr_acf(qr_acf)
  plot_qr_acf(qrsquare_acf)

  if(plot_indstats == TRUE) {
    # Obtain tests statistics
    qrtest <- quantileResidualTests(gsmar, lagsAC=1:nlags, lagsCH=1:nlags, nsimu=nsimu, printRes=FALSE)

    plot_inds <- function(which_ones) { # ac_res or ch_res
      res <- qrtest[[which(names(qrtest) == which_ones)]]
      inds_normalized <- res$indStat/res$stdError
      plot(0, 0, type="n", yaxt="n", xlim=c(0, nlags + 0.1), ylim=c(-3, 3), xlab="", ylab="",
           main=ifelse(which_ones == "ac_res", "Autocovariance statistics", "Cond. h.sked. statistics"))
      abline(h=0, lty=1)
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

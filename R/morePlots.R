#' @title Conditional mean or variance plot for GMAR, StMAR, and G-StMAR models
#'
#' @description \code{condmomentPlot} plots the in-sample conditional means/variances of the model along with
#'  the time series contained in the model (e.g. the time series the model was fitted to). Also plots
#'  the regimiwise conditional means/variances multiplied with the mixing weights.
#'
#' @inheritParams simulateGSMAR
#' @param which_moment should conditional means or variances be plotted?
#' @details The conditional mean plot works best if the data contains positive values only.
#' @return \code{condmomentPlot} only plots to a graphical device and does not return anything. Numerical values
#'  of the conditional means/variances can be extracted from the model with the dollar sign.
#' @inherit simulateGSMAR references
#' @seealso \code{\link{profile_logliks}}, \code{\link{diagnosticPlot}}, \code{\link{fitGSMAR}}, \code{\link{GSMAR}}, \code{\link{quantileResidualTests}},
#'  \code{\link{quantileResidualPlot}}
#' @examples
#' \donttest{
#' # GMAR model
#' fit12 <- fitGSMAR(simudata, p=1, M=2, model="GMAR")
#' condmomentPlot(fit12, which_moment="mean")
#' condmomentPlot(fit12, which_moment="variance")
#'
#' # Restricted StMAR model: plot also the individual statistics with
#' # their approximate critical bounds using the given data
#' fit42r <- fitGSMAR(T10Y1Y, p=4, M=2, model="StMAR", restricted=TRUE)
#' condmomentPlot(fit42r, which_moment="mean")
#' condmomentPlot(fit42r, which_moment="variance")
#'
#' # G-StMAR model with one GMAR type and one StMAR type regime
#' fit42g <- fitGSMAR(T10Y1Y, p=4, M=c(1, 1), model="G-StMAR")
#' condmomentPlot(fit42g, which_moment="mean")
#' condmomentPlot(fit42g, which_moment="variance")
#' }
#' @export

condmomentPlot <- function(gsmar, which_moment=c("mean", "variance")) {
  check_gsmar(gsmar)
  check_data(gsmar)
  which_moment <- match.arg(which_moment)
  p <- gsmar$model$p
  M <- sum(gsmar$model$M)
  data <- gsmar$data

  if(which_moment == "mean") {
    total_moments <- gsmar$total_cmeans
    mw_x_reg <- gsmar$mixing_weights*gsmar$regime_cmeans
    vals <- c(as.vector(mw_x_reg), data)
  } else {
    total_moments <- gsmar$total_cvars
    mw_x_reg <- gsmar$mixing_weights*gsmar$regime_cvars
    vals <- mw_x_reg
  }

  old_par <- par(no.readonly=TRUE)
  on.exit(par(old_par))
  par(mar=c(2.6, 2.6, 2.6, 2.6))
  ymin <- floor(min(vals))
  ymax <- ceiling(max(vals))
  make_ts <- function(dat) ts(c(rep(NA, p), dat), start=start(data), frequency=frequency(data))

  if(which_moment == "mean") {
    plot(data, ylim=c(ymin, ymax), xlab="", ylab="", main="Conditional means")
    lines(make_ts(total_moments), col="grey", lty=2, lwd=2)
  } else {
    plot(data, xlab="", ylab="", main="Conditional variances")
    par(new=TRUE)
    plot(make_ts(total_moments), ylim=c(ymin, ymax), col="grey", lty=2, lwd=2, xlab="", ylab="", yaxt="n", xaxt="n")
    axis(side=4, col="grey", lwd=2)
  }
  colpal <- grDevices::colorRampPalette(c("blue", "turquoise1", "green", "red"))(M)
  for(i1 in 1:M) {
    lines(make_ts(mw_x_reg[,i1]), col=colpal[i1], lty=3)
  }
  legend("topleft", legend=c("total", paste0("regime ", 1:M)), bty="n", col=c("grey", colpal),
         lty=c(2, rep(3, M)), lwd=2, text.font=2, cex=0.65, x.intersp=0.5, y.intersp=1)
}

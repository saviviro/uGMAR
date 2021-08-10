#' @title Conditional mean or variance plot for GMAR, StMAR, and G-StMAR models
#'
#' @description \code{cond_moment_plot} plots the one-step in-sample conditional means/variances of the model along with
#'  the time series contained in the model (e.g. the time series the model was fitted to). Also plots
#'  the regimewise conditional means/variances multiplied with the mixing weights.
#'
#' @inheritParams simulate.gsmar
#' @param which_moment should conditional means or variances be plotted?
#' @details The conditional mean plot works best if the data contains positive values only.
#' @return \code{cond_moment_plot} only plots to a graphical device and does not return anything. Numerical values
#'  of the conditional means/variances can be extracted from the model with the dollar sign.
#' @inherit simulate.gsmar references
#' @seealso \code{\link{profile_logliks}}, \code{\link{diagnostic_plot}}, \code{\link{fitGSMAR}}, \code{\link{GSMAR}}, \code{\link{quantile_residual_tests}},
#'  \code{\link{quantile_residual_plot}}
#' @examples
#' # GMAR model
#' params12 <- c(1.70, 0.85, 0.30, 4.12, 0.73, 1.98, 0.63)
#' gmar12 <- GSMAR(data=simudata, p=1, M=2, params=params12, model="GMAR")
#' cond_moment_plot(gmar12, which_moment="mean")
#' cond_moment_plot(gmar12, which_moment="variance")
#'
#' # G-StMAR model
#' params42gs <- c(0.04, 1.34, -0.59, 0.54, -0.36, 0.01, 0.06, 1.28, -0.36,
#'                 0.2, -0.15, 0.04, 0.19, 9.75)
#' gstmar42 <- GSMAR(data=M10Y1Y, p=4, M=c(1, 1), params=params42gs,
#'                   model="G-StMAR")
#' cond_moment_plot(gstmar42, which_moment="mean")
#' cond_moment_plot(gstmar42, which_moment="variance")
#' @export

cond_moment_plot <- function(gsmar, which_moment=c("mean", "variance")) {
  check_gsmar(gsmar)
  check_data(gsmar)
  which_moment <- match.arg(which_moment)
  p <- gsmar$model$p
  M <- sum(gsmar$model$M)
  data <- gsmar$data

  if(which_moment == "mean") { # Obtain the conditional means
    total_moments <- gsmar$total_cmeans
    mw_x_reg <- gsmar$mixing_weights*gsmar$regime_cmeans
    vals <- c(as.vector(mw_x_reg), data)
    ymin <- floor(min(vals))
  } else { # Obtain the conditional varainces
    total_moments <- gsmar$total_cvars
    mw_x_reg <- gsmar$mixing_weights*gsmar$regime_cvars
    vals <- mw_x_reg
    ymin <- 0
  }

  # Graphical settings
  old_par <- par(no.readonly=TRUE)
  on.exit(par(old_par))
  par(mar=c(2.6, 2.6, 2.6, 2.6))
  ymax <- max(vals)
  make_ts <- function(dat) ts(c(rep(NA, p), dat), start=start(data), frequency=frequency(data))

  if(which_moment == "mean") { # Plot the conditional means
    plot(data, ylim=c(ymin, ymax), xlab="", ylab="", main="Conditional means")
    lines(make_ts(total_moments), col="grey", lty=2, lwd=2)
  } else { # Plot the conditional variances
    plot(data, xlab="", ylab="", main="Conditional variances")
    par(new=TRUE)
    plot(make_ts(total_moments), ylim=c(ymin, ymax), col="grey", lty=2, lwd=2, xlab="", ylab="", yaxt="n", xaxt="n")
    axis(side=4, col="grey", lwd=2)
  }
  colpal <- grDevices::colorRampPalette(c("blue", "turquoise1", "green", "red"))(M)
  for(i1 in 1:M) { # Plot the scaled regimewise conditional moments
    lines(make_ts(mw_x_reg[,i1]), col=colpal[i1], lty=3)
  }
  legend("topleft", legend=c("total", paste0("regime ", 1:M)), bty="n", col=c("grey", colpal),
         lty=c(2, rep(3, M)), lwd=2, text.font=2, cex=0.65, x.intersp=0.5, y.intersp=1)
}

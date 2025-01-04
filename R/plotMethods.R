#' @import graphics
#'
#' @title Plot method for class 'gsmarpred' objects
#'
#' @description \code{plot.gsmarpred} is plot method for class 'gsmarpred' objects
#'
#' @inheritParams predict.gsmar
#' @param x object of class \code{'gsmarpred'} created with \code{predict.gsmar}.
#' @param add_grid should grid a be added to the plots?
#' @param ... arguments passed to function \code{grid}.
#' @details This method is intended for plotting forecasts of GSMAR processes.
#' @export

plot.gsmarpred <- function(x, ..., nt, mix_weights=TRUE, add_grid=TRUE) {

  # Pick statistics etc
  gsmarpred <- x
  data <- as.ts(gsmarpred$gsmar$data)
  n_obs <- length(data)
  q <- gsmarpred$q
  M <- sum(gsmarpred$gsmar$model$M)
  mix_weights <- gsmarpred$mix_weights & M > 1
  mixing_weights <- gsmarpred$gsmar$mixing_weights
  n_mix <- nrow(mixing_weights)

  # Graphical parameters
  old_par <- par(no.readonly=TRUE) # Save old settings
  on.exit(par(old_par)) # Restore the settings before quitting
  if(mix_weights) {
    par(mfrow=c(2, 1), mar=c(2.1, 2.6, 2.1, 1.1))
  } else {
    par(mfrow=c(1, 1), mar=c(2.6, 2.6, 2.1, 1.1))
  }

  # How many observations to plot preceding the prediction
  if(missing(nt)) {
    nt <- round(length(data)*0.15)
  } else {
    stopifnot(nt > 0 & nt %% 1 == 0)
    if(nt > length(data)) {
      warning("nt > length(data); using nt = length(data)")
      nt <- length(data)
    }
  }

  # Function that creates ts objects to plot
  make_ts <- function(a, mix=FALSE, m) { # a is vector of values; if mix == TRUE, use component m mixing weight last observation as the first value
    last_obs <- ifelse(mix, mixing_weights[n_mix, m], data[n_obs])
    ts(c(last_obs, a), start=time(data)[n_obs], frequency=frequency(data))
  }

  # Prediction time series and all values
  ts_pred <- make_ts(gsmarpred$pred)
  ts_dat <- ts(data[(n_obs - nt + 1):n_obs], start=time(data)[n_obs - nt + 1], frequency=frequency(data))
  t0 <- time(ts_pred)
  all_val <- c(ts_dat, ts_pred, gsmarpred$pred_ints)

  # Prediction invervals
  if(gsmarpred$pi_type == "two-sided") {
    # Functions that create prediction interval time series
    ts1_lapply <- function(pred_ints, mix=FALSE, m) lapply(1:(length(q)/2), function(i1) make_ts(pred_ints[,i1], mix, m)) # Lower bounds
    ts2_lapply <- function(pred_ints, mix=FALSE, m) lapply((length(q)/2 + 1):length(q), function(i1) make_ts(pred_ints[,i1], mix, m)) # Upper bounds

    ints1 <- gsmarpred$pred_ints # Process pi
    ints1_mix <- gsmarpred$mix_pred_ints # Mixing weight pi

  } else { # pi_type = upper, lower, or none
    # If upper/lower/none pi, at least one of the "bounds" is just constants
    what_to_rep <- ifelse(gsmarpred$pi_type == "upper", round(min(all_val)) - 3, round(max(all_val)) + 3) # Otherwise "lower" or "none"
    what_to_rep_mix <- ifelse(gsmarpred$pi_type == "upper", 0, 1)

    # Functions that create prediction interval time series
    ts1_lapply <- function(pred_ints, mix=FALSE, m) lapply(seq_along(q), function(i1) make_ts(pred_ints[,i1], mix, m)[-1]) # Upper or lower graphical device box bound (redundant argument)
    ts2_lapply <- function(pred_ints, mix=FALSE, m) lapply(seq_along(q), function(i1) make_ts(pred_ints[,i1], mix, m)) # Make upper or lower prediction bound

    # Constant lower or upper "bound" for one-sided prediction intervals
    ints1 <- matrix(rep(what_to_rep, times=length(q)*length(ts_pred)), ncol=length(q))
    ints1_mix <- array(rep(what_to_rep_mix, times=length(ts_pred)*length(q)*M), dim=c(length(ts_pred), length(q), M))
  }

  # Create the ts objects for confidence intervals
  ts1 <- ts1_lapply(ints1)
  ts2 <- ts2_lapply(gsmarpred$pred_ints)

  # Mixing weight prediction intervals
  if(mix_weights) {
    ts1_mw <- vector(mode="list", length=M) # Sublist for each regime
    ts2_mw <- vector(mode="list", length=M)
    for(m in 1:M) {
      ts1_mw[[m]] <- ts1_lapply(matrix(ints1_mix[, , m], ncol=length(q)), mix=TRUE, m=m)
      ts2_mw[[m]] <- ts2_lapply(matrix(gsmarpred$mix_pred_ints[, , m], ncol=length(q)), mix=TRUE, m=m)
    }
  }

  # Plot the point predictions
  ts.plot(ts_dat, ts_pred, gpars=list(col=c("black", "blue"), lty=1:2,
                                      ylim=c(round(min(all_val)) - 1,
                                             round(max(all_val)) + 1),
                                      main=paste("Forecast", gsmarpred$n_ahead, "steps ahead")))
  if(add_grid) grid(...)

  # Plot the prediction intervals
  draw_poly <- function(ts1_or_ts2, pred_ts, col) polygon(x=c(t0, rev(t0)), y=c(ts1_or_ts2, rev(pred_ts)), col=col, border=NA)
  col_pred <- grDevices::rgb(0, 0, 1, 0.2)
  if(gsmarpred$pi_type %in% c("two-sided", "upper", "lower")) {
    for(i1 in 1:length(gsmarpred$pi)) { # Go through the prediction intervals
      draw_poly(ts1[[i1]], ts_pred, col=col_pred)
      draw_poly(ts2[[i1]], ts_pred, col=col_pred)
    }
  }

  # Plot mixing weight predictions
  if(mix_weights) {

    # Plot point predictions
    colpal_mw <- grDevices::colorRampPalette(c("blue", "turquoise1", "green", "red"))(M)
    colpal_mw2 <- grDevices::adjustcolor(colpal_mw, alpha.f=0.5)
    mix_ts <- ts(mixing_weights[(n_mix - nt + 1):n_mix,], start=time(data)[n_obs - nt + 1],
                 frequency=frequency(data))
    mix_pred_ts <- ts(rbind(mix_ts[nrow(mix_ts),], gsmarpred$mix_pred), start=time(data)[n_obs], frequency=frequency(data))
    ts.plot(mix_ts, mix_pred_ts, gpars=list(col=c(colpal_mw2, colpal_mw), ylim=c(0, 1), lty=c(rep(1, M), rep(2, M))),
            main="Mixing weights")
    legend("topleft", legend=paste0("regime ", 1:M), bty="n", col=colpal_mw, lty=1, lwd=2,
           text.font=2, cex=0.70, x.intersp=0.3, y.intersp=1, seg.len=1, inset=c(-0.01, -0.035))
    if(add_grid) grid(...)

    # Plot mixing weight prediction intervals
    colpal_mw3 <- grDevices::adjustcolor(colpal_mw, alpha.f=0.2)
    if(gsmarpred$pi_type %in% c("two-sided", "upper", "lower")) { # Don't plot if pi_type == "none"
      for(m in 1:M) { # Go through regimes
        for(i1 in 1:length(gsmarpred$pi)) { # Go through the prediction intervals
          draw_poly(ts1_mw[[m]][[i1]], mix_pred_ts[, m], col=colpal_mw3[m])
          draw_poly(ts2_mw[[m]][[i1]], mix_pred_ts[, m], col=colpal_mw3[m])
        }
      }
    }
  }

  invisible(gsmarpred)
}


#' @import graphics
#' @describeIn quantile_residual_tests Plot p-values of the autocorrelation and conditional
#'  heteroskedasticity tests.
#' @inheritParams print.qrtest
#' @export

plot.qrtest <- function(x, ...) {
  # Graphical settings
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  qrtest <- x
  par(mfrow=c(1, 2), mar=c(5.1, 3.1, 3.1, 1.1))

  # Function to plot p-values
  plot_pvalues <- function(which_ones) { # ac_res, ch_res
    res <- qrtest[[which(names(qrtest) == which_ones)]]
    pvals <- res$pvalue
    seq_pvals <- seq_along(pvals)
    plot(pvals, ylim=c(0, 1), xlim=c(min(seq_pvals) - 0.2, max(seq_pvals) + 0.2), ylab="", xlab="lags",
         main=ifelse(which_ones == "ac_res", "Autocorrelation", "Cond. h.skedasticity"),
         xaxt="n", yaxt="n", pch=16, col="blue")
    axis(side=1, at=seq_pvals, labels=res$lags)
    levels <- c(0.01, 0.05, 0.10, seq(from=0.20, to=1.00, by=0.20))
    axis(side=2, at=levels, las=1, cex.axis=0.8)
    abline(h=0, lwd=2)
    abline(h=c(0.01, 0.05, 0.10, 1.00), lty=2, col=c("red", "red", "red", "darkgreen"))
    segments(x0=seq_pvals, y0=0, y1=pvals, x1=seq_pvals, ...)
    points(pvals)
  }

  # Plot the p-values
  plot_pvalues("ac_res") # Autocorrelation tests
  plot_pvalues("ch_res") # Conditional heteroskedasticity tests
  invisible(qrtest)
}


#' @import graphics
#' @describeIn GSMAR Plot method for class 'gsmar'
#' @param x object of class \code{'gsmar'} created with \code{fitGSMAR} or \code{GSMAR}.
#' @param ... in the plot method: arguments passed to the function \code{density} which
#'   calculates the kernel density estimate of the data.
#' @param include_dens Plot also kernel density estimate of the data and model implied stationary
#'  density with regimewise densities? See the details.
#' @details If \code{include_dens == TRUE}, the kernel density estimate of the data is calculated
#'   with the function \code{density} from the package \code{stats}. By default, the default
#'   settings of that function are used, including the usage of Gaussian kernel. Use the dot
#'   parameters to adjust the settings if desired.
#'
#'   By the model implied stationary density we mean the stationary one-dimensional mixture
#'   density of M regimes (see KMS 2015, Theorem 1 and the discussion following it for the Gaussian
#'   case and Theorem 2 in PMS 2018 for the Student's t case). The regimewise densities (i.e. each
#'   density 1,...,M in the stationary mixture density) are multiplied with the mixing weight
#'   parameters accordingly.
#'
#'   In the density plot black represents the kernel density estimate of the data, grey dashed line the
#'   model implied density, and the colored dotted lines the regime wise densities.
#' @export

plot.gsmar <- function(x, ..., include_dens=TRUE) {

  # Pick the relevant statistics etc
  gsmar <- x
  check_data(gsmar)
  data <- as.ts(gsmar$data)
  n_obs <- length(data)
  p <- gsmar$model$p
  M <- sum(gsmar$model$M)
  ts_mw <- ts(rbind(matrix(NA, nrow=p, ncol=M), as.matrix(gsmar$mixing_weights)),
              start=start(data), frequency=frequency(data)) # First p observations are starting values

  # Graphical settings
  old_par <- par(no.readonly=TRUE)
  on.exit(par(old_par))
  if(include_dens) {
    layout(mat=matrix(c(1, 3, 2, 3), nrow=2, byrow=TRUE), widths=c(2, 1))
    par(mar=c(2.5, 2.3, 1.5, 0.3))
  } else {
    par(mfrow=c(2, 1), mar=c(2.5, 2.5, 2.1, 1))
  }
  colpal_mw <- grDevices::colorRampPalette(c("blue", "turquoise1", "green", "red"))(M)

  # Plot the time series and mixing weights
  ts.plot(data, gpars=list(main="Time series"))
  ts.plot(ts_mw, gpars=list(main="Mixing weights", ylim=c(0, 1), col=colpal_mw, lty=2))
  legend("topleft", legend=paste0("regime ", 1:M), bty="n", col=colpal_mw, lty=1, lwd=2,
         text.font=2, cex=0.75, x.intersp=0.2, y.intersp=1, seg.len=1, inset=c(-0.03, -0.035))

  # Plot kernel density estimate of the data with the model implied density
  if(include_dens) {

    # Collect the required statistics
    params <- gsmar$params
    M_orig <- gsmar$model$M
    model <- gsmar$model$model
    restricted <- gsmar$model$restricted
    constraints <- gsmar$model$constraints
    alphas <- pick_alphas(p=p, M=M_orig, params=params, model=model, restricted=restricted, constraints=constraints)
    means <- get_regime_means(gsmar)
    vars <- get_regime_vars(gsmar)
    if(model == "GMAR") {
      M1 <- M
      M2 <- 0
    } else if(model == "StMAR") {
      M1 <- 0
      M2 <- M
    } else { # model == "G-StMAR"
      M1 <- M_orig[1]
      M2 <- M_orig[2]
    }

    # Degrees of freedom and Student's t density
    if(model %in% c("StMAR", "G-StMAR")) {
      dfs <- pick_dfs(p=p, M=M_orig, params=params, model=model)

      # The student's density function as in the model definition
      my_dt <- function(y, mean, var, df) {
        tmp <- lgamma(0.5*(1 + df)) - lgamma(0.5*df) - 0.5*log(pi*(df - 2)) - 0.5*log(var) -
          0.5*(1 + df)*log(1 + (y - mean)^2/(var*(df - 2)))
        exp(tmp)
      }
    }

    # Regime density multiplied with the corresponding mixing weight parameter
    reg_dens <- function(m, xx) {
      if(m <= M1) {
        return(alphas[m]*dnorm(xx, mean=means[m], sd=sqrt(vars[m])))
      } else {
        return(alphas[m]*my_dt(xx, mean=means[m], var=vars[m], df=dfs[m - M1]))
      }
    }

    # The model density
    mod_dens_f <- function(xx) {
      rowSums(vapply(1:M, function(m) reg_dens(m, xx), numeric(length(xx))))
    }

    # Densities and figure bounds
    data_dens <- density(data)
    mod_mean <- gsmar$uncond_moments$uncond_mean
    mod_sd <- sqrt(gsmar$uncond_moments$uncond_var)
    x0 <- min(mod_mean - 3*mod_sd, min(data_dens$x))
    x1 <- max(mod_mean + 3*mod_sd, max(data_dens$x))
    xpp <- seq(from=x0, to=x1, length.out=500)
    mod_dens <- mod_dens_f(xpp)
    y0 <- 0
    y1 <- max(c(data_dens$y, mod_dens))

    # Plot the densities
    col_seriesdens <- "black"
    col_modeldens <- "darkgrey"
    plot(x=data_dens$x, y=data_dens$y, xlim=c(x0, x1), ylim=c(y0, y1), main="Density",
         ylab="", xlab="", cex.axis=0.8, font.axis=2, lwd=1, type="l", col=col_seriesdens)
    lines(x=xpp, y=mod_dens, type="l", lty=2, lwd=2, col=col_modeldens)
    for(m in 1:M) {
      lines(x=xpp, y=reg_dens(m, xx=xpp), type="l", lty=3, lwd=2, col=colpal_mw[m])
    }
    legend("topleft", legend=c("series", "model", paste0("regime ", 1:M)),
           bty="n", col=c(col_seriesdens, col_modeldens, colpal_mw), lty=c(1, 2, rep(3, times=M)), lwd=2,
           text.font=2, cex=0.75, x.intersp=0.2, y.intersp=1, seg.len=1, inset=c(-0.075, -0.01))
  }

  invisible(gsmar)
}

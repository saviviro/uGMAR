#' @import graphics
#'
#' @title Plot method for class 'gsmarpred' objects
#'
#' @description \code{plot.gsmarpred} is plot method for class 'gsmarpred' objects
#'
#' @param x object of class \code{'gsmarpred'} created with \code{predict.gsmar}.
#' @param nt a positive integer specifying the number of observations to be plotted
#'   along with the prediction. Default is \code{round(length(data)*0.2)}.
#' @param add_grid should grid a be added to the plots?
#' @param ... arguments passed to function \code{grid}.
#' @details This method is intended for plotting forecasts of GSMAR processes.
#' @examples
#'  \donttest{
#'  # GMAR-model
#'  params12 <- c(0.18281409, 0.92657275, 0.00214552,
#'   0.85725129, 0.68210294, 0.01900299, 0.88342018)
#'  gmar12 <- GSMAR(logVIX, 1, 2, params12)
#'  pred <- predict(gmar12, n_ahead=10, plotRes=FALSE, pi=c(0.8, 0.9, 0.99), pi_type="two-sided")
#'  plot(pred, nt=50)
#'  }
#' @export

plot.gsmarpred <- function(x, ..., nt, add_grid=TRUE) {
  gsmarpred <- x
  data <- as.ts(gsmarpred$gsmar$data)
  n_obs <- length(data)
  q <- gsmarpred$q
  if(missing(nt)) {
    nt <- round(length(data)*0.2)
  } else {
    stopifnot(nt > 0 & nt %% 1 == 0)
    if(nt > length(data)) {
      warning("nt > length(data); using nt = length(data)")
      nt <- length(data)
    }
  }
  make_ts <- function(a) ts(c(data[n_obs], a), start=time(data)[n_obs], frequency=frequency(data))
  ts_pred <- make_ts(gsmarpred$pred)
  ts_dat <- ts(data[(n_obs - nt + 1):n_obs], start=time(data)[n_obs - nt + 1], frequency=frequency(data))
  t0 = time(ts_pred)
  all_val <- c(ts_dat, ts_pred, gsmarpred$pred_ints)

  if(gsmarpred$pi_type == "two-sided") {
    ts1 <- lapply(1:(length(q)/2), function(i1) make_ts(gsmarpred$pred_ints[,i1]))
    ts2 <- lapply((length(q)/2 + 1):length(q), function(i1) make_ts(gsmarpred$pred_ints[,i1]))

  } else {
    what_to_rep <- ifelse(gsmarpred$pi_type == "upper", round(min(all_val)) - 3, round(max(all_val)) + 3) # Otherwise "lower" or "none"
    ts0 <- rep(what_to_rep, times=length(ts_pred))
    ts1 <- lapply(seq_along(q), function(i1) make_ts(ts0)[-1]) # Upper or lower graphical device box bound
    ts2 <- lapply(seq_along(q), function(i1) make_ts(gsmarpred$pred_ints[,i1])) # Upper or lower prediction bound
  }

  ts.plot(ts_dat, ts_pred, gpars=list(col=c("black", "blue"), lty=1:2,
                                      ylim=c(round(min(all_val)) - 1,
                                             round(max(all_val)) + 1),
                                      main=paste("Forecast", gsmarpred$n_ahead, "steps ahead")))
  if(add_grid) grid(...)

  if(gsmarpred$pi_type %in% c("two-sided", "upper", "lower")) {
    for(i1 in 1:length(gsmarpred$pi)) { # Go through the prediction intervals
      polygon(x=c(t0, rev(t0)), y=c(ts1[[i1]], rev(ts_pred)), col=grDevices::rgb(0, 0, 1, 0.2), border=NA)
      polygon(x=c(t0, rev(t0)), y=c(ts2[[i1]], rev(ts_pred)), col=grDevices::rgb(0, 0, 1, 0.2), border=NA)
    }
  }
  invisible(gsmarpred)
}


#' @import graphics
#' @describeIn quantileResidualTests Plot p-values of the autocorrelation and conditional
#'  heteroskedasticity tests.
#' @inheritParams print.qrtest
#' @export

plot.qrtest <- function(x, ...) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  qrtest <- x
  par(mfrow=c(1, 2), mar=c(5.1, 3.1, 3.1, 1.1))

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

  plot_pvalues("ac_res")
  plot_pvalues("ch_res")
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
#'   case and Theorem 2 in PMS 2018 for the Studen's t case). The regimewise densities (i.e. each
#'   density 1,...,M in the stationary mixture density) are multiplied with the mixing weight
#'   parameters accordingly.
#'
#'   In the density plot black represents the kernel density estimate of the data, grey dashed line the
#'   model implied density, and the colored dotter lines the regime wise densities.
#' @export

plot.gsmar <- function(x, ..., include_dens=TRUE) {
  gsmar <- x
  check_data(gsmar)
  data <- as.ts(gsmar$data)
  n_obs <- length(data)
  p <- gsmar$model$p
  M <- sum(gsmar$model$M)
  ts_mw <- ts(rbind(matrix(NA, nrow=p, ncol=M), as.matrix(gsmar$mixing_weights)),
              start=start(data), frequency=frequency(data)) # First p observations are starting values

  old_par <- par(no.readonly=TRUE)
  on.exit(par(old_par))
  if(include_dens) {
    layout(mat=matrix(c(1, 3, 2, 3), nrow=2, byrow=TRUE), widths=c(2, 1))
    par(mar=c(2.5, 2.3, 1.5, 0.3))
  } else {
    par(mfrow=c(2, 1), mar=c(2.5, 2.5, 2.1, 1))
  }
  colpal_mw <- grDevices::colorRampPalette(c("blue", "turquoise1", "green", "red"))(M)

  ts.plot(data, gpars=list(main="Time series"))
  ts.plot(ts_mw, gpars=list(main="Mixing weights", ylim=c(0, 1), col=colpal_mw, lty=2))
  legend("topleft", legend=paste0("regime ", 1:M), bty="n", col=colpal_mw, lty=1, lwd=2,
         text.font=2, cex=0.65, x.intersp=0.5, y.intersp=1)

  # Plot kernel density estimate of the data with the model implied density
  if(include_dens) {

    # Collected the required statistics
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
    } else { # Model == "G-StMAR"
      M1 <- M_orig[1]
      M2 <- M_orig[2]
    }

    # Degrees of freedom and Student's t density
    if(model %in% c("StMAR", "G-StMAR")) {
      dfs <- pick_dfs(p=p, M=M_orig, params=params, model=model)

      # The student's density as in the model definition
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
    plot(x=data_dens$x, y=data_dens$y, xlim=c(x0, x1), ylim=c(y0, y1), main="Density",
         ylab="", xlab="", cex.axis=0.8, font.axis=2, type="l")
    lines(x=xpp, y=mod_dens, type="l", lty=2, lwd=2, col="darkgrey")
    for(m in 1:M) {
      lines(x=xpp, y=reg_dens(m, xx=xpp), type="l", lty=3, col=colpal_mw[m])
    }
  }

  invisible(gsmar)
}

#' @import stats
#' @import graphics
#' @importFrom grDevices rgb
#'
#' @title Forecast GMAR, StMAR or G-StMAR process
#'
#' @description \code{forecastGMAR} forecasts the specified GMAR, StMAR or G-StMAR process by using the given data to simulate
#'  its possible future values. For one-step forecasts use of the exact formula for the conditional mean is supported.
#'
#' @inheritParams loglikelihood
#' @param data a numeric vector or column matrix or an (univariate) time series object containing the data. \code{NA} values are not supported.
#' @param nsteps a positive integer specifying how many steps in the future should be forecasted.
#' @param conflevel an (optional) numeric vector whose elements are in the open interval \eqn{(0,1)},
#'  specifying the confidence levels for the confidence intervals that should be calculated. Default is \code{c(0.95, 0.8)}.
#' @param nsimu an (optional) positive integer specifying to how many simulations the forecast should be based on. Default is 10000.
#' @param printRes an (optional) logical argument defining wether results should be printed or not. Default is \code{TRUE}.
#' @param plotRes an (optional) logical argument defining wether the forecast should be plotted or not. Default is \code{TRUE}.
#' @param nt an (optional) positive integer specifying the number of observations to be plotted along with the prediction. Default is \code{round(length(data)*0.2)}.
#' @param useMean set \code{TRUE} if the prediction should be based on sample mean instead of median. Default is FALSE.
#' @param oneStepCondMean set \code{TRUE} for optimal one-step forecast using the exact formula for the conditional mean. Default is \code{FALSE}.
#' @param upperOnly set \code{TRUE} for upper confidence interval only. Default is \code{FALSE}.
#' @param lowerOnly set \code{TRUE} for lower confidence interval only. Default is \code{FALSE}.
#' @details \code{forecastGMAR} uses the last \code{p} values of the given data to simulate \code{nsimu} possible future values for each step.
#'   The best prediction is then obtained by calculating the sample median (or mean) of each step and the confidence intervals are obtained from the empirical fractiles.
#'
#' @return Returns a data frame containing the empirical best predicton and confidence intervals accordingly to \code{conflevel}.
#'   Or if \code{oneStepCondMean==TRUE} returns the optimal prediction as (1x1) numeric vector.
#' @inherit simulateGMAR references
#' @examples
#' \donttest{
#' # GMAR model
#' params12 <- c(1.1, 0.9, 0.3, 4.5, 0.7, 3.2, 0.8)
#' pred12 <- forecastGMAR(VIX, 1, 2, params12, nsteps=10)
#'
#' # Restricted GMAR model
#' params12r <- c(1.4, 1.8, 0.9, 0.3, 3.2, 0.8)
#' pred12r <- forecastGMAR(VIX, 1, 2, params12r, restricted=TRUE, nsteps=20,
#'                         conflevel=c(0.9, 0.8, 0.6), nt=200, upperOnly=TRUE)
#'
#' # StMAR model (one-step predictions by simulation and by exact formula)
#' params12t <- c(1.1, 0.9, 0.3, 4.5, 0.7, 3.2, 0.8, 5, 8)
#' pred12t <- forecastGMAR(VIX, 1, 2, params12t, StMAR=TRUE, nsteps=1, useMean=TRUE)
#' pred12te <- forecastGMAR(VIX, 1, 2, params12t, StMAR=TRUE, oneStepCondMean=TRUE)
#'
#' # Non-mixture version of StMAR model with data as (fictional) time series object
#' params11t <- c(0.76, 0.93, 1.4, 2.4)
#' pred11t <- forecastGMAR(ts(VIX, start=1900, freq=12), 1, 1, params11t,
#'                         StMAR=TRUE, nsteps=5, conflevel=c(0.99))
#'
#' # G-StMAR model
#' params13gs <- c(2.0, 0.83, 0.36, 1.14, 0.90, 0.06, 4.23, 0.72, 3.85, 0.6, 0.20, 3.3)
#' pred13gs <- forecastGMAR(VIX, 1, c(2, 1), params13gs, GStMAR=TRUE, nsteps=10)
#'
#' # G-StMAR model: optimal one-step forecast using the exact formula of conditional mean
#' params12gs <- c(3.98, 0.68, 0.36, 0.70, 0.94, 11.75, 0.25, 2.03)
#' pred12gs <- forecastGMAR(VIX, 1, c(1, 1), params12gs, GStMAR=TRUE, oneStepCondMean=TRUE)
#'
#' # GMAR model as a mixture of AR(2) and AR(1) models
#' R <- list(diag(1, ncol=2, nrow=2), as.matrix(c(1, 0)))
#' params22c <- c(1.2, 0.8, 0.1, 0.3, 3.3, 0.8, 2.8, 0.8)
#' pred22c <- forecastGMAR(VIX, 2, 2, params22c, constraints=TRUE, R=R,
#'                         nsteps=15, conflevel=c(0.99, 0.9, 0.8))
#'
#' # Such StMAR(3,2) that the AR coefficients are restricted to be
#' # the same for both regimes and that the second AR coefficients are
#' # constrained to zero.
#' params32trc <- c(2.2, 1.8, 0.88, -0.03, 2.4, 0.27, 0.40, 3.9, 1000)
#' pred32trc <- forecastGMAR(VIX, 3, 2, params32trc, StMAR=TRUE, restricted=TRUE, constraints=TRUE,
#'                           R=matrix(c(1, 0, 0, 0, 0, 1), ncol=2), nsteps=5)
#' }
#' @export

forecastGMAR <- function(data, p, M, params, StMAR=FALSE, GStMAR=FALSE, restricted=FALSE, constraints=FALSE, R, nsteps, conflevel=c(0.95, 0.8), nsimu=10000, printRes=TRUE, plotRes=TRUE, nt, useMean=FALSE, oneStepCondMean=FALSE, upperOnly=FALSE, lowerOnly=FALSE) {

  checkLogicals(StMAR=StMAR, GStMAR=GStMAR)
  checkPM(p, M, GStMAR=GStMAR)
  n_obs = length(data)
  if(upperOnly==TRUE & lowerOnly==TRUE) {
    stop("Arguments upperOnly and lowerOnly are both set TRUE, which doesn't make any sense. Set both FALSE (default) for two-sided confidence intervals.")
  }
  if(missing(nt)) {
    nt = round(n_obs*0.2)
  }
  if(!is.ts(data)) {
    data = as.ts(data)
  }
  if(oneStepCondMean==TRUE) {
    if(missing(nsteps)) {
      nsteps = 1
    } else if(nsteps!=1) {
      print("Exact conditional expectation is supported for one-step forecasts only!")
      nsteps = 1
    }
  }
  if(any(conflevel>=1) | any(conflevel<=0)) {
    stop("Confidence level has to be in the open interval (0,1)")
  } else if(nsteps<1) {
    stop("Number of steps has to be equal or larger than one")
  } else if(nsimu<1) {
    stop("Number of simulations has to be equal or largen than one")
  } else if(length(conflevel)<1) {
    stop("The lengths of vector 'conflevel' has to be equal or greater than one")
  } else if(n_obs<p) {
    stop("The data should contain p values at minimum")
  }
  # Simulate future values of the process
  initvalues = data[(n_obs-p+1):n_obs]
  if(any(is.na(initvalues))) {
    stop("The last p values of the given data contain NA values, which is not supported")
  }
  if(oneStepCondMean==TRUE) { # Exact optimal forecast
    mw = mixingWeights(data, p, M, params, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted, constraints=constraints, R=R)
    pars = reformParameters(p, M, params, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted)$pars
    if(p==1) {
      bestcast = sum(mw[nrow(mw),]*(pars[1,] + initvalues*pars[2:(2+p-1),]))
    } else {
      bestcast = sum(mw[nrow(mw),]*(pars[1,] + t(rev(initvalues))%*%pars[2:(2+p-1),]))
    }
    if(printRes==TRUE) {
      print(paste0("Optimal one-step prediction: ", round(bestcast, 3)))
    }
    if(plotRes==TRUE) {
      pred = ts(c(data[n_obs], bestcast), start=time(data)[n_obs], frequency=frequency(data))
      dat = ts(data[(n_obs-nt):n_obs], start=time(data)[n_obs-nt], frequency=frequency(data))
      ts.plot(dat, pred, gpars=list(col=c("black", "blue"), lty=c(1, 2)), ylim=c(round(min(min(dat), bestcast)-max(0.1*bestcast, 1)),
                                                                                 round(max(max(dat), bestcast)+max(0.1*bestcast, 1))))
    }
    return(bestcast)
  } else { # If forecast by simulation is considered
    res = simulateGMAR(p, M, params, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted, constraints=constraints,
                       R=R, nsimu=nsteps, initvalues=initvalues, ntimes=nsimu)

    # Calculate forecast and confidence intervals
    if(useMean==TRUE) {
      bestcast = rowMeans(res)
    } else {
      bestcast = vapply(1:nsteps, function(i1) median(res[i1,]), numeric(1))
    }
    if(upperOnly==TRUE) {
      confq = conflevel
    } else if(lowerOnly==TRUE) {
      confq = 1-conflevel
    } else {
      confq = as.vector(vapply(1:length(conflevel), function(i1) c((1-conflevel[i1])/2, 1-(1-conflevel[i1])/2), numeric(2)))
    }
    confq = confq[order(confq, decreasing=FALSE)]
    if(is.vector(res)) {
      confints = as.matrix(quantile(res, confq))
    } else {
      confints = vapply(1:nsteps, function(i1) quantile(res[i1,], confq), numeric(length(confq)))
    }
    if(upperOnly==TRUE) {
      uppers = confints
      lowers = NULL
    } else if(lowerOnly==TRUE) {
      uppers = NULL
      lowers = confints
    } else {
      lowers = confints[1:(length(confq)/2),]
      uppers = confints[(length(confq)/2+1):length(confq),]
    }
    if(length(conflevel)==1) {
      if(!is.null(lowers)) {
        lowers = t(lowers)
      }
      if(!is.null(uppers)) {
        uppers = t(uppers)
      }
    }
    if(upperOnly==TRUE) {
      fcast = data.frame(bestcast, t(uppers))
    } else if(lowerOnly==TRUE) {
      fcast = data.frame(t(lowers), bestcast)
    } else {
      fcast = data.frame(t(lowers), bestcast, t(uppers))
    }
    if(useMean==TRUE) {
      if(upperOnly==TRUE) {
        colnames(fcast) = c("MEAN", confq[1:length(confq)])
      } else if(lowerOnly==TRUE) {
        colnames(fcast) = c(confq[1:length(confq)], "MEAN")
      } else {
        colnames(fcast) = c(confq[1:(length(confq)/2)], "MEAN", confq[(length(confq)/2+1):length(confq)])
      }
    } else { # Uses median
      if(upperOnly==TRUE) {
        colnames(fcast) = c(0.5, confq[1:length(confq)])
      } else if(lowerOnly==TRUE) {
        colnames(fcast) = c(confq[1:length(confq)], 0.5)
      } else {
        colnames(fcast) = c(confq[1:(length(confq)/2)], 0.5, confq[(length(confq)/2+1):length(confq)])
      }
    }
    if(printRes==TRUE) {
      print(round(fcast, 2))
    }
    if(plotRes==TRUE) {
      if(nsteps==1) {
        if(!is.null(lowers)) {
          lowers=matrix(t(lowers))
        }
        if(!is.null(uppers)) {
          uppers=matrix(t(uppers))
        }
      }
      if(upperOnly==TRUE) {
        ts_uppers = lapply(1:length(confq), function(i1) ts(c(data[n_obs], uppers[i1,]), start=time(data)[n_obs], frequency=frequency(data)))
        pred = ts(c(data[n_obs], bestcast), start=time(data)[n_obs], frequency=frequency(data))
        dat = ts(data[(n_obs-nt):n_obs], start=time(data)[n_obs-nt], frequency=frequency(data))
        t0 = time(pred); values0 = c(as.vector(confints), as.vector(dat))
        ts.plot(dat, pred, gpars=list(col=c("black", "blue"), lty=c(1, 2)), ylim=c(round(min(values0))-1, round(max(values0))+1))
        ts_lowers = ts(rep(round(min(values0))-2, times=nsteps+1), start=time(data)[n_obs], frequency=frequency(data))
        for(i1 in 1:length(conflevel)) {
          polygon(x=c(t0, rev(t0)), y=c(ts_lowers, rev(pred)), col=rgb(0, 0, 1, 0.2), border=NA)
          polygon(x=c(t0, rev(t0)), y=c(ts_uppers[[i1]], rev(pred)), col=rgb(0, 0, 1, 0.2), border=NA)
        }
      } else if(lowerOnly==TRUE) {
        ts_lowers = lapply(1:length(confq), function(i1) ts(c(data[n_obs], lowers[i1,]), start=time(data)[n_obs], frequency=frequency(data)))
        pred = ts(c(data[n_obs], bestcast), start=time(data)[n_obs], frequency=frequency(data))
        dat = ts(data[(n_obs-nt):n_obs], start=time(data)[n_obs-nt], frequency=frequency(data))
        t0 = time(pred); values0 = c(as.vector(confints), as.vector(dat))
        ts.plot(dat, pred, gpars=list(col=c("black", "blue"), lty=c(1, 2)), ylim=c(round(min(values0))-1, round(max(values0))+1))
        ts_uppers = ts(rep(round(max(values0))+2, times=nsteps+1), start=time(data)[n_obs], frequency=frequency(data))
        for(i1 in 1:length(conflevel)) {
          polygon(x=c(t0, rev(t0)), y=c(ts_lowers[[i1]], rev(pred)), col=rgb(0, 0, 1, 0.2), border=NA)
          polygon(x=c(t0, rev(t0)), y=c(ts_uppers, rev(pred)), col=rgb(0, 0, 1, 0.2), border=NA)
        }
      } else { # If two-sided
        ts_lowers = lapply(1:(length(confq)/2), function(i1) ts(c(data[n_obs], lowers[i1,]), start=time(data)[n_obs], frequency=frequency(data)))
        ts_uppers = lapply(1:(length(confq)/2), function(i1) ts(c(data[n_obs], uppers[i1,]), start=time(data)[n_obs], frequency=frequency(data)))
        pred = ts(c(data[n_obs], bestcast), start=time(data)[n_obs], frequency=frequency(data))
        dat = ts(data[(n_obs-nt):n_obs], start=time(data)[n_obs-nt], frequency=frequency(data))
        t0 = time(pred)
        values0 = c(as.vector(confints), as.vector(dat))
        ts.plot(dat, pred, gpars=list(col=c("black", "blue"), lty=c(1, 2)), ylim=c(round(min(values0))-1, round(max(values0))+1))
        for(i1 in 1:length(conflevel)) {
          polygon(x=c(t0, rev(t0)), y=c(ts_lowers[[i1]], rev(pred)), col=rgb(0, 0, 1, 0.2), border=NA)
          polygon(x=c(t0, rev(t0)), y=c(ts_uppers[[i1]], rev(pred)), col=rgb(0, 0, 1, 0.2), border=NA)
        }
      }
    }
    return(fcast)
  }
}


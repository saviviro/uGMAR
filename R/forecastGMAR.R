#' @import stats
#' @import graphics
#' @importFrom grDevices rgb
#'
#' @title Forecast GMAR pr StMAR process
#'
#' @description \code{forecastGMAR} forecasts the specified GMAR or StMAR process by using the given data to simulate
#'  its possible future values.
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
#' @param useMean set TRUE if the prediction should be based on sample mean instead of median. Default is FALSE.
#' @details \code{forecastGMAR} uses the last \code{p} values of the given data to simulate \code{nsimu} possible future values for each step.
#'   The best prediction is then obtained by calculating the sample median (or mean) of each step and the confidence intervals are obtained from the empirical fractiles.
#'
#' @return Returns a data frame containing the empirical best predicton and confidence intervals accordingly to \code{conflevel}.
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
#'                         conflevel=c(0.9, 0.8, 0.6), nt=200)
#'
#' # StMAR model
#' params12t <- c(1.1, 0.9, 0.3, 4.5, 0.7, 3.2, 0.8, 5, 8)
#' pred12t <- forecastGMAR(VIX, 1, 2, params12t, StMAR=TRUE, nsteps=1)
#'
#' # Non-mixture version of StMAR model with data as (fictional) time series object
#' params11t <- c(0.76, 0.93, 1.4, 2.4)
#' pred11t <- forecastGMAR(ts(VIX, start=1900, freq=12), 1, 1, params11t,
#'                         StMAR=TRUE, nsteps=5, useMean=TRUE)
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

forecastGMAR <- function(data, p, M, params, StMAR=FALSE, restricted=FALSE, constraints=FALSE, R, nsteps, conflevel=c(0.95, 0.8), nsimu=10000, printRes=TRUE, plotRes=TRUE, nt, useMean=FALSE) {

  checkPM(p, M)
  n_obs = length(data)
  if(missing(nt)) {
    nt = round(n_obs*0.2)
  }
  if(!is.ts(data)) {
    data = as.ts(data)
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
    stop("the data should contain p values at minimum")
  }

  # Simulate future values of the process
  initvalues = data[(n_obs-p+1):n_obs]
  if(any(is.na(initvalues))) {
    stop("The last p values of the given data contain NA values, which is not supported")
  }
  res = vapply(1:nsimu, function(x) simulateGMAR(p, M, params, StMAR=StMAR, restricted=restricted, constraints=constraints,
                                                 R=R, nsimu=nsteps, initvalues=initvalues)$sample, numeric(nsteps))
  if(nsteps==1) {
    res = t(res)
  }
  # Calculate forecast and confidence intervals
  if(useMean==TRUE) {
    bestcast = rowMeans(res)
  } else {
    bestcast = vapply(1:nsteps, function(i1) median(res[i1,]), numeric(1))
  }
  confq = as.vector(vapply(1:length(conflevel), function(i1) c((1-conflevel[i1])/2, 1-(1-conflevel[i1])/2), numeric(2)))
  confq = confq[order(confq, decreasing=FALSE)]
  if(is.vector(res)) {
    confints = as.matrix(quantile(res, confq))
  } else {
    confints = vapply(1:nsteps, function(i1) quantile(res[i1,], confq), numeric(length(confq)))
  }
  lowers = confints[1:(length(confq)/2),]
  uppers = confints[(length(confq)/2+1):length(confq),]

  fcast = data.frame(t(lowers), bestcast, t(uppers))
  if(useMean==TRUE) {
    colnames(fcast) = c(confq[1:(length(confq)/2)], "MEAN", confq[(length(confq)/2+1):length(confq)])
  } else {
    colnames(fcast) = c(confq[1:(length(confq)/2)], 0.5, confq[(length(confq)/2+1):length(confq)])
  }
  if(printRes==TRUE) {
    print(round(fcast, 2))
  }
  if(plotRes==TRUE) {
    if(nsteps==1) {
      lowers=matrix(t(lowers))
      uppers=matrix(t(uppers))
    }
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
  return(fcast)
}


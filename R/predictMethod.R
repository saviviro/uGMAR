#' @import stats
#'
#' @title Forecast GMAR, StMAR or G-StMAR process
#'
#' @description \code{predict.gsmar} forecasts the specified GMAR, StMAR or G-StMAR process by using the given data to simulate
#'  its possible future values. For one-step forecasts using the exact formula of conditional mean is supported.
#'
#' @param object object of class \code{'gsmar'} created with function \code{fitGSMAR} or \code{GSMAR}.
#' @param ... additional arguments passed to \code{grid} (ignored if \code{plot_res==FALSE}).
#' @param n_ahead a positive integer specifying how many steps in the future should be forecasted.
#' @param nsimu a positive integer specifying to how many simulations the forecast should be based on.
#' @param ci a numeric vector specifying the confidence levels of the confidence intervals.
#' @param pred_type should the prediction be based on sample "mean" or "median"? Or should it
#'   be one-step-ahead forecast based on conditional mean (\code{"cond_mean"})? Confidence intervals
#'   won't be calculated if conditional mean is used.
#' @param ci_type should the confidence intervals be "two-sided", "upper" or "lower"?
#' @param nt a positive integer specifying the number of observations to be plotted
#'   along with the prediction. Default is \code{round(length(data)*0.2)}.
#' @param plotRes a logical argument defining wether the forecast should be plotted or not.
#' @details \code{predict.gsmar} uses the last \code{p} values of the data to simulate \code{nsimu} possible
#'  future values for each step. The best prediction is then obtained by calculating the sample median (or mean)
#'  of each step and the confidence intervals are obtained from the empirical fractiles.
#'
#'  We encourage directly using the function \code{simulateGSMAR} for quantile based forecasting. With \code{simulateGSMAR}
#'  it's easy to forecast the mixing weights too.
#'
#' @return Returns a data frame containing the empirical best predicton and confidence intervals accordingly to \code{ci}.
#'   Or if \code{pred_type=="cond_mean"} returns the optimal prediction as (1x1) numeric vector.
#' @inherit simulateGSMAR references
#' @seealso \code{\link{simulateGSMAR}}, \code{\link{condMoments}}, \code{\link{fitGSMAR}}, \code{\link{GSMAR}},
#'  \code{\link{quantileResidualTests}}, \code{\link{diagnosticPlot}}
#' @examples
#' \donttest{
#' # GMAR model
#' params12 <- c(1.12, 0.91, 0.29, 4.53, 0.70, 3.21, 0.84)
#' gmar12 <- GSMAR(VIX, 1, 2, params12)
#' pred12 <- predict(gmar12, n_ahead=10)
#' pred12
#'
#' # Restricted GMAR model, one-step conditional mean prediction
#' params12r <- c(1.4, 1.8, 0.88, 0.29, 3.18, 0.84)
#' gmar12r <- GSMAR(data=VIX, p=1, M=2, params=params12r, model="GMAR",
#'  restricted=TRUE)
#' pred12r <- predict(gmar12r, pred_type="cond_mean", plotRes=FALSE)
#' pred12r
#'
#' # StMAR model, upper confidence intervals
#' params12t <- c(1.38, 0.88, 0.27, 3.8, 0.74, 3.15, 0.8, 100, 3.6)
#' stmar12 <- GSMAR(data=VIX, p=1, M=2, params=params12t, model="StMAR")
#' predict(stmar12, n_ahead=10, ci_type="upper", ci=c(0.99, 0.95, 0.9))
#'
#' # G-StMAR model, no confidence intervals
#' params12gs <- c(1.38, 0.88, 0.27, 3.8, 0.74, 3.15, 0.8, 3.6)
#' gstmar12 <- GSMAR(data=VIX, p=1, M=c(1, 1), params=params12gs,
#'  model="G-StMAR")
#' pred12gs <- predict(gstmar12, n_ahead=10, pred_type="median",
#'  ci_type="none", plotRes=FALSE)
#' pred12gs
#' plot(pred12gs)
#'
#' # GMAR model as a mixture of AR(2) and AR(1) models
#' constraints <- list(diag(1, ncol=2, nrow=2), as.matrix(c(1, 0)))
#' params22c <- c(1.2, 0.85, 0.04, 0.3, 3.3, 0.77, 2.8, 0.77)
#' gmar22c <- GSMAR(data=VIX, p=2, M=2, params=params22c,
#'  model="GMAR", constraints=constraints)
#' predict(gmar22c, n_ahead=5, nsimu=10000, nt=10)
#'
#' # Such StMAR(3,2) that the AR coefficients are restricted to be
#' # the same for both regimes and that the second AR coefficients are
#' # constrained to zero.
#' params32trc <- c(2.2, 1.8, 0.88, -0.03, 2.4, 0.27, 0.40, 3.9, 1000)
#' stmar32rc <- GSMAR(data=VIX, p=3, M=2, params=params32trc, model="StMAR",
#'  restricted=TRUE, constraints=matrix(c(1, 0, 0, 0, 0, 1), ncol=2))
#' predict(stmar32rc, n_ahead=3, ci_type="lower")
#' }
#' @export

predict.gsmar <- function(object, ..., n_ahead, nsimu=5000, ci=c(0.95, 0.8), pred_type=c("mean", "median", "cond_mean"),
                         ci_type=c("two-sided", "upper", "lower", "none"), nt, plotRes=TRUE) {
  gsmar <- object
  pred_type <- match.arg(pred_type)
  ci_type <- match.arg(ci_type)
  stopifnot(pred_type %in% c("mean", "median", "cond_mean"))
  stopifnot(ci_type %in% c("two-sided", "upper", "lower", "none"))
  check_gsmar(gsmar)
  check_data(gsmar)
  data <- gsmar$data
  n_obs <- length(data)

  p <- gsmar$model$p
  M <- gsmar$model$M
  params <- gsmar$params
  model <- gsmar$model$model
  restricted <- gsmar$model$restricted
  constraints <- gsmar$model$constraints

  if(pred_type == "cond_mean") {
    if(missing(n_ahead)) {
      n_ahead <- 1
    } else if(n_ahead != 1) {
      warning("Exact conditional expectation is supported for one-step forecasts only! Using n_ahead=1.")
      n_ahead <- 1
    }
  }
  if(missing(nt)) {
    nt <- round(length(data)*0.2)
  } else {
    stopifnot(nt > 0 & nt %% 1 == 0)
    if(nt > length(data)) {
      warning("nt > length(data), using nt = length(data)")
      nt <- length(data)
    }
  }
  if(!all_pos_ints(c(n_ahead, nsimu))) stop("Arguments n_ahead and nsimu must be positive integers")
  if(any(ci >= 1) | any(ci <= 0)) stop("Each confidence level has to be in the open interval (0, 1)")
  if(!is.null(constraints)) checkConstraintMat(p=p, M=M, restricted=restricted, constraints=constraints)

  # Simulate future values of the process
  if(pred_type == "cond_mean") { # Exact optimal forecast

    # Collect parameter values and calculate mixing weights
    if(gsmar$model$parametrization == "mean") {
      params <- change_parametrization(p=p, M=M, params=params, model=model, restricted=restricted,
                                       constraints=constraints, change_to="intercept")
    }
    mw <- mixingWeights_int(data, p, M, params, model=model, restricted=restricted, constraints=constraints,
                            parametrization="intercept", checks=TRUE, to_return="mw_tplus1")
    pars <- pick_pars(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints)

    # Calculate the conditional mean
    pred <- sum(mw[nrow(mw),]*(pars[1,] + t(rev(data[(n_obs - p + 1):n_obs]))%*%pars[2:(2 + p - 1),]))
    conf_ints <- NULL
    ci <- NULL
    ci_type <- "none"
    q_tocalc <- numeric(0)
  } else {

    # Simulations
    res <- simulateGSMAR(gsmar, nsimu=n_ahead, initvalues=data, ntimes=nsimu)$sample

    # Predictions
    if(pred_type == "mean") {
      pred <- rowMeans(res)
    } else {
      pred <- vapply(1:n_ahead, function(i1) median(res[i1,]), numeric(1))
    }

    # Confidence intercals
    if(ci_type == "upper") {
      q_tocalc <- ci
    } else if(ci_type == "lower") {
      q_tocalc <- 1 - ci
    } else if(ci_type == "two-sided") {
      lower <- (1 - ci)/2
      upper <- rev(1 - lower)
      q_tocalc <- c(lower, upper)
    } else {  # If ci_type == "none"
      q_tocalc <- numeric(0)
      ci <- NULL
    }
    q_tocalc <- sort(q_tocalc, decreasing=FALSE)

    if(is.vector(res)) {
      conf_ints <- as.matrix(quantile(res, probs=q_tocalc))
    } else {
      conf_ints <- t(vapply(1:n_ahead, function(i1) quantile(res[i1,], probs=q_tocalc), numeric(length(q_tocalc))))
    }
  }

  ret <- structure(list(gsmar=gsmar,
                        pred=pred,
                        conf_ints=conf_ints,
                        n_ahead=n_ahead,
                        nsimu=nsimu,
                        ci=ci,
                        ci_type=ci_type,
                        pred_type=pred_type,
                        q=q_tocalc),
                   class="gsmarpred")
  if(plotRes == TRUE) plot.gsmarpred(x=ret, nt=nt, ...)
  ret
}


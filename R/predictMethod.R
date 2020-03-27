#' @import stats
#'
#' @title Forecast GMAR, StMAR, or G-StMAR process
#'
#' @description \code{predict.gsmar} forecasts the specified GMAR, StMAR, or G-StMAR process by using the given
#'  data to simulate its possible future values. For one-step forecasts using the exact formula for conditional
#'  mean is supported.
#'
#' @param object object of class \code{'gsmar'} created with function \code{fitGSMAR} or \code{GSMAR}.
#' @param ... additional arguments passed to \code{grid} (ignored if \code{plot_res==FALSE}).
#' @param n_ahead a positive integer specifying how many steps in the future should be forecasted.
#' @param nsimu a positive integer specifying to how many simulations the forecast should be based on.
#' @param pi a numeric vector specifying confidence levels for the prediction intervals.
#' @param pred_type should the prediction be based on sample "median" or "mean"? Or should it
#'   be one-step-ahead forecast based on the exact conditional mean (\code{"cond_mean"})? prediction
#'   intervals won't be calculated if the exact conditional mean is used.
#' @param pi_type should the prediction intervals be "two-sided", "upper", or "lower"?
#' @param plotRes a logical argument defining whether the forecast should be plotted or not.
#' @param mix_weights \code{TRUE} if forecasts for mixing weights should be plotted, \code{FALSE} in not.
#' @param nt a positive integer specifying the number of observations to be plotted
#'   along with the prediction. Default is \code{round(length(data)*0.15)}.
#' @details \code{predict.gsmar} uses the last \code{p} values of the data to simulate \code{nsimu}
#'  possible future values for each step-ahead. The point prediction is then obtained by calculating
#'  the sample median or mean for each step and the prediction intervals are obtained from the
#'  empirical fractiles.
#'
#'  The function \code{simulateGSMAR} can also be used directly for quantile based forecasting.
#'
#' @return Returns a class \code{'gsmarpred'} object containing, among the specifications,...
#'  \describe{
#'    \item{$pred}{Point forecasts}
#'    \item{$pred_ints}{Prediction intervals}
#'    \item{$mix_pred}{Point forecasts for mixing weights}
#'    \item{mix_pred_ints}{Individual prediction intervals for mixing weights, as \code{[, , m]}, m=1,..,M.}
#'  }
#' @inherit simulateGSMAR references
#' @seealso \code{\link{simulateGSMAR}}, \code{\link{condMoments}}, \code{\link{fitGSMAR}}, \code{\link{GSMAR}},
#'  \code{\link{quantileResidualTests}}, \code{\link{diagnosticPlot}}
#' @examples
#' \donttest{
#' # GMAR model
#' fit12 <- fitGSMAR(data=logVIX, p=1, M=2, model="GMAR")
#' pred12 <- predict(fit12, n_ahead=10, pi=c(0.95, 0.8))
#' pred12
#'
#' # Non-mixture StMAR model, upper prediction intervals
#' fit11t <- fitGSMAR(logVIX, 1, 1, model="StMAR", ncores=1, ncalls=1)
#' predict(fit11t, n_ahead=10, pi_type="upper", pi=0.9)
#'
#' # G-StMAR model, no prediction intervals
#' fit12gs <- fitGSMAR(logVIX, 1, M=c(1, 1), model="G-StMAR")
#' pred12gs <- predict(fit12gs, n_ahead=2, pred_type="median",
#'  pi_type="none", plotRes=FALSE)
#' pred12gs
#' plot(pred12gs)
#'
#' # Restricted GMAR model, one-step conditional mean prediction
#' fit12r <- fitGSMAR(logVIX, 1, 2, model="GMAR", restricted=TRUE)
#' pred12r <- predict(fit12r, pred_type="cond_mean", plotRes=FALSE)
#' pred12r
#'
#' # Such StMAR(3,2) that the AR coefficients are restricted to be
#' # the same for both regimes and that the second AR coefficients are
#' # constrained to zero.
#' fit32rc <- fitGSMAR(logVIX, 3, 2, model="StMAR", restricted=TRUE,
#'  constraints=matrix(c(1, 0, 0, 0, 0, 1), ncol=2))
#' predict(fit32rc, n_ahead=3, pi_type="lower")
#' }
#' @export

predict.gsmar <- function(object, ..., n_ahead, nsimu=10000, pi=c(0.95, 0.8), pred_type=c("median", "mean", "cond_mean"),
                         pi_type=c("two-sided", "upper", "lower", "none"), plotRes=TRUE, mix_weights=TRUE, nt) {
  gsmar <- object
  pred_type <- match.arg(pred_type)
  pi_type <- match.arg(pi_type)
  stopifnot(pred_type %in% c("mean", "median", "cond_mean"))
  stopifnot(pi_type %in% c("two-sided", "upper", "lower", "none"))
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
    nt <- round(length(data)*0.15)
  } else {
    stopifnot(nt > 0 & nt %% 1 == 0)
    if(nt > length(data)) {
      warning("nt > length(data); using nt = length(data)")
      nt <- length(data)
    }
  }
  if(!all_pos_ints(c(n_ahead, nsimu))) stop("Arguments n_ahead and nsimu must be positive integers")
  if(any(pi >= 1) | any(pi <= 0)) stop("Each confidence level has to be in the open interval (0, 1)")
  if(!is.null(constraints)) checkConstraintMat(p=p, M=M, restricted=restricted, constraints=constraints)

  # Calculate the prediction
  if(pred_type == "cond_mean") { # Exact conditional mean

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
    pred_ints <- NULL
    pi <- NULL
    pi_type <- "none"
    q_tocalc <- numeric(0)
    mix_weights <- FALSE
    mix_pred <- NULL
    mix_pred_ints <- NULL
  } else { # pred_type != cond_mean: Simulate future values of the process

    # Simulations
    sim <- simulateGSMAR(gsmar, nsimu=n_ahead, initvalues=data, ntimes=nsimu, drop=FALSE)
    sample <- sim$sample
    alpha_mt <- sim$mixing_weights
    colnames(alpha_mt) <- vapply(1:sum(M), function(m) paste("regime", m), character(1))

    # Point forecasts
    myFUN <- ifelse(pred_type == "mean", mean, median)
    pred <- apply(sample, 1, FUN=myFUN)
    mix_pred <- apply(alpha_mt, MARGIN=1:2, FUN=mean)

    # Prediction intervals
    if(pi_type == "upper") {
      q_tocalc <- pi
    } else if(pi_type == "lower") {
      q_tocalc <- 1 - pi
    } else if(pi_type == "two-sided") {
      lower <- (1 - pi)/2
      upper <- rev(1 - lower)
      q_tocalc <- c(lower, upper)
    } else {  # If pi_type == "none"
      q_tocalc <- numeric(0)
      pi <- NULL
    }
    q_tocalc <- sort(q_tocalc, decreasing=FALSE)
    pred_ints <- apply(sample, 1, FUN=quantile, probs=q_tocalc)
    mix_pred_ints <- apply(alpha_mt, MARGIN=1:2, FUN=quantile, probs=q_tocalc)

    if(pi_type != "none") {
      if(length(q_tocalc) == 1) {
        pred_ints <- as.matrix(pred_ints)
        mix_pred_ints <- array(mix_pred_ints, dim=c(n_ahead, sum(M), length(q_tocalc)), dimnames=list(NULL, colnames(alpha_mt), q_tocalc))
        mix_pred_ints <- aperm(mix_pred_ints, perm=c(1, 3, 2))
      } else {
        pred_ints <- t(pred_ints)
        mix_pred_ints <- aperm(mix_pred_ints, perm=c(2, 1, 3)) # So that for each [, , i1] the dimensions match with point forecasts
      }
      colnames(pred_ints) <- q_tocalc
    }
  }

  ret <- structure(list(gsmar=gsmar,
                        pred=pred,
                        pred_ints=pred_ints,
                        mix_pred=mix_pred,
                        mix_pred_ints=mix_pred_ints,
                        n_ahead=n_ahead,
                        nsimu=nsimu,
                        pi=pi,
                        pi_type=pi_type,
                        pred_type=pred_type,
                        q=q_tocalc,
                        mix_weights=mix_weights),
                   class="gsmarpred")
  if(plotRes) plot.gsmarpred(x=ret, nt=nt, mix_weights=mix_weights, ...)
  ret
}


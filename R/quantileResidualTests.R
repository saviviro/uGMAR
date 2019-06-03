#' @import stats
#'
#' @title Quantile residual tests for GMAR, StMAR or G-StMAR model
#'
#' @description \code{quantileResidualTests} performs quantile residual tests for GMAR, StMAR or G-StMAR model,
#'  testing normality, autocorrelation and conditional heteroscedasticy.
#'
#' @inheritParams simulateGSMAR
#' @param lagsAC a numeric vector of positive integers specifying the lags for which autocorrelation is tested.
#' @param lagsCH a numeric vector of positive integers specifying the lags for which conditional heteroscedasticy is tested.
#' @param nsimu a positive integer specifying to how many simulated observations the covariance matrix Omega
#'  (see Kalliovirta (2012)) should be based on. If smaller than data size, then omega will be based on the given data.
#' @param printRes a logical argument defining whether results should be printed or not.
#' @details For details about the quantile residual tests see the cited article by \emph{Kalliovirta (2012)}.
#' @return Returns an object of class \code{'qrtest'} containing the test results in data frames. In the cases
#'   of autocorrelation and conditional heteroscedasticity tests, the returned object also contains the
#'   associated individual statistics and their standard errors, discussed by \emph{Kalliovirta (2012)} at
#'   the pages 369-370.
#' @section Suggested packages:
#'   Install the suggested package "gsl" for faster evaluations in the cases of StMAR and G-StMAR models.
#'   For large StMAR and G-StMAR models with large data the evaluations may take significantly long time without
#'   the package "gsl".
#' @inherit quantileResiduals_int references
#' @seealso \code{\link{fitGSMAR}}, \code{\link{GSMAR}}, \code{\link{diagnosticPlot}}, \code{\link{predict.gsmar}},
#'  \code{\link{getOmega}},
#' @examples
#' \donttest{
#' # GMAR model
#' params12 <- c(1.12, 0.91, 0.29, 4.53, 0.70, 3.21, 0.84)
#' gmar12 <- GSMAR(VIX, 1, 2, params12)
#' qrtest12 <- quantileResidualTests(gmar12)
#' plot(qrtest12)
#'
#' # Restricted GMAR model
#' params12r <- c(1.4, 1.8, 0.88, 0.29, 3.18, 0.84)
#' gmar12r <- GSMAR(data=VIX, p=1, M=2, params=params12r, model="GMAR",
#'  restricted=TRUE)
#' qrtest12r <- quantileResidualTests(gmar12r, lagsAC=1:10, nsimu=1)
#' plot(qrtest12r)
#'
#' # StMAR model
#' params12t <- c(1.38, 0.88, 0.27, 3.8, 0.74, 3.15, 0.8, 100, 3.6)
#' stmar12 <- GSMAR(data=VIX, p=1, M=2, params=params12t, model="StMAR")
#' quantileResidualTests(stmar12, lagsAC=c(1, 2, 5), nsimu=1, printRes=FALSE)
#'
#' # G-StMAR model
#' params12gs <- c(1.38, 0.88, 0.27, 3.8, 0.74, 3.15, 0.8, 3.6)
#' gstmar12 <- GSMAR(data=VIX, p=1, M=c(1, 1), params=params12gs,
#'  model="G-StMAR")
#' quantileResidualTests(gstmar12, lagsAC=c(1, 3), lagsCH=1:2,
#'  nsimu=1, printRes=FALSE)
#'
#' # Such StMAR(3,2) that the AR coefficients are restricted to be
#' # the same for both regimes and that the second AR coefficients are
#' # constrained to zero.
#' params32trc <- c(2.2, 1.8, 0.88, -0.03, 2.4, 0.27, 0.40, 3.9, 1000)
#' stmar32rc <- GSMAR(data=VIX, p=3, M=2, params=params32trc, model="StMAR",
#'  restricted=TRUE, constraints=matrix(c(1, 0, 0, 0, 0, 1), ncol=2))
#' quantileResidualTests(stmar32rc, lagsAC=c(1, 3), nsimu=1, printRes=FALSE)
#' }
#' @export

quantileResidualTests <- function(gsmar, lagsAC=c(1, 2, 5, 10), lagsCH=lagsAC, nsimu=2000, printRes=TRUE) {
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
  T_obs <- length(data) - p

  if(!requireNamespace("gsl", quietly=TRUE) & model %in% c("StMAR", "G-StMAR")) {
    message("Suggested package 'gsl' not found and StMAR or G-StMAR model is considered: performing the tests may take a while")
  }
  if(max(c(lagsAC, lagsCH)) >= T_obs) stop("The lags are too large compared to the data size")
  qresiduals <- quantileResiduals_int(data=data, p=p, M=M, params=params, model=model, restricted=restricted,
                                      constraints=constraints, parametrization=parametrization)
  if(nsimu > length(data)) {
    omegaData <- as.matrix(simulateGSMAR(gsmar, nsimu=nsimu)$sample)
  } else {
    omegaData <- data
  }

  try_to_get_omega <- function(g, dim_g, which_test, which_lag=NA) {
    print_message <- function(because_of) {
      if(which_test == "norm") {
        message(paste("Can't perform normality test", because_of))
      } else if(which_test == "ac") {
        message(paste("Can't perform autocorrelation test for lag", which_lag, because_of))
      } else if(which_test == "ch") {
        message(paste("Can't perform conditional heteroskedasticity test for lag", which_lag, because_of))
      }
    }
    num_string <- "because of numerical problems."
    omg <- tryCatch(getOmega(data=omegaData, p=p, M=M, params=params, model=model, restricted=restricted,
                             constraints=constraints, parametrization=parametrization, g=g, dim_g=dim_g),
                    error=function(e) {
                      if(model == "StMAR") {
                        dfs <- pick_dfs(p=p, M=M, params=params, model=model)
                        if(any(dfs > 30)) {
                          print_message("- possibly because some degrees of freedom parameter is very large. Consider estimating a G-StMAR model.")
                        } else {
                          print_message(num_string)
                        }
                      } else if(model == "G-StMAR") {
                        dfs <- pick_dfs(p=p, M=M, params=params, model=model)
                        if(any(dfs > 30)) {
                          print_message("- possibly because some degrees of freedom parameter is very large. Consider changing one StMAR-type component into GMAR-type and re-estimate.")
                        }
                      } else {
                        print_message(num_string)
                      }
                      return(NA)
                    })
    if(is.matrix(omg) & anyNA(omg)) {
      print_message("- possibly because the model fits too poorly")
    } else if(length(omg) == 1) {
      if(is.na(omg)) return(matrix(NA, nrow=dim_g, ncol=dim_g))
    }
    omg
  }

  format_value0 <- format_valuef(0)
  format_value3 <- format_valuef(3)
  print_resf <- function(lag, p_val) {
    if(lag < 10) {
      cat(" ", format_value0(lag), " | ", format_value3(p_val), "\n")
    } else {
      cat(" ", format_value0(lag), "| ", format_value3(p_val), "\n")
    }
  }

  ####################
  ## Test normality ## (Kalliovirta 2012 sec. 3.3)
  ####################

  g <- function(r) {
    cbind(r^2 - 1, r^3, r^4 - 3)
  }
  dim_g <- 3

  # Omega (Kalliovirta 2013 eq.(2.4))
  Omega <- try_to_get_omega(g=g, dim_g=dim_g, which_test="norm", which_lag=NA)

  # Test statistics and p-value
  sumg <- as.matrix(rowSums(t(g(qresiduals))))
  N <- crossprod(sumg, solve(Omega, sumg))/T_obs
  pvalue <- 1 - pchisq(N, df=dim_g)

  # Results
  if(printRes == TRUE) cat(paste0("Normality test p-value: ", format_value3(pvalue)), "\n\n")
  norm_res <- data.frame(testStat=N, df=dim_g, pvalue=pvalue, row.names=NULL)

  ##########################
  ## Test autocorrelation ## (Kalliovirta 2012 sec. 3.1)
  ##########################
  tmp <- rep(NA, length(lagsAC))
  ac_res <- data.frame(lags=lagsAC, testStat=tmp, df=tmp, pvalue=tmp, indStat=tmp, stdError=tmp)

  # Calculate autocorrelations  FUN = r[i1]*r[i1 - i2]
  get_g <- function(lag, FUN) {
    FUN <- match.fun(FUN)
    function(r) {
      res <- vapply((1 + lag):length(r), function(i1) vapply(1:lag, function(i2) FUN(r, i1, i2), numeric(1)), numeric(lag))
      if(lag == 1) {
        return(as.matrix(res))
      } else {
        return(t(res))
      }
    }
  } # Returns (T - lag x lag) matrix (lag = dim_g)

  if(printRes == TRUE) cat("Autocorrelation tests:\nlags | p-value\n")
  j1 <- 1
  for(lag in lagsAC) {
    g <- get_g(lag, FUN=function(r, i1, i2) r[i1]*r[i1 - i2])

    # Omega (Kalliovirta 2013 eq.(2.4))
    Omega <- try_to_get_omega(g=g, dim_g=lag, which_test="ac", which_lag=lag)

    # Test statistics, sample autocorrelation for of the current lag and p-value
    sumg <- as.matrix(colSums(g(qresiduals)))  # Unscaled and uncentered sample autocovariances of the quantile residuals
    A <- crossprod(sumg, solve(Omega, sumg))/(T_obs - lag)
    sampleAC <- sumg[lag]/(T_obs - lag)
    stdError <- sqrt(Omega[lag, lag]/T_obs)
    pvalue <- 1 - pchisq(A, df=lag)

    # Results
    if(printRes == TRUE) print_resf(lag=lag, p_val=pvalue)
    ac_res$testStat[j1] <- A
    ac_res$df[j1] <- lag
    ac_res$pvalue[j1] <- pvalue
    ac_res$indStat[j1] <- sampleAC
    ac_res$stdError[j1] <- stdError
    j1 <- j1 + 1
  }

  #########################################
  ## Test conditional heteroscedasticity ## (Kalliovirta 2012 sec. 3.2)
  #########################################
  tmp <- rep(NA, length(lagsCH))
  ch_res <- data.frame(lags=lagsCH, testStat=tmp, df=tmp, pvalue=tmp, indStat=tmp, stdError=tmp)

  # Calculate c-heteroskedasticity statistics FUN = (r[i1]^2 - 1)*r[i1 - i2]^2
  if(printRes == TRUE) cat("\nConditional heteroskedasticity tests:\nlags | p-value\n")
  j1 <- 1
  for(lag in lagsCH) {
    g <- get_g(lag, FUN=function(r, i1, i2) (r[i1]^2 - 1)*r[i1 - i2]^2)

    # Omega (Kalliovirta 2012 eq.(2.4))
    Omega <- try_to_get_omega(g=g, dim_g=lag, which_test="ch", which_lag=lag)

    # Test statistics,individual statisics d_k and p-value
    sumg <- as.matrix(colSums(g(qresiduals)))
    H <- crossprod(sumg, solve(Omega, sumg))/(T_obs - lag)
    indStat <- sumg[lag]/(T_obs - lag)
    stdError <- sqrt(Omega[lag, lag]/T_obs)
    pvalue <- 1 - pchisq(H, df=lag)

    # Results
    if(printRes == TRUE) print_resf(lag=lag, p_val=pvalue)
    ch_res$testStat[j1] <- H
    ch_res$df[j1] <- lag
    ch_res$pvalue[j1] <- pvalue
    ch_res$indStat[j1] <- indStat
    ch_res$stdError[j1] <- stdError
    j1 = j1 + 1
  }

  structure(list(norm_res=norm_res,
                 ac_res=ac_res,
                 ch_res=ch_res),
            class="qrtest")
}



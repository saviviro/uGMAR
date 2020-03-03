#' @import stats
#'
#' @title Quantile residual tests for GMAR, StMAR , and G-StMAR models
#'
#' @description \code{quantileResidualTests} performs quantile residual tests for GMAR, StMAR,
#'  and G-StMAR models, testing normality, autocorrelation, and conditional heteroscedasticity
#'  of the quantile residuals.
#'
#' @inheritParams simulateGSMAR
#' @param lagsAC a numeric vector of positive integers specifying the lags for which autocorrelation is tested.
#' @param lagsCH a numeric vector of positive integers specifying the lags for which conditional heteroscedasticity
#'  is tested.
#' @param nsimu a positive integer specifying to how many simulated observations the covariance matrix Omega
#'  (see Kalliovirta (2012)) should be based on. If smaller than data size, then omega will be based on the
#'  given data and not on simulated data. Having the covariance matrix omega based on a large simulated sample
#'  might improve the tests size properties.
#' @param printRes a logical argument defining whether the results should be printed or not.
#' @details
#'   For a correctly specified GSMAR model employing the maximum likelihood estimator, the quantile residuals
#'   are asymptotically independent with standard normal distribution. They can hence be used in a similar
#'   manner to conventional Pearson's residuals. For more details about quantile residual based diagnostics,
#'   and in particular, about the quantile residual tests, see the cited article by \emph{Kalliovirta (2012)}.
#' @return Returns an object of class \code{'qrtest'} containing the test results in data frames. In the cases
#'   of autocorrelation and conditional heteroscedasticity tests, the returned object also contains the
#'   associated individual statistics and their standard errors, discussed in \emph{Kalliovirta (2012)} at
#'   the pages 369-370.
#' @section Suggested packages:
#'   Install the suggested package "gsl" for faster evaluations in the cases of StMAR and G-StMAR models.
#'   For large StMAR and G-StMAR models with large data, the evaluations may take significantly long time
#'   without the package "gsl".
#' @inherit quantileResiduals_int references
#' @seealso \code{\link{profile_logliks}}, \code{\link{fitGSMAR}}, \code{\link{GSMAR}}, \code{\link{diagnosticPlot}},
#'  \code{\link{predict.gsmar}}, \code{\link{getOmega}},
#' @examples
#' \donttest{
#' # GMAR model
#' fit12 <- fitGSMAR(data=logVIX, p=1, M=2, model="GMAR")
#' qrtest12 <- quantileResidualTests(fit12, nsimu=1)
#' plot(qrtest12)
#'
#' # Restricted GMAR model
#' fit12r <- fitGSMAR(logVIX, 1, 2, model="GMAR", restricted=TRUE)
#' qrtest12r <- quantileResidualTests(fit12r, lagsAC=1:10, nsimu=1)
#' plot(qrtest12r)
#'
#' # Non-mixture version of StMAR model
#' fit11t <- fitGSMAR(logVIX, 1, 1, model="StMAR", ncores=1, ncalls=1)
#' quantileResidualTests(fit11t, lagsAC=c(1, 2, 5), nsimu=1, printRes=FALSE)
#'
#' # G-StMAR model
#' fit12gs <- fitGSMAR(logVIX, 1, M=c(1, 1), model="G-StMAR")
#' quantileResidualTests(fit12gs, lagsAC=c(1, 3), lagsCH=1:2,
#'  nsimu=1, printRes=FALSE)
#'
#' # GMAR(p=2, M=2) model such that the second AR coefficient of the
#' # second regime is constrained to zero.
#' constraints <- list(diag(1, ncol=2, nrow=2), as.matrix(c(1, 0)))
#' fit22c <- fitGSMAR(logVIX, 2, 2, constraints=constraints)
#' quantileResidualTests(fit22c, lagsAC=c(1, 3), nsimu=1, printRes=FALSE)
#' }
#' @export

quantileResidualTests <- function(gsmar, lagsAC=c(1, 2, 5, 10), lagsCH=lagsAC, nsimu=1, printRes=TRUE) {
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
    message("Suggested package 'gsl' not found and a StMAR or G-StMAR model is considered: performing the tests may take a while")
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
                      if(model %in% c("StMAR", "G-StMAR")) {
                        dfs <- pick_dfs(p=p, M=M, params=params, model=model)
                        if(any(dfs > 100)) {
                          print_message("- possibly because some degrees of freedom parameter is very large. Consider changing the corresponding StMAR type regimes into GMAR type with the function 'stmar_to_gstmar'.")
                        } else {
                          print_message(num_string)
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
    sep <- ifelse(lag < 10, " | ", "| ")
    cat(" ", format_value0(lag), sep, format_value3(p_val), "\n")
  }

  ####################
  ## Test normality ## (Kalliovirta 2012 sec. 3.3)
  ####################

  g <- function(r) {
    cbind(r^2 - 1, r^3, r^4 - 3)
  }
  dim_g <- 3

  # Omega (Kalliovirta 2012, Lemma 2.2)
  Omega <- try_to_get_omega(g=g, dim_g=dim_g, which_test="norm", which_lag=NA)

  # Test statistics and p-value
  sumg <- as.matrix(colSums(g(qresiduals)))
  N <- crossprod(sumg, solve(Omega, sumg))/T_obs
  pvalue <- 1 - pchisq(N, df=dim_g)

  # Results
  if(printRes) cat(paste0("Normality test p-value: ", format_value3(pvalue)), "\n\n")
  norm_res <- data.frame(testStat=N, df=dim_g, pvalue=pvalue, row.names=NULL)

  #############################################################
  ## Test autocorrelation and conditional heteroscedasticity ## (Kalliovirta 2012 sec. 3.1 and 3.2)
  #############################################################

  # Storage for the results
  tmp <- rep(NA, length(lagsAC))
  ac_res <- data.frame(lags=lagsAC, testStat=tmp, df=tmp, pvalue=tmp, indStat=tmp, stdError=tmp)
  tmp <- rep(NA, length(lagsCH))
  ch_res <- data.frame(lags=lagsCH, testStat=tmp, df=tmp, pvalue=tmp, indStat=tmp, stdError=tmp)
  ret <- list(norm_res=norm_res,
              ac_res=ac_res,
              ch_res=ch_res)

  # Returns the function 'g' for a given lag and function FUN to be applied.
  get_g <- function(lag, FUN) {
    FUN <- match.fun(FUN)
    function(r) {  # Takes in quantile residuals vector r, returns a (T - lag x lag) matrix. (lag = dim_g)
      res <- vapply((1 + lag):length(r), function(i1) vapply(1:lag, function(i2) FUN(r, i1, i2), numeric(1)), numeric(lag))
      if(lag == 1) {
        return(as.matrix(res))
      } else {
        return(t(res))
      }
    }
  }

  # Apart from the function 'g', the test statistics are similar for AC and CH tests
  test_ac_or_ch <- function(which_test=c("ac", "ch")) {
    which_test <- match.arg(which_test)

    # Which lags to go through
    if(which_test == "ac") {
      lags_to_loop <- lagsAC
    } else {
      lags_to_loop <- lagsCH
    }

    j1 <- 1 # Count iterations
    for(lag in lags_to_loop) {

      # The function 'g'
      if(which_test == "ac") {
        g <- get_g(lag, FUN=function(r, i1, i2) r[i1]*r[i1 - i2]) # FUN = r[i1]*r[i1 - i2]
      } else { # to_test == "ch"
        g <- get_g(lag, FUN=function(r, i1, i2) (r[i1]^2 - 1)*r[i1 - i2]^2) # FUN = (r[i1]^2 - 1)*r[i1 - i2]^2
      }

      # Omega (Kalliovirta 2013 eq.(2.4))
      Omega <- try_to_get_omega(g=g, dim_g=lag, which_test=which_test, which_lag=lag)

      # Test statistics: sample autocorrelation c_k/h.sked statistic d_k for of the current lag, standard error, and p-value
      sumg <- as.matrix(colSums(g(qresiduals)))
      AorH <- crossprod(sumg, solve(Omega, sumg))/(T_obs - lag) # See A ("ac") and H ("ch") test statistics in Kalliovirta 2012, pp.369-370
      indStat <- sumg[lag]/(T_obs - lag) # c_k ("ac") or d_k ("ch")
      stdError <- sqrt(Omega[lag, lag]/T_obs) # See Kalliovirta 2012, pp.369-370
      pvalue <- 1 - pchisq(AorH, df=lag)

      # Store the results
      if(printRes) print_resf(lag=lag, p_val=pvalue)
      index <- ifelse(which_test == "ac", 2, 3) # Index in the list 'ret'
      ret[[index]]$testStat[j1] <<- AorH
      ret[[index]]$df[j1] <<- lag
      ret[[index]]$pvalue[j1] <<- pvalue
      ret[[index]]$indStat[j1] <<- indStat
      ret[[index]]$stdError[j1] <<- stdError
      j1 <- j1 + 1
    }
  }

  if(printRes) cat("Autocorrelation tests:\nlags | p-value\n")
  test_ac_or_ch("ac")

  if(printRes) cat("\nConditional heteroskedasticity tests:\nlags | p-value\n")
  test_ac_or_ch("ch")

  structure(ret, class="qrtest")
}



#' @import stats
#'
#' @title Quantile residual tests for GMAR, StMAR , and G-StMAR models
#'
#' @description \code{quantile_residual_tests} performs quantile residual tests for GMAR, StMAR,
#'  and G-StMAR models, testing normality, autocorrelation, and conditional heteroscedasticity
#'  of the quantile residuals.
#'
#' @inheritParams add_data
#' @param lags_ac a numeric vector of positive integers specifying the lags for which autocorrelation is tested.
#' @param lags_ch a numeric vector of positive integers specifying the lags for which conditional heteroscedasticity
#'  is tested.
#' @param nsimu a positive integer specifying to how many simulated observations the covariance matrix Omega
#'  (see Kalliovirta (2012)) should be based on. If smaller than data size, then omega will be based on the
#'  given data and not on simulated data. Having the covariance matrix omega based on a large simulated sample
#'  might improve the tests size properties.
#' @param print_res a logical argument defining whether the results should be printed or not.
#' @details
#'   For a correctly specified GSMAR model employing the maximum likelihood estimator, the quantile residuals
#'   are asymptotically independent with standard normal distribution. They can hence be used in a similar
#'   manner to conventional Pearson's residuals. For more details about quantile residual based diagnostics,
#'   and in particular, about the quantile residual tests, see the cited article by \emph{Kalliovirta (2012)}.
#' @return Returns an object of class \code{'qrtest'} containing the test results in data frames. In the cases
#'   of autocorrelation and conditional heteroscedasticity tests, the returned object also contains the
#'   associated individual statistics and their standard errors, discussed in \emph{Kalliovirta (2012)} at
#'   the pages 369-370.
#' @inherit quantile_residuals_int references
#' @seealso \code{\link{profile_logliks}}, \code{\link{fitGSMAR}}, \code{\link{GSMAR}}, \code{\link{diagnostic_plot}},
#'  \code{\link{predict.gsmar}}, \code{\link{get_test_Omega}},
#' @examples
#' \donttest{
#' ## The below examples take approximately 30 seconds to run.
#'
#' # G-StMAR model with one GMAR type and one StMAR type regime
#' fit42gs <- fitGSMAR(data=M10Y1Y, p=4, M=c(1, 1), model="G-StMAR",
#'                     ncalls=1, seeds=4)
#'
#' # Tests based on the observed data (without simulation procedure) with the
#' # default lags:
#' qrt1 <- quantile_residual_tests(fit42gs)
#'
#' # Tests based on the simulation procedure using sample size 10000 and with
#' # the lags specified by hand:
#' set.seed(1)
#' qrt2 <- quantile_residual_tests(fit42gs, lags_ac=c(1, 6), nsimu=10000)
#'
#' # GMAR model
#' fit12 <- fitGSMAR(data=simudata, p=1, M=2, model="GMAR", ncalls=1, seeds=1)
#' qrt3 <- quantile_residual_tests(fit12, lags_ac=c(1, 5, 10, 15))
#' }
#' @export

quantile_residual_tests <- function(gsmar, lags_ac=c(1, 3, 6, 12), lags_ch=lags_ac, nsimu=1, print_res=TRUE) {

  # Checks + collect the relevant statistics etc
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
  if(max(c(lags_ac, lags_ch)) >= T_obs) stop("The lags are too large compared to the data size")

  # Obtain the quantile residuals
  qresiduals <- quantile_residuals_int(data=data, p=p, M=M, params=params, model=model, restricted=restricted,
                                       constraints=constraints, parametrization=parametrization)

  # Sample used in calculation of Omega: either simulated or the data
  if(nsimu > length(data)) {
    omegaData <- as.matrix(simulate.gsmar(gsmar, nsim=nsimu)$sample)
  } else {
    omegaData <- data
  }

  # Function to calculate the Omega matrix (that is presented in Kalliovirta 2012, Lemma 2.2); calls the function get_test_Omega
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
    omg <- tryCatch(get_test_Omega(data=omegaData, p=p, M=M, params=params, model=model, restricted=restricted,
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

  # Functions to format values and print the results
  format_value0 <- format_valuef(0)
  format_value3 <- format_valuef(3)
  print_resf <- function(lag, p_val) {
    sep <- ifelse(lag < 10, " | ", "| ")
    cat(" ", format_value0(lag), sep, format_value3(p_val), "\n")
  }

  ####################
  ## Test normality ## (Kalliovirta 2012, Section 3.3)
  ####################

  # The function 'g' of Kalliovirta 2012, Section 3.3
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
  if(print_res) cat(paste0("Normality test p-value: ", format_value3(pvalue)), "\n\n")
  norm_res <- data.frame(testStat=N, df=dim_g, pvalue=pvalue, row.names=NULL)

  #############################################################
  ## Test autocorrelation and conditional heteroscedasticity ## (Kalliovirta 2012, Sections 3.1 and 3.2)
  #############################################################

  # Storage for the results
  tmp <- rep(NA, length(lags_ac))
  ac_res <- data.frame(lags=lags_ac, testStat=tmp, df=tmp, pvalue=tmp, indStat=tmp, stdError=tmp)
  tmp <- rep(NA, length(lags_ch))
  ch_res <- data.frame(lags=lags_ch, testStat=tmp, df=tmp, pvalue=tmp, indStat=tmp, stdError=tmp)
  ret <- list(norm_res=norm_res,
              ac_res=ac_res,
              ch_res=ch_res)

  # Returns the function 'g' for a given lag and function FUN to be applied (see Kalliovirta 2012, Sections 3.1 and 3.2 for the 'g').
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

  # Apart from the function 'g', the test statistics are similar for AC and CH tests.
  # Function to calculate ac and ch tests
  test_ac_or_ch <- function(which_test=c("ac", "ch")) {
    which_test <- match.arg(which_test)

    # Which lags to go through
    if(which_test == "ac") {
      lags_to_loop <- lags_ac
    } else {
      lags_to_loop <- lags_ch
    }

    j1 <- 1 # Count iterations
    for(lag in lags_to_loop) {

      # The function 'g'
      if(which_test == "ac") {
        g <- get_g(lag, FUN=function(r, i1, i2) r[i1]*r[i1 - i2]) # FUN = r[i1]*r[i1 - i2]
      } else { # to_test == "ch"
        g <- get_g(lag, FUN=function(r, i1, i2) (r[i1]^2 - 1)*r[i1 - i2]^2) # FUN = (r[i1]^2 - 1)*r[i1 - i2]^2
      }

      # Omega (Kalliovirta 2012, Lemma 2.2)
      Omega <- try_to_get_omega(g=g, dim_g=lag, which_test=which_test, which_lag=lag)

      # Test statistics: sample autocorrelation c_k/h.sked statistic d_k for of the current lag, standard error, and p-value
      sumg <- as.matrix(colSums(g(qresiduals)))
      AorH <- crossprod(sumg, solve(Omega, sumg))/(T_obs - lag) # See A ("ac") and H ("ch") test statistics in Kalliovirta 2012, pp.369-370
      indStat <- sumg[lag]/(T_obs - lag) # c_k ("ac") or d_k ("ch") of Kalliovirta 2012
      stdError <- sqrt(Omega[lag, lag]/T_obs) # See Kalliovirta 2012, pp.369-370
      pvalue <- 1 - pchisq(AorH, df=lag)

      # Store the results
      if(print_res) print_resf(lag=lag, p_val=pvalue)
      index <- ifelse(which_test == "ac", 2, 3) # Index in the list 'ret'
      ret[[index]]$testStat[j1] <<- AorH
      ret[[index]]$df[j1] <<- lag
      ret[[index]]$pvalue[j1] <<- pvalue
      ret[[index]]$indStat[j1] <<- indStat
      ret[[index]]$stdError[j1] <<- stdError
      j1 <- j1 + 1
    }
  }

  # Calculate the tests and print the results
  if(print_res) cat("Autocorrelation tests:\nlags | p-value\n")
  test_ac_or_ch("ac")

  if(print_res) cat("\nConditional heteroskedasticity tests:\nlags | p-value\n")
  test_ac_or_ch("ch")

  structure(ret, class="qrtest")
}



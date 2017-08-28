#' @import stats
#'
#' @title Quantile residual tests for GMAR or StMAR model
#'
#' @description \code{quantileResidualTests} performs quantile residual tests for GMAR or StMAR model, testing normality, autocorrelation and conditional heteroscedasticy.
#'
#' @inheritParams loglikelihood
#' @param lagsAC an (optional) numeric vector of positive integers specifying the lags for which autocorrelation is tested. Default is \code{c(1, 2, 4, 6, 8, 10)}.
#' @param lagsCH an (optional) numeric vector of positive integers specifying the lags for which conditional heteroscedasticy is tested. Default is \code{c(1, 2, 4, 6, 8, 10)}.
#' @param nsimu an (optional) positive integer specifying to how many simulated observations the covariance matrix Omega should be based on. If smaller than data size, then omega will be based on the given data. Default is 2000.
#' @param printRes an (optional) logical argument defining wether results should be printed or not. Default is \code{TRUE}.
#' @details \code{NA} values mean that it was not (numerically) possible for the code to calculate all the necessary estimates for the tests.
#'   This may suggest that the model is misspecified.
#'
#' For details about the quantile residual tests see \emph{Kalliovirta 2012} from the references.
#' @return a list of data frames containing the test results:
#'   \describe{
#'   \item{\code{$normality}}{A data frame containing results from the normality test.
#'     \code{$testStat} containing the test statistic and \code{$pvalue} containing the corresponding p-value.}
#'   \item{\code{$autocorrelation}}{A data frame containing results from the autocorrelation tests.
#'     \code{$testStat} containing the test statistics and \code{$pvalue} containing the corresponding p-values.
#'     \code{$indStat} containing uncentered sample autocovariances of the quantile residuals for each lag,
#'     and \code{$stdError} containing their approximate standard errors.}
#'   \item{\code{$cond.heteroscedasticity}}{A data frame containing results from the conditional heteroskedasticity tests.
#'     \code{$testStat} containing the test statistics and \code{$pvalue} containing the corresponding p-values.
#'     \code{$indStat} containing the individual statistics assiociated with each lag,
#'     and \code{$stdError} containing their approximate standard errors.}
#'   }
#' @section Printed results:
#'   The results from the tests are printed so that the letter "N" means normality test, "A" autocorrelation test
#'   and "H" conditional heteroscedasticity test. The numbers right next to "A" and "H" indicate the number of lags used
#'   in each test. The statistics following them are the corresponding test statistics and p-values.
#' @section Suggested packages:
#'   Install the suggested package "gsl" for faster evaluations in the cases of StMAR models.
#'   For large StMAR models with large data the evaluations may take significantly long time without
#'   the package "gsl".
#' @inherit quantileResiduals_int references
#' @examples
#' \donttest{
#' # GMAR model
#' params13 <- c(1.4, 0.88, 0.26, 2.46, 0.82, 0.74, 5.0, 0.68, 5.2, 0.72, 0.2)
#' tests13 <- quantileResidualTests(VIX, 1, 3, params13)
#'
#' # Restricted GMAR model, using the given data instead of simulated data
#' params12r <- c(1.4, 1.8, 0.88, 0.29, 3.18, 0.84)
#' tests12r <- quantileResidualTests(VIX, 1, 2, params12r, restricted=TRUE, nsimu=1)
#'
#' # StMAR model
#' params12t <- c(1.38, 0.88, 0.27, 3.8, 0.74, 3.15, 0.8, 120, 3.6)
#' tests12t <- quantileResidualTests(VIX, 1, 2, params12t, StMAR=TRUE,
#'                                   lagsAC=c(1, 5, 10), lagsCH=c(1, 5, 10))
#'
#' # GMAR model as a mixture of AR(2) and AR(1) models
#' R <- list(diag(1, ncol=2, nrow=2), as.matrix(c(1, 0)))
#' params22c <- c(1.2, 0.85, 0.04, 0.3, 3.3, 0.77, 2.8, 0.77)
#' tests22c <- quantileResidualTests(VIX, 2, 2, params22c, constraints=TRUE, R=R)
#'
#' # Such StMAR(3,2) that the AR coefficients are restricted to be
#' # the same for both regimes and that the second AR coefficients are
#' # constrained to zero.
#' params32trc <- c(2.2, 1.8, 0.88, -0.03, 2.4, 0.27, 0.40, 3.9, 1000)
#' qr32trc <- quantileResidualTests(VIX, 3, 2, params32trc, StMAR=TRUE,
#'                                   restricted=TRUE, constraints=TRUE,
#'                                   R=matrix(c(1, 0, 0, 0, 0, 1), ncol=2),
#'                                   lagsAC=c(1, 5, 10), lagsCH=c(1, 5, 10))
#' }
#' @export

quantileResidualTests <- function(data, p, M, params, StMAR=FALSE, restricted=FALSE, constraints=FALSE, R, lagsAC=c(1, 2, 4, 6, 8, 10), lagsCH=c(1, 2, 4, 6, 8, 10), nsimu=2000, printRes=TRUE) {

  checkPM(p, M)
  if(length(params)!=nParams(p=p, M=M, StMAR=StMAR, restricted=restricted, constraints=constraints, R=R)) {
    stop("the parameter vector is wrong dimension")
  }
  data = checkAndCorrectData(data, p)
  T0 = length(data)-p
  if(max(c(lagsAC, lagsCH))>=T0) {
    stop("the lags are too large compared to the data size")
  }
  qresiduals = quantileResiduals_int(data, p, M, params, StMAR=StMAR, restricted=restricted, constraints=constraints, R=R)
  if(nsimu > T0) {
    simuData = simulateGMAR(p, M, params, StMAR=StMAR, restricted=restricted, constraints=constraints, R=R, nsimu=nsimu+200)
    simuData = as.matrix(simuData$sample[201:length(simuData$sample)])
  } else {
    simuData = data
  }
  qresiduals_simuData = quantileResiduals_int(simuData, p, M, params, StMAR=StMAR, restricted=restricted, constraints=constraints, R=R)
  results = list()

  ####################
  ## Test normality ## (Kalliovirta 2012 sec. 3.3)
  ####################

  g <- function(r) {
    cbind(r^2-1, r^3, r^4-3)
  }

  # Omega (Kalliovirta 2013 eq.(2.4))
  Omega = getOmega(simuData, p, M, params, StMAR=StMAR, restricted=restricted, constraints=constraints, R=R, g=g, dim_g=3, qresiduals=qresiduals_simuData)

  # Test statistics and p-value
  sumg = as.matrix(rowSums(t(g(qresiduals))))
  N = t(sumg)%*%solve(Omega, sumg)/T0
  p_norm = 1 - pchisq(N, df=3)

  # Results
  if(printRes==TRUE) {
    print(sprintf("N: %.2f, p-value: %.2f", N, p_norm))
  }
  normRes = data.frame(t(c(N, p_norm)), row.names=NULL)
  colnames(normRes) = c("testStat", "pvalue")
  results[[1]] = normRes

  ##########################
  ## Test autocorrelation ## (Kalliovirta 2012 sec. 3.1)
  ##########################
  acRes = matrix(nrow=length(lagsAC), ncol=5)
  acRes[,1] = lagsAC

  # Calculate autocorrelations
  g0 <- function(r, lag) {
    sapply((1+lag):length(r), function(i1) sapply(1:lag, function(i2) r[i1]*r[i1-i2]))
  }
  j1 = 1
  for(lag in lagsAC) {
    g <- function(r) {
      if(lag>1) {
        return(t(g0(r, lag)))
      } else {
        return(as.matrix(g0(r, lag)))
      }
    }

    # Omega (Kalliovirta 2013 eq.(2.4))
    Omega = getOmega(simuData, p, M, params, StMAR=StMAR, restricted=restricted, constraints=constraints, R=R, g=g, dim_g=lag, qresiduals=qresiduals_simuData)

    # Test statistics, sample autocorrelation for of the current lag and p-value
    sumg = as.matrix(colSums(g(qresiduals)))  # Unscaled and uncentered sample autocovariances of the quantile residuals
    A = t(sumg)%*%solve(Omega, sumg)/(T0-lag)
    sampleAC = sumg[lag]/(T0-lag)
    stdError = sqrt(Omega[lag,lag]/T0)
    p_ac = 1 - pchisq(A, df=lag)

    # Results
    if(printRes==TRUE) {
      print(sprintf("A%.0f: %.2f, p-value: %.2f", lag, A, p_ac))
    }
    acRes[j1, 2] = A
    acRes[j1, 3] = p_ac
    acRes[j1, 4] = sampleAC
    acRes[j1, 5] = stdError
    j1 = j1 + 1
  }
  acRes = data.frame(acRes, row.names=NULL)
  colnames(acRes) = c("lag", "testStat", "pvalue", "indStat", "stdError")
  results[[2]] = acRes

  #########################################
  ## Test conditional heteroscedasticity ## (Kalliovirta 2012 sec. 3.2)
  #########################################

  chRes = matrix(nrow=length(lagsCH), ncol=5)
  chRes[,1] = lagsCH

  # Calculate autocorrelations
  g0 <- function(r, lag) {
    sapply((1+lag):length(r), function(i1) sapply(1:lag, function(i2) (r[i1]^2-1)*r[i1-i2]^2))
  }
  j1 = 1
  for(lag in lagsCH) {
    g <- function(r) {
      if(lag>1) {
        return(t(g0(r, lag)))
      } else {
        return(as.matrix(g0(r, lag)))
      }
    }

    # Omega (Kalliovirta 2012 eq.(2.4))
    Omega = getOmega(simuData, p, M, params, StMAR=StMAR, restricted=restricted, constraints=constraints, R=R, g=g, dim_g=lag, qresiduals=qresiduals_simuData)

    # Test statistics,individual statisics d_k and p-value
    sumg = as.matrix(colSums(g(qresiduals)))
    H = t(sumg)%*%solve(Omega, sumg)/(T0-lag)
    indStat = sumg[lag]/(T0-lag)
    stdError = sqrt(Omega[lag,lag]/T0)
    p_ch = 1 - pchisq(H, df=lag)

    # Results
    if(printRes==TRUE) {
      print(sprintf("H%.0f: %.2f, p-value: %.2f", lag, H, p_ch))
    }
    chRes[j1, 2] = H
    chRes[j1, 3] = p_ch
    chRes[j1, 4] = indStat
    chRes[j1, 5] = stdError
    j1 = j1 + 1
  }
  chRes = data.frame(chRes, row.names=NULL)
  colnames(chRes) = c("lag", "testStat", "pvalue", "indStat", "stdError")
  results[[3]] = chRes
  names(results) = c("normality", "autocorrelation", "cond.heteroscedasticity")
  return(results)
}



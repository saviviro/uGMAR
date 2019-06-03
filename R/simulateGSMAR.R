#' @import stats
#'
#' @title Simulate values from GMAR, StMAR or G-StMAR process
#'
#' @description \code{simulateGSMAR} simulates values from the specified GMAR, StMAR or G-StMAR process. Can be used for
#'  forecasting future values of the process.
#'
#' @inheritParams loglikelihood
#' @param gsmar object of class \code{'gsmar'} created with the function \code{fitGSMAR} or \code{GSMAR}.
#' @param nsimu a positive integer specifying how many values will be simulated.
#' @param initvalues a numeric vector with length \code{>=p} specifying the initial values for the simulation. The \strong{last}
#'  element will be used as the initial value for the first lag, the second last element will be initial value for the second lag, etc.
#'  If not specified, initial values will be simulated from the process's stationary distribution.
#' @param ntimes a positive integer specifying how many sets of simulations should be performed.
#' @param drop if \code{TRUE} (default) then the components of the returned list are coerced to the lowest possible dimension, i.e, for
#'  \code{ntimes=1}, \code{$sample} and \code{$component} will be vectors and \code{$mixing_weights} will be matrix.
#' @details The argument \code{ntimes} is intended for forecasting: a GSMAR process can be forecasted by simulating it's possible future values.
#'  One can easily perform a large number simulations and calculate the sample quantiles from the simulated values to obtain prediction
#'  intervals (see the forecasting example).
#' @return If \code{drop==TRUE} and \code{ntimes=1} (default): \code{$sample} and \code{$component} are vectors
#'  and is \code{$mixing_weights} is matrix. Otherwise, returns a list with...
#'   \describe{
#'     \item{\code{$sample}}{a size (\code{nsimu}\eqn{x}\code{ntimes}) matrix containing the simulated values.}
#'     \item{\code{$component}}{a size (\code{nsimu}\eqn{x}\code{ntimes}) matrix containing the information from which component each
#'      value was generated from.}
#'     \item{\code{$mixing_weights}}{a size (\code{nsimu}\eqn{xMx}\code{ntimes}) array containing the mixing weights corresponding to the
#'      sample: the dimension \code{[i, , ]} is the time index, the dimension \code{[, i, ]} indicates the regime, and the dimension
#'      \code{[, , i]} indicates the i:th set of simulations.}
#'   }
#' @seealso \code{\link{fitGSMAR}}, \code{\link{GSMAR}}, \code{\link{predict.gsmar}},
#'  \code{\link{add_data}}, \code{\link{condMoments}}, \code{\link{mixingWeights}}
#' @inherit loglikelihood references
#' @examples
#'  \donttest{
#'  # GMAR model
#'  params13 <- c(1.4, 0.88, 0.26, 2.46, 0.82, 0.74, 5.0, 0.68, 5.2, 0.72, 0.2)
#'  gmar13 <- GSMAR(p=1, M=3, params=params13, model="GMAR")
#'  sim13 <- simulateGSMAR(gmar13, nsimu=500)
#'  ts.plot(sim13$sample)
#'  ts.plot(sim13$component)
#'  ts.plot(sim13$mixing_weights, col=rainbow(3), lty=2)
#'
#'  # FORECASTING EXAMPLE:
#'  # Restricted GMAR model, 10000 sets of simulations with initial values 6 and 6.2.
#'  params22r <- c(1.4, 1.8, 0.8, -0.1, 0.29, 3.18, 0.84)
#'  gmar22r <- GSMAR(p=2, M=2, params=params22r, model="GMAR",
#'   restricted=TRUE)
#'  sim22r <- simulateGSMAR(gmar22r, nsimu=5, initval=c(6, 6.2), ntimes=10000)
#'  apply(sim22r$sample, 1, median) # Point forecast
#'  apply(sim22r$sample, 1, quantile, probs=c(0.025, 0.975)) # 95% interval
#'  apply(sim22r$mixing_weights, MARGIN=1:2, FUN=median) # mix.weight point forecast
#'  apply(sim22r$mixing_weights, MARGIN=1:2, FUN=quantile,
#'   probs=c(0.025, 0.975)) # mix.weight 95% intervals
#'
#'
#'  # G-StMAR model, with initial values
#'  params12gs <- c(1.38, 0.88, 0.27, 3.8, 0.74, 3.15, 0.8, 3.6)
#'  gstmar12 <- GSMAR(p=1, M=c(1, 1), params=params12gs,
#'  model="G-StMAR")
#'  sim12gs <- simulateGSMAR(gstmar12, nsimu=500, initvalues=5:6)
#'  ts.plot(sim12gs$sample)
#'  ts.plot(sim12gs$component)
#'  ts.plot(sim12gs$mixing_weights, col=rainbow(3), lty=2)
#'
#'  # Such StMAR(3,2) that the AR coefficients are restricted to be
#'  # the same for both regimes and that the second AR coefficients are
#'  # constrained to zero.
#'  params32trc <- c(2.2, 1.8, 0.88, -0.03, 2.4, 0.27, 0.40, 3.9, 1000)
#'  stmar32rc <- GSMAR(data=VIX, p=3, M=2, params=params32trc, model="StMAR",
#'   restricted=TRUE, constraints=matrix(c(1, 0, 0, 0, 0, 1), ncol=2))
#'  sim32trc <- simulateGSMAR(stmar32rc, nsimu=500)
#'  ts.plot(sim32trc$sample)
#'  ts.plot(sim32trc$component)
#'  ts.plot(sim32trc$mixing_weights, col=rainbow(3), lty=2)
#'
#'  # Mixture version of Heterogenuous Autoregressive (HAR) model (without data)
#'  paramsHAR <- c(1, 0.1, 0.2, 0.3, 1, 1, 0.15, 0.25, 0.35, 1, 0.55)
#'  r1 = c(1, rep(0, 21)); r2 = c(rep(0.2, 5), rep(0, 17)); r3 = rep(1/22, 22)
#'  R0 = cbind(r1, r2, r3)
#'  mixhar <- GSMAR(p=22, M=2, params=paramsHAR, model="GMAR", constraints=list(R0, R0))
#'  simhar <- simulateGSMAR(mixhar, nsimu=1000)
#'  ts.plot(simhar$sample)
#'  ts.plot(simhar$component)
#'  ts.plot(simhar$mixing_weights, col=rainbow(3), lty=2)
#' }
#' @export

simulateGSMAR <- function(gsmar, nsimu, initvalues, ntimes=1, drop=TRUE) {
  epsilon <- round(log(.Machine$double.xmin) + 10)
  check_gsmar(gsmar)
  p <- gsmar$model$p
  M <- gsmar$model$M
  model <- gsmar$model$model
  M_orig <- M
  if(model == "G-StMAR") {
    M1 <- M[1]
    M2 <- M[2]
    M <- sum(M)
  }

  if(!missing(initvalues)) {
    if(length(initvalues) < p) {
      stop("The length of initial values vector has to be at least p")
    }
  }
  if(nsimu < 1) {
    stop("The number of simulations has to be equal or larger than one")
  } else if(nsimu%%1 != 0) {
    stop("Argument nsimu has to be positive integer")
  } else if(ntimes < 1 | ntimes%%1 != 0) {
    stop("The argument ntimes should be a positive integer")
  }

  # Reform and collect parameters
  params <- removeAllConstraints(p=p, M=M_orig, params=gsmar$params, model=model, restricted=gsmar$model$restricted,
                                 constraints=gsmar$model$constraints)
  if(gsmar$model$parametrization == "mean") {
    params <- change_parametrization(p=p, M=M_orig, params=params, model=model, restricted=FALSE, constraints=NULL, change_to="intercept")
  }
  pars <- pick_pars(p=p, M=M_orig, params=params, model=model, restricted=FALSE, constraints=NULL)
  alphas <- pick_alphas(p=p, M=M_orig, params=params, model=model, restricted=FALSE, constraints=NULL)
  dfs <- pick_dfs(p=p, M=M_orig, params=params, model=model)
  sigmas <- pars[p + 2,]
  parameterChecks(p=p, M=M_orig, params=params, model=model, restricted=FALSE, constraints=NULL)

  # Create container for the simulated values and initial values. First row row initival values vector, and t+1:th row for (Y_t,Y_[t-1],...,Y_[t-p+1])
  Y <- matrix(nrow=nsimu + 1, ncol=p)

  # Calculate inverses of covariance matrices Gamma_m and their determinants
  invG <- array(dim=c(p, p, M))
  detG <- numeric(M)
  if(p == 1) { # Inverse formula by Galbraith, R., Galbraith, J., (1974)
    for(i1 in 1:M) {
      invG[, , i1] <- (1 - pars[p + 1, i1]^2)/sigmas[i1]
      detG[i1] <- 1/invG[, , i1, drop=TRUE]
     }
  } else {
     for(i1 in 1:M) {
      ARcoefs <- pars[2:(p + 1), i1]
      U <- diag(1, p, p)
      V <- diag(ARcoefs[p], p, p)
      for(i2 in 1:(p - 1)) {
        U[(i2 + 1):p, i2] <- -ARcoefs[1:(p - i2)]
        V[(i2 + 1):p, i2] <- rev(ARcoefs[i2:(p - 1)])
      }
      invG[, , i1] <- (crossprod(U, U) - crossprod(V, V))/sigmas[i1]
      detG[i1] <- 1/det(invG[, , i1])
     }
  }

  # Calculate expected values mu_m (Kalliovirta 2015, s.250)
  mu <- vapply(1:M, function(i1) pars[1, i1]/(1 - sum(pars[2:(p + 1), i1])), numeric(1))
  mu_mp <- matrix(vapply(1:M, function(i1) mu[i1]*matrix(rep(1, p)), numeric(p)), ncol=M) # Column for each component

  # If initial values are missing simulate them from the processes stationary distribution
  if(missing(initvalues)) {
    mv_samples <- numeric(p)
    m <- sample.int(M, size=1, prob=alphas)
    if(model == "GMAR") {
      Gamma_m <- solve(invG[, , m])
      L <- t(chol(Gamma_m))
      mv_samples <- mu_mp[,m] + L%*%rnorm(p)
    } else if(model == "StMAR") {
        Gamma_m <- solve(invG[, , m])
        L <- t(chol((dfs[m] - 2)/dfs[m]*Gamma_m))
        Z <- L%*%rnorm(p) # Sample from N(0,K)
        U <- rchisq(p, df=dfs[m]) # Sample from chisq(v)
        mv_samples <- mu_mp[,m] + Z*sqrt(dfs[m]/U) # sample from student's t
    } else { # If model == "G-StMAR
        if(m <= M1) { # p-dimensional Gaussian samples
          Gamma_m <- solve(invG[, , m])
          L <- t(chol(Gamma_m))
          mv_samples = mu_mp[,m] + L%*%rnorm(p)
        } else { # p-dimensional Student samples
          Gamma_m <- solve(invG[, , m])
          L <- t(chol((dfs[m - M1] - 2)/dfs[m - M1]*Gamma_m))
          Z <- L%*%rnorm(p) # Sample from N(0,K)
          U <- rchisq(p, df=dfs[m - M1]) # Sample from chisq(v)
          mv_samples <- mu_mp[,m] + Z*sqrt(dfs[m - M1]/U) # sample from Student's t
        }
    }
    Y[1,] <- mv_samples # Initial values sampled from the stationary distribution
  } else {
    initvalues <- initvalues[(length(initvalues) - p + 1):length(initvalues)]
    Y[1,] <- rev(initvalues)
  }

  # Initialize data structures
  sample <- matrix(nrow=nsimu, ncol=ntimes)
  component <- matrix(nrow=nsimu, ncol=ntimes)
  mixing_weights <- array(dim=c(nsimu, M, ntimes))

  for(j1 in 1:ntimes) {

    ### Start simulation ###
    for(i1 in 1:nsimu) {
      # Calculate log multinormal values (Kalliovirta 2015, s250 eq.(7)) or log student-t values
      matProd <- vapply(1:M, function(i2) crossprod(Y[i1,] - mu_mp[,i2], as.matrix(invG[, , i2]))%*%(Y[i1,] - mu_mp[,i2]), numeric(1))
      if(model == "GMAR") {
        logmv_values <- vapply(1:M, function(i2) -0.5*p*log(2*pi) - 0.5*log(detG[i2]) - 0.5*matProd[i2], numeric(1))
      } else if(model == "StMAR") {
        logmv_values <- vapply(1:M, function(i2) lgamma(0.5*(p + dfs[i2])) - 0.5*p*log(pi) - 0.5*p*log(dfs[i2] - 2) - lgamma(0.5*dfs[i2]) -
                               0.5*log(detG[i2]) - 0.5*(p + dfs[i2])*log(1 + matProd[i2]/(dfs[i2] - 2)), numeric(1))
      } else { # If model == "G-StMAR
        logmv_values <- c(vapply(1:M1, function(i2) - 0.5*p*log(2*pi) - 0.5*log(detG[i2]) - 0.5*matProd[i2], numeric(1)),
                          vapply((M1 + 1):M, function(i2) lgamma(0.5*(p + dfs[i2 - M1])) - 0.5*p*log(pi) - 0.5*p*log(dfs[i2 - M1] - 2) -
                                 lgamma(0.5*dfs[i2 - M1]) - 0.5*log(detG[i2]) - 0.5*(p + dfs[i2 - M1])*log(1 + matProd[i2]/(dfs[i2 - M1] - 2)),
                                 numeric(1)))
      }

      # Calculate the alpha_mt mixing weights (Kalliovirta 2015, s.250 eq.(8)). Close to zero values handled with Brobdingnag
      if(M == 1) {
        alpha_mt <- 1
      } else if(any(logmv_values < epsilon)) {
        numerators <- lapply(1:M, function(i2) alphas[i2]*Brobdingnag::as.brob(exp(1))^logmv_values[i2])
        denominator <- Reduce("+", numerators)
        alpha_mt <- vapply(1:M, function(i2) as.numeric(numerators[[i2]]/denominator), numeric(1))
      } else {
        mv_values <- exp(logmv_values)
        denominator <- sum(alphas*mv_values)
        alpha_mt <- alphas*mv_values/denominator
      }

      # Draw the component and store the values
      m <- sample.int(M, size=1, prob=alpha_mt)
      component[i1, j1] <- m
      mixing_weights[i1, , j1] <- alpha_mt

      # Draw sample and store it
      mu_mt <- pars[1, m] + sum(pars[2:(p + 1), m]*Y[i1,])
      if(model == "GMAR") {
        sample[i1, j1] <- mu_mt + sqrt(sigmas[m])*rnorm(1)
      } else if(model == "StMAR") {
        sigma_mt <- sigmas[m]*(dfs[m] - 2 + matProd[m])/(dfs[m] - 2 + p)
        sample[i1, j1] <- mu_mt + sqrt(sigma_mt*(dfs[m] + p - 2)/(dfs[m] + p))*rt(1, df=dfs[m] + p)
      } else { # model == "G-StMAR"
        if(m <= M1) {
          sample[i1, j1] <- mu_mt + sqrt(sigmas[m])*rnorm(1)
        } else {
          sigma_mt <- sigmas[m]*(dfs[m - M1] - 2 + matProd[m])/(dfs[m - M1] - 2 + p)
          sample[i1, j1] <- mu_mt + sqrt(sigma_mt*(dfs[m - M1] + p - 2)/(dfs[m - M1] + p))*rt(1, df=dfs[m - M1] + p)
        }
      }

      if(p == 1) {
        Y[i1 + 1] <- sample[i1, j1]
      } else {
        Y[i1 + 1,] <- c(sample[i1, j1], Y[i1, 1:(p - 1)])
      }
    }
  }

  if(ntimes == 1 & drop == TRUE) {
    sample <- as.vector(sample)
    component <- as.vector(component)
    mixing_weights <- matrix(mixing_weights, ncol=M, nrow=nsimu, byrow=FALSE)
  }
  list(sample=sample,
       component=component,
       mixing_weights=mixing_weights)
}



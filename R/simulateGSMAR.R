#' @import stats
#'
#' @title Simulate values from GMAR, StMAR, and G-StMAR processes
#'
#' @description \code{simulateGSMAR} simulates values from the specified GMAR, StMAR, or G-StMAR process.
#'  Can be utilized for forecasting future values of the process.
#'
#' @param gsmar object of class \code{'gsmar'} created with the function \code{fitGSMAR} or \code{GSMAR}.
#' @param nsimu a positive integer specifying how many values (ahead from \code{init_values}) will be simulated.
#' @param init_values a numeric vector with length \code{>=p} specifying the initial values for the simulation. The \strong{last}
#'  element will be used as the initial value for the first lag, the second last element will be initial value for the second lag, etc.
#'  If not specified, initial values will be simulated from the process's stationary distribution.
#' @param ntimes a positive integer specifying how many sets of simulations should be performed.
#' @param drop if \code{TRUE} (default) then the components of the returned list are coerced to lower dimension if \code{ntimes==1},
#'   i.e., \code{$sample} and \code{$component} will be vectors and \code{$mixing_weights} will be matrix.
#' @details The argument \code{ntimes} is intended for forecasting: a GSMAR process can be forecasted by simulating its
#'  possible future values. One can perform a large number of sets of simulations and calculate the sample quantiles from
#'  the simulated values to obtain prediction intervals. See the forecasting example below for a hand-on demonstration.
#' @return If \code{drop==TRUE} and \code{ntimes==1} (default): \code{$sample} and \code{$component} are vectors
#'  and \code{$mixing_weights} is a (\code{nsimu}\eqn{xM}) matrix. Otherwise, returns a list with...
#'   \describe{
#'     \item{\code{$sample}}{a size (\code{nsimu}\eqn{x}\code{ntimes}) matrix containing the simulated values.}
#'     \item{\code{$component}}{a size (\code{nsimu}\eqn{x}\code{ntimes}) matrix containing the information from which
#'      mixture component each value was generated from.}
#'     \item{\code{$mixing_weights}}{a size (\code{nsimu}\eqn{xMx}\code{ntimes}) array containing the mixing weights corresponding to the
#'      sample: the dimension \code{[i, , ]} is the time index, the dimension \code{[, i, ]} indicates the regime, and the dimension
#'      \code{[, , i]} indicates the i:th set of simulations.}
#'   }
#' @seealso \code{\link{fitGSMAR}}, \code{\link{GSMAR}}, \code{\link{predict.gsmar}},
#'  \code{\link{add_data}}, \code{\link{cond_moments}}, \code{\link{mixing_weights}}
#' @inherit loglikelihood references
#' @examples
#'  \donttest{
#'  # GMAR model:
#'  params12 <- c(0.18, 0.93, 0.01, 0.86, 0.68, 0.02, 0.88)
#'  gmar12 <- GSMAR(p=1, M=2, params=params12, model="GMAR")
#'  sim12 <- simulateGSMAR(gmar12, nsimu=500)
#'  ts.plot(sim12$sample)
#'  ts.plot(sim12$component)
#'  ts.plot(sim12$mixing_weights, col=rainbow(2), lty=2)
#'
#'
#'  # G-StMAR model, with initial values:
#'  params12gs <- c(1.38, 0.88, 0.27, 3.8, 0.74, 3.15, 0.8, 3.6)
#'  gstmar12 <- GSMAR(p=1, M=c(1, 1), params=params12gs,
#'  model="G-StMAR")
#'  sim12gs <- simulateGSMAR(gstmar12, nsimu=500, init_values=5:6)
#'  ts.plot(sim12gs$sample)
#'  ts.plot(sim12gs$component)
#'  ts.plot(sim12gs$mixing_weights, col=rainbow(2), lty=2)
#'
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
#' }
#' @export

simulateGSMAR <- function(gsmar, nsimu, init_values, ntimes=1, drop=TRUE) {
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
  } else if(model == "StMAR") {
    M1 <- 0
  } else { # model == "GMAR"
    M1 <- M
  }

  if(!missing(init_values)) {
    if(length(init_values) < p) {
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
  params <- remove_all_constraints(p=p, M=M_orig, params=gsmar$params, model=model, restricted=gsmar$model$restricted,
                                 constraints=gsmar$model$constraints)
  if(gsmar$model$parametrization == "mean") { # For simplicity switch to use intercept parametrization in all cases
    params <- change_parametrization(p=p, M=M_orig, params=params, model=model, restricted=FALSE, constraints=NULL, change_to="intercept")
  }
  pars <- pick_pars(p=p, M=M_orig, params=params, model=model, restricted=FALSE, constraints=NULL)
  alphas <- pick_alphas(p=p, M=M_orig, params=params, model=model, restricted=FALSE, constraints=NULL)
  dfs <- pick_dfs(p=p, M=M_orig, params=params, model=model)
  sigmas <- pars[p + 2,] # sigma^2
  parameter_checks(p=p, M=M_orig, params=params, model=model, restricted=FALSE, constraints=NULL)

  # Create a container for the simulated values and initial values.
  # First row for initival values vector, and t+1:th row for (Y_t,Y_[t-1],...,Y_[t-p+1]), t=1,...,nsimu
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

  # Calculate expected values mu_m (KMS 2015, p.250, MPS eq.(4))
  mu <- pars[1, ]/(1 - colSums(pars[2:(p + 1), , drop=FALSE]))
  mu_mp <- tcrossprod(rep(1, p), mu) # \mu_m*1_p, column for each component

  # If initial values are missing simulate them from the processes stationary distribution
  if(missing(init_values)) {
    mv_samples <- numeric(p)
    m <- sample.int(M, size=1, prob=alphas) # From which mixture component the p-dimensional initial value is drawn from
    if(model == "GMAR" || (model == "G-StMAR" && m <= M1)) { # Draw from GMAR type component
      Gamma_m <- solve(invG[, , m])
      L <- t(chol(Gamma_m))
      mv_samples <- mu_mp[,m] + L%*%rnorm(p) # sample from N(mu_mp, Gamma_m)
    } else { # model == "StMAR" || (model == "G-StMAR" && m > M1); Draw from StMAR type component
      # Note that we use Student's t parametrization: mean, covariance (=v/(v - 2)*Scale), dfs (see MPS 2018, Supplementary material).
      Gamma_m <- solve(invG[, , m])
      L <- t(chol((dfs[m - M1] - 2)/dfs[m - M1]*Gamma_m)) # M1 = 0 for StMAR models
      Z <- L%*%rnorm(p) # Sample from N(0, ((v - 2)/v)/Gamma_m)
      U <- rchisq(1, df=dfs[m - M1]) # Sample from chisq(v)
      mv_samples <- mu_mp[,m] + Z*sqrt(dfs[m - M1]/U) # sample from St(mu_mp, Sigma_m, v)
    }
    Y[1,] <- mv_samples # Initial values sampled from the stationary distribution
  } else {
    init_values <- init_values[(length(init_values) - p + 1):length(init_values)]
    Y[1,] <- rev(init_values) # Take the last value to be in the first column, the second last in the second column, etc.
  }

  # Initialize data storages
  sample <- matrix(nrow=nsimu, ncol=ntimes)
  component <- matrix(nrow=nsimu, ncol=ntimes)
  mixing_weights <- array(dim=c(nsimu, M, ntimes))

  for(j1 in 1:ntimes) { # For each set of simulations...

    ### Start simulation ###
    for(i1 in 1:nsimu) {
      # Calculate log multinormal values (KMS 2015, eq.(7)) or log Student's t values (PMS 2018, p.5); for each regime 1,...,M.
      matProd <- vapply(1:M, function(i2) crossprod(Y[i1,] - mu_mp[,i2], as.matrix(invG[, , i2]))%*%(Y[i1,] - mu_mp[,i2]), numeric(1))
      if(model == "GMAR" || model == "G-StMAR") { # GMAR type regimes, M1 = M for GMAR models
        logmv_valuesM1 <- vapply(1:M1, function(i2) -0.5*p*log(2*base::pi) - 0.5*log(detG[i2]) - 0.5*matProd[i2], numeric(1))
      } else {
        logmv_valuesM1 <- numeric(0)
      }
      if(model == "StMAR" || model == "G-StMAR") { # StMAR type regimes, M1 = 0 for StMAR models
        logmv_valuesM2 <- vapply((M1 + 1):M, function(i2) lgamma(0.5*(p + dfs[i2 - M1])) - 0.5*p*log(pi) - 0.5*p*log(dfs[i2 - M1] - 2) -
                                 lgamma(0.5*dfs[i2 - M1]) - 0.5*log(detG[i2]) - 0.5*(p + dfs[i2 - M1])*log(1 + matProd[i2]/(dfs[i2 - M1] - 2)),
                                 numeric(1))
      } else {
        logmv_valuesM2 <- numeric(0)
      }
      logmv_values <- c(logmv_valuesM1, logmv_valuesM2)

      # Calculate the alpha_mt mixing weights (KMS 2015, eq.(8), PMS 2018, eq.(11)).
      alpha_mt <- get_alpha_mt(M=M, log_mvnvalues=logmv_values, alphas=alphas,
                               epsilon=epsilon, also_l_0=FALSE)

      # Draw the component and store the values
      m <- sample.int(M, size=1, prob=alpha_mt)
      component[i1, j1] <- m
      mixing_weights[i1, , j1] <- alpha_mt

      # Draw sample and store it
      mu_mt <- pars[1, m] + sum(pars[2:(p + 1), m]*Y[i1,])
      if(model == "GMAR" || (model == "G-StMAR" && m <= M1)) { # Draw from GMAR type regime
        sample[i1, j1] <- mu_mt + sqrt(sigmas[m])*rnorm(1) # MPS 2018, eq.(3)
      } else { # model == "StMAR" || (model == "G-StMAR" && m > M1)
        sigma_mt <- sigmas[m]*(dfs[m - M1] - 2 + matProd[m])/(dfs[m - M1] - 2 + p)
        sample[i1, j1] <- mu_mt + sqrt(sigma_mt*(dfs[m - M1] + p - 2)/(dfs[m - M1] + p))*rt(1, df=dfs[m - M1] + p) # MPS 2018, eq.(7)
      }

      # Setup for the next round
      if(p == 1) {
        Y[i1 + 1] <- sample[i1, j1]
      } else {
        Y[i1 + 1,] <- c(sample[i1, j1], Y[i1, 1:(p - 1)])
      }
    }
  }

  if(ntimes == 1 & drop) {
    sample <- as.vector(sample)
    component <- as.vector(component)
    mixing_weights <- matrix(mixing_weights, ncol=M, nrow=nsimu, byrow=FALSE)
  }
  list(sample=sample,
       component=component,
       mixing_weights=mixing_weights)
}



#' @import stats
#'
#' @title Simulate values from GMAR, StMAR or G-StMAR process
#'
#' @description \code{simulateGMAR} simulates values from the specified GMAR, StMAR or G-StMAR process.
#'
#' @inheritParams loglikelihood
#' @param nsimu a positive integer specifying how many values will be simulated.
#' @param initvalues an (optional) size \eqn{(px1)} vector (where \eqn{p} is the order of AR coefficients) specifying the initial
#'  values for the simulation. If not specified, initial values will be simulated from the process's stationary distribution.
#' @param ntimes a positive integer specifying how many sets of simulations should be performed. If larger than one then
#'   only the samples are returned. Uses the same initial values for each set. Default is one.
#' @return Returns a list with...
#'   \describe{
#'     \item{\code{$sample}}{a numeric vector containing the simulated values.}
#'     \item{\code{$component}}{a numeric vector the containing the information from which component each value is generated from.}
#'     \item{\code{$mixingWeights}}{a size (\code{nsimu}\eqn{xM}) matrix containing the mixing weights corresponding to the
#'      sample so that \eqn{i}:th column is for \eqn{i}:th component.}
#'   }
#'   Or if \strong{\code{ntimes>1}} returns a matrix containing the samples so that the \code{i}:th
#'   sample can be obtained at \code{i}:th column.
#' @inherit loglikelihood references
#' @examples
#' # GMAR process
#' params12 <- c(-0.3, 0.9, 0.5, 0.1, -0.2, 0.1, 0.6)
#' simu12 <- simulateGMAR(1, 2, params12, nsimu=100)
#'
#' # Restricted GMAR process: 10 realizations, samples only
#' params12r <- c(1.4, 1.8, 0.9, 0.3, 3.2, 0.8)
#' simu12r <- simulateGMAR(1, 2, params12r, restricted=TRUE, nsimu=100, ntimes=10)
#'
#' # StMAR process with initial values
#' params22t <- c(0.1, 0.7, 0.1, 0.5, -0.1, 0.5, -0.2, 0.3, 0.6, 3, 10)
#' simu22t <- simulateGMAR(2, 2, params22t, StMAR=TRUE, nsimu=100, initvalues=c(0.1, 0.2))
#'
#' # Restricted StMAR process
#' params13tr <- c(0.1, 0.2, 0.3, -0.9, 0.1, 0.2, 0.3, 0.45, 0.35, 3, 9, 27)
#' simu13tr <- simulateGMAR(1, 3, params13tr, StMAR=TRUE, restricted=TRUE, nsimu=100)
#'
#' # G-StMAR process
#' params13gs <- c(0, 0.9, 1, -1, 0.5, 0.8, 1, -0.5, 0.5, 0.3, 0.4, 5, 7)
#' simu13gs <- simulateGMAR(1, c(1, 2), params13gs, GStMAR=TRUE, nsimu=100)
#'
#' # Restricted G-StMAR process
#' params22gsr <- c(-1, 1, -0.7, 0.2, 2, 1, 0.5, 7)
#' simu22gsr <- simulateGMAR(2, c(1, 1), params22gsr, GStMAR=TRUE, restricted=TRUE, nsimu=100)
#'
#' # Restricted GMAR process with p=4, where the first three AR coefficients are restricted to be zero
#' R <- as.matrix(c(0, 0, 0, 1))
#' params42rc <- c(0.4, 0.5, 0.9, 0.5, 0.6, 0.6)
#' simu42rc <-simulateGMAR(4, 2, params42rc, restricted=TRUE, constraints=TRUE, R=R, nsimu=100)
#'
#' # Mixture version of Heterogenuous Autoregressive (HAR) process
#' paramsHAR2 <- c(1, 0.3, 0.2, 0.1, 1, 1.5, 0.3, 0.25, -0.1, 0.6, 0.55)
#' r1 = c(1, rep(0, 21)); r2 = c(rep(0.2, 5), rep(0, 17)); r3 = rep(1/22, 22)
#' R0 = cbind(r1, r2, r3)
#' simuHAR2 <- simulateGMAR(22, 2, paramsHAR2, constraints=TRUE, R=list(R0, R0), nsimu=100)
#' @export

simulateGMAR <- function(p, M, params, StMAR=FALSE, GStMAR=FALSE, restricted=FALSE, constraints=FALSE, R, nsimu, initvalues, ntimes=1) {

  checkLogicals(StMAR=StMAR, GStMAR=GStMAR)
  checkPM(p, M, GStMAR=GStMAR)
  M_orig = M
  if(GStMAR==TRUE) {
    M1 = M[1]
    M2 = M[2]
    M = sum(M)
  }
  if(length(params)!=nParams(p=p, M=M_orig, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted, constraints=constraints, R=R)) {
    stop("The parameter vector is wrong dimension")
  }
  if(!missing(initvalues)) {
    if(length(initvalues)!=p) {
      stop("The length of initial values vector has to be equal to the order of AR coefficients p")
    }
  }
  if(nsimu<1) {
    stop("The number of simulations has to be equal or larger than one")
  } else if(nsimu%%1!=0) {
    stop("Argument nsimu has to be positive integer")
  } else if(ntimes<1 | ntimes%%1!=0) {
    stop("The argument ntimes should be a positive integer")
  }
  epsilon = round(log(.Machine$double.xmin)+10)

  # Reform and collect parameters
  if(constraints==TRUE) {
    checkConstraintMat(p, M, R, restricted=restricted)
    params = reformConstrainedPars(p, M_orig, params, R=R, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted)
  }
  pars_tmp = reformParameters(p, M_orig, params, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted)
  params = pars_tmp$params
  pars = pars_tmp$pars
  alphas = pars_tmp$alphas
  sigmas = pars[p+2,]
  if(StMAR==TRUE | GStMAR==TRUE) {
    dfs = pars_tmp$dfs
  }
  parameterChecks(p, M_orig, params, pars, alphas, StMAR=StMAR, GStMAR=GStMAR, constraints=constraints)

  # Create container for the simulated values and initial values. First row row initival values vector, and t+1:th row for (Y_t,Y_[t-1],...,Y_[t-p+1])
  Y = matrix(nrow=nsimu+1, ncol=p)

  # Calculate inverses of covariance matrices Gamma_m and their determinants
  invG = array(dim=c(p, p, M))
  detG = numeric(M)
  if(p==1) { # Lutkepohl 2005, s.15-29
    for(i1 in 1:M) {
      A = pars[p+1, i1]
      Sigma = as.matrix(sigmas[i1])
      VecGamma = solve(1-kronecker(A,A), Sigma)
      invG[,,i1] = as.matrix(1/VecGamma)
      detG[i1] = VecGamma
     }
  } else { # Galbraith and Galbraith 1974
     for(i1 in 1:M) {
      ARcoefs = pars[2:(p+1), i1]
      U = diag(1, p, p)
      V = diag(ARcoefs[p], p, p)
      for(i2 in 1:(p-1)) {
        U[(i2+1):p, i2] <- -ARcoefs[1:(p-i2)]
        V[(i2+1):p, i2] <- rev(ARcoefs[i2:(p-1)])
      }
      invG[,,i1] = (crossprod(U, U) - crossprod(V, V))/sigmas[i1]
      detG[i1] = 1/det(invG[,,i1])
     }
  }

  # Calculate expected values mu_m (Kalliovirta 2015, s.250)
  mu = vapply(1:M, function(i1) pars[1, i1]/(1-sum(pars[2:(p+1), i1])), numeric(1))
  mu_mp = matrix(vapply(1:M, function(i1) mu[i1]*matrix(rep(1,p)), numeric(p)), ncol=M) # Column for each component

  # If initial values are missing simulate them from the processes stationary distribution
  if(missing(initvalues)) {
    mv_samples = matrix(ncol=M, nrow=p)
    if(StMAR==FALSE & GStMAR==FALSE) {
      for(i1 in 1:M) {
        Gamma_m = solve(invG[,,i1])
        L = t(chol(Gamma_m))
        mv_samples[,i1] = mu_mp[,i1] + L%*%rnorm(p)
      }
    } else if(StMAR==TRUE) {
       for(i1 in 1:M) {
        Gamma_m = solve(invG[,,i1])
        L = t(chol((dfs[i1]-2)/dfs[i1]*Gamma_m))
        Z = L%*%rnorm(p) # Sample from N(0,K)
        U = rchisq(p, df=dfs[i1]) # Sample from chisq(v)
        mv_samples[,i1] = mu_mp[,i1] + Z*sqrt(dfs[i1]/U) # sample from student's t
      }
    } else { # If GStMAR==TRUE
      for(i1 in 1:M) {
        if(i1<=M1) { # p-dimensional Gaussian samples
          Gamma_m = solve(invG[,,i1])
          L = t(chol(Gamma_m))
          mv_samples[,i1] = mu_mp[,i1] + L%*%rnorm(p)
        } else { # p-dimensional Student samples
          Gamma_m = solve(invG[,,i1])
          L = t(chol((dfs[i1-M1]-2)/dfs[i1-M1]*Gamma_m))
          Z = L%*%rnorm(p) # Sample from N(0,K)
          U = rchisq(p, df=dfs[i1-M1]) # Sample from chisq(v)
          mv_samples[,i1] = mu_mp[,i1] + Z*sqrt(dfs[i1-M1]/U) # sample from student's t
        }
      }
    }
    if(p==1) {
      Y[1,] = sum(sapply(1:M, function(i1) alphas[i1]*mv_samples[i1])) # Initial values sampled from the stationary distribution
    } else {
      Y[1,] = rowSums(vapply(1:M, function(i1) alphas[i1]*mv_samples[,i1], numeric(p))) # Initial values sampled from the stationary distribution
    }
  } else {
    Y[1,] = rev(initvalues)
  }

  ret0 = matrix(nrow=nsimu, ncol=ntimes)
  for(j1 in 1:ntimes) {
    sample = numeric(nsimu)
    sampleWeights = matrix(nrow=nsimu, ncol=M)
    whichComp = numeric(nsimu)

    ### Start simulation ###
    for(i1 in 1:nsimu) {

      # Calculate log multinormal values (Kalliovirta 2015, s250 eq.(7)) or log student-t values
      matProd = vapply(1:M, function(i2) crossprod(Y[i1,]-mu_mp[,i2], as.matrix(invG[,,i2]))%*%(Y[i1,]-mu_mp[,i2]), numeric(1))
      if(StMAR==FALSE & GStMAR==FALSE) {
        logmv_values = vapply(1:M, function(i2) -0.5*p*log(2*pi)-0.5*log(detG[i2])-0.5*matProd[i2], numeric(1))
      } else if(StMAR==TRUE) {
        logmv_values = vapply(1:M, function(i2) lgamma(0.5*(p+dfs[i2]))-0.5*p*log(pi)-0.5*p*log(dfs[i2]-2)-lgamma(0.5*dfs[i2]) -
                                0.5*log(detG[i2]) - 0.5*(p+dfs[i2])*log(1 + matProd[i2]/(dfs[i2]-2)), numeric(1))
      } else { # If GStMAR==TRUE
        logmv_values = c(vapply(1:M1, function(i2) -0.5*p*log(2*pi)-0.5*log(detG[i2])-0.5*matProd[i2], numeric(1)),
                         vapply((M1+1):M, function(i2) lgamma(0.5*(p+dfs[i2-M1]))-0.5*p*log(pi)-0.5*p*log(dfs[i2-M1]-2)-lgamma(0.5*dfs[i2-M1]) -
                                  0.5*log(detG[i2]) - 0.5*(p+dfs[i2-M1])*log(1 + matProd[i2]/(dfs[i2-M1]-2)), numeric(1)))
      }

      # Calculate the alpha_mt mixing weights (Kalliovirta 2015, s.250 eq.(8)). Close to zero values handled with Brobdingnag
      if(M==1) {
        alpha_mt = 1
      } else if(any(logmv_values < epsilon)) {
        numerators = lapply(1:M, function(i2) alphas[i2]*Brobdingnag::as.brob(exp(1))^logmv_values[i2])
        denominator = Reduce("+", numerators)
        alpha_mt = vapply(1:M, function(i2) as.numeric(numerators[[i2]]/denominator), numeric(1))
      } else {
        mv_values = exp(logmv_values)
        denominator = sum(alphas*mv_values)
        alpha_mt = alphas*mv_values/denominator
      }
      sampleWeights[i1,] = alpha_mt

      # Draw the component
      m = sample.int(M, size=1, prob=alpha_mt)
      whichComp[i1] = m

      # Draw sample and store it
      mu_mt = pars[1, m] + sum(pars[2:(p+1), m]*Y[i1,])
      if(StMAR==FALSE & GStMAR==FALSE) { # If GMAR
        sample[i1] = mu_mt + sqrt(sigmas[m])*rnorm(1)
      } else if (StMAR==TRUE) {
        sigma_mt = sigmas[m]*(dfs[m] - 2 + matProd[m])/(dfs[m] - 2 + p)
        sample[i1] = mu_mt + sqrt(sigma_mt*(dfs[m]+p-2)/(dfs[m]+p))*rt(1, df=dfs[m]+p)
      } else { # If GStMAR==TRUE
        if(m<=M1) {
          sample[i1] = mu_mt + sqrt(sigmas[m])*rnorm(1)
        } else {
          sigma_mt = sigmas[m]*(dfs[m-M1] - 2 + matProd[m])/(dfs[m-M1] - 2 + p)
          sample[i1] = mu_mt + sqrt(sigma_mt*(dfs[m-M1]+p-2)/(dfs[m-M1]+p))*rt(1, df=dfs[m-M1]+p)
        }
      }

      if(p==1) {
        Y[i1+1] = sample[i1]
      } else {
        Y[i1+1,] = c(sample[i1], Y[i1,1:(p-1)])
      }
    }
    if(ntimes==1) {
      ret = list(sample, whichComp, sampleWeights)
      names(ret) = c("sample", "component", "mixingWeights")
      return(ret)
    } else {
      ret0[,j1] = sample
    }
  }
  return(ret0)
}



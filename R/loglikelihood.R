#' @import stats
#'
#' @title Compute the log-likelihood of Gaussian or/and Student's t Mixture Autoregressive model
#'
#' @description FOR INTERNAL USE. \code{loglikelihood_int} computes the log-likelihood value of the specified GMAR, StMAR or G-StMAR model for the given data.
#'
#' @param data a numeric vector or column matrix containing the data. \code{NA} values are not supported.
#' @param p a positive integer specifying the order of AR coefficients.
#' @param M a positive integer specifying the number of mixture components. Except for G-StMAR model a size (2x1) vector specifying the number of GMAR-components M1 in
#'  the first element and StMAR-components M2 in the second - then the total number of components is M=M1+M2.
#' @param params a real valued parameter vector specifying the model.
#'  \describe{
#'    \item{For \strong{non-restricted} models:}{
#'      \describe{
#'        \item{For \strong{GMAR} model:}{Size \eqn{(M(p+3)-1x1)} vector \strong{\eqn{\theta}}\eqn{=}(\strong{\eqn{\upsilon_{1}}},...,\strong{\eqn{\upsilon_{M}}},
#'          \eqn{\alpha_{1},...,\alpha_{M-1}}), where \strong{\eqn{\upsilon_{m}}}\eqn{=(\phi_{m,0},}\strong{\eqn{\phi_{m}}}\eqn{,
#'          \sigma_{m}^2)} and \strong{\eqn{\phi_{m}}}=\eqn{(\phi_{m,1},...,\phi_{m,p}), m=1,...,M}.}
#'        \item{For \strong{StMAR} model:}{Size \eqn{(M(p+4)-1x1)} vector (\strong{\eqn{\theta, \nu}})\eqn{=}(\strong{\eqn{\upsilon_{1}}},...,\strong{\eqn{\upsilon_{M}}},
#'          \eqn{\alpha_{1},...,\alpha_{M-1}, \nu_{1},...,\nu_{M}}).}
#'        \item{For \strong{G-StMAR} model:}{Size \eqn{(M(p+3)+M2-1x1)} vector (\strong{\eqn{\theta, \nu}})\eqn{=}(\strong{\eqn{\upsilon_{1}}},...,\strong{\eqn{\upsilon_{M}}},
#'          \eqn{\alpha_{1},...,\alpha_{M-1}, \nu_{M1+1},...,\nu_{M}}).}
#'        \item{With \strong{linear constraints}:}{Replace the vectors \strong{\eqn{\phi_{m}}} with vectors \strong{\eqn{\psi_{m}}} and provide a  list of constraint
#'          matrices \strong{R} that satisfy \strong{\eqn{\phi_{m}}}\eqn{=}\strong{\eqn{R_{m}\psi_{m}}} for all \eqn{m=1,...,M}, where
#'          \strong{\eqn{\psi_{m}}}\eqn{=(\psi_{m,1},...,\psi_{m,q_{m}})}.}
#'      }
#'    }
#'    \item{For \strong{restricted} models:}{
#'      \describe{
#'        \item{For \strong{GMAR} model:}{Size \eqn{(3M+p-1x1)} vector \strong{\eqn{\theta}}\eqn{=(\phi_{1,0},...,\phi_{M,0},}\strong{\eqn{\phi}}\eqn{,
#'          \sigma_{1}^2,...,\sigma_{M}^2,\alpha_{1},...,\alpha_{M-1})}, where \strong{\eqn{\phi}}=\eqn{(\phi_{1},...,\phi_{M})}.}
#'        \item{For \strong{StMAR} model:}{Size \eqn{(4M+p-1x1)} vector (\strong{\eqn{\theta, \nu}})\eqn{=(\phi_{1,0},...,\phi_{M,0},}\strong{\eqn{\phi}}\eqn{,
#'          \sigma_{1}^2,...,\sigma_{M}^2,\alpha_{1},...,\alpha_{M-1}, \nu_{1},...,\nu_{M})}.}
#'        \item{For \strong{G-StMAR} model:}{Size \eqn{(3M+M2+p-1x1)} vector (\strong{\eqn{\theta, \nu}})\eqn{=(\phi_{1,0},...,\phi_{M,0},}\strong{\eqn{\phi}}\eqn{,
#'          \sigma_{1}^2,...,\sigma_{M}^2,\alpha_{1},...,\alpha_{M-1}, \nu_{M1+1},...,\nu_{M})}.}
#'        \item{With \strong{linear constraints}:}{Replace the vector \strong{\eqn{\phi}} with vector \strong{\eqn{\psi}} and provide a constraint matrix
#'          \strong{\eqn{R}} that satisfies \strong{\eqn{\phi}}\eqn{=}\strong{\eqn{R\psi}}, where
#'          \strong{\eqn{\psi}}\eqn{=(\psi_{1},...,\psi_{q})}.}
#'      }
#'    }
#'  }
#'  Symbol \eqn{\phi} denotes an AR coefficient, \eqn{\sigma^2} a variance, \eqn{\alpha} a mixing weight and \eqn{\nu} a degrees of
#'  freedom parameter. In the \strong{G-StMAR} model the first M1 components are GMAR-type and the rest M2 components are StMAR-type.
#'  Note that in the case \strong{M=1} the parameter \eqn{\alpha} is dropped, and in the case of \strong{StMAR} or \strong{G-StMAR} model
#'  the degrees of freedom parameters \eqn{\nu_{m}} have to be larger than \eqn{2}.
#' @param restricted an (optional) logical argument stating whether the AR coefficients \eqn{\phi_{m,1},...,\phi_{m,p}} are restricted
#'  to be the same for all regimes. Default is \code{FALSE}.
#' @param StMAR an (optional) logical argument stating whether StMAR model should be considered instead of GMAR model. Default is \code{FALSE}.
#' @param GStMAR an (optional) logical argument stating whether G-StMAR model should be considered instead of GMAR model. In G-StMAR model the first M1 components
#'  are GMAR-type and the rest M2 components are StMAR-type. Default is \code{FALSE}.
#' @param constraints an (optional) logical argument stating whether general linear constraints should be applied to the model. Default is \code{FALSE}.
#' @param R Specifies the linear constraints.
#'   \describe{
#'   \item{For \strong{non-restricted} models:}{a list of size \eqn{(pxq_{m})} constraint matrices \strong{\eqn{R_{m}}} of full column rank
#'     satisfying \strong{\eqn{\phi_{m}}}\eqn{=}\strong{\eqn{R_{m}\psi_{m}}} for all \eqn{m=1,...,M}, where
#'     \strong{\eqn{\phi_{m}}}\eqn{=(\phi_{m,1},...,\phi_{m,p})} and \strong{\eqn{\psi_{m}}}\eqn{=(\psi_{m,1},...,\psi_{m,q_{m}})}.}
#'   \item{For \strong{restricted} models:}{a size \eqn{(pxq)} constraint matrix \strong{\eqn{R}} of full column rank satisfying
#'     \strong{\eqn{\phi}}\eqn{=}\strong{\eqn{R\psi}}, where \strong{\eqn{\phi}}\eqn{=(\phi_{1},...,\phi_{p})} and
#'     \strong{\eqn{\psi}}\eqn{=\psi_{1},...,\psi_{q}}.}
#'   }
#'   Symbol \eqn{\phi} denotes an AR coefficient. Note that regardless of any constraints, the nominal order of AR coefficients is alway \code{p} for all regimes.
#'   This argument is ignored if \code{constraints==FALSE}.
#' @param conditional an (optional) logical argument specifying wether the conditional or exact log-likehood function should be used. Default is \code{TRUE}.
#' @param boundaries an (optional) logical argument. If \code{TRUE} then \code{loglikelihood} returns \code{minval} if...
#' \itemize{
#'   \item any component variance is not larger than zero,
#'   \item any parametrized mixing weight \eqn{\alpha_{1},...,\alpha_{M-1}} is not larger than zero,
#'   \item sum of the parametrized mixing weights is not smaller than one,
#'   \item if the model is not stationary,
#'   \item or if \code{StMAR=TRUE} or \code{GStMAR=TRUE} and any degrees of freedom parameter \eqn{\nu_{m}} is not larger than two or is larger than 342-p.
#' }
#' Default is \code{FALSE}.
#' @param checks an (optional) logical argument defining whether argument checks are made. If \code{FALSE} then no argument checks
#' such as stationary checks etc are made. The default is \code{TRUE}.
#' @param returnTerms set \code{TRUE} if the terms \eqn{l_{t}: t=1,..,T} in the log-likelihood function (see \emph{KMS 2015, eq.(13)})
#'  should not be summed to calculate the log-likelihood value, but to be returned as a numeric vector. Default is \code{FALSE}.
#' @param epsilon an (optional) negative real number specifying the logarithm of the smallest positive non-zero number that will be
#'  handled without external packages. Too small value may lead to a failure or biased results and too large value will make the code
#'  run significantly slower. Default is \code{round(log(.Machine$double.xmin)+10)} and should not be adjusted too much.
#' @param minval a negative real number defining the log-likelihood value that will be returned with \code{boundaries==TRUE}
#'  when the parameter is out of the boundaries. Ignored if \code{boundaries==FALSE}.
#' @return
#'  \describe{
#'   \item{By default:}{log-likelihood value of the specified GMAR, StMAR or G-StMAR model,}
#'   \item{If \code{returnTerms==TRUE:}}{size \eqn{Tx1} numeric vector containing the terms \eqn{l_{t}}.}
#'  }
#' @references
#'  \itemize{
#'    \item Kalliovirta L., Meitz M. and Saikkonen P. (2015) Gaussian Mixture Autoregressive model for univariate time series.
#'            \emph{Journal of Time Series Analysis}, \strong{36}, 247-266.
#'    \item Lutkepohl H. New Introduction to Multiple Time Series Analysis,
#'            \emph{Springer}, 2005.
#'    \item Galbraith, R., Galbraith, J., (1974). On the inverses of some patterned matrices arising
#'            in the theory of stationary time series. \emph{Journal of Applied Probability} \strong{11}, 63-71.
#'    \item References regarding StMAR and G-StMAR models will be updated after they are published.
#'  }

loglikelihood_int <- function(data, p, M, params, StMAR=FALSE, GStMAR=FALSE, restricted=FALSE, constraints=FALSE, R, conditional=TRUE, boundaries=FALSE, checks=TRUE, returnTerms=FALSE, epsilon, minval) {
  M_orig = M
  if(GStMAR==TRUE) {
    M1 = M[1]
    M2 = M[2]
    M = sum(M)
  }
  if(missing(epsilon)) {
    epsilon = round(log(.Machine$double.xmin)+10)
  }

  # Reform and collect parameters
  if(constraints==TRUE) {
    if(checks==TRUE) {
      checkConstraintMat(p=p, M=M_orig, R=R, restricted=restricted)
    }
    params = reformConstrainedPars(p=p, M=M_orig, params=params, R=R, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted)
  }
  if(M==1) {
    pars = matrix(params[1:(M*(p+2))], ncol=M)
    alphas = c(1)
    if(StMAR==TRUE) {
      dfs = params[p+3]
    }
  } else {
    if(restricted==FALSE) {
      pars = matrix(params[1:(M*(p+2))], ncol=M) # Component parameters by column (except alphas and dfs)
      alphas = params[(M*(p+2)+1):(M*(p+3)-1)]
      alphas = c(alphas, 1-sum(alphas))
      if(StMAR==TRUE) {
        dfs = params[(M*(p+3)):(M*(p+4)-1)] # degrees of freedom
      } else if(GStMAR==TRUE) {
        dfs = params[(M*(p+3)):(M*(p+3)+M2-1)]
      }
    } else {
      # If restricted TRUE: transform the restricted parameter vector into the standard form. Then do everything like for non-restricted.
      phi0 = params[1:M]
      arcoefs = matrix(rep(params[(M+1):(M+p)], M), ncol=M)
      variances = params[(M+p+1):(p+2*M)]
      pars = rbind(phi0, arcoefs, variances)
      alphas = params[(p+2*M+1):(3*M+p-1)]
      if(StMAR==TRUE) {
        dfs = params[(3*M+p):(4*M+p-1)] # degrees of freedom
        params = c(as.vector(pars), alphas, dfs)
      } else if(GStMAR==TRUE) {
        dfs = params[(3*M+p):(3*M+p-1+M2)]
        params = c(as.vector(pars), alphas, dfs)
      } else {
        params = c(as.vector(pars), alphas)
      }
      alphas = c(alphas, 1-sum(alphas))
    }
  }
  rownames(pars) = NULL
  sigmas = pars[p+2,]

  # Return minval if parameters are out of their boundaries.
  if(boundaries==TRUE) {
    if(any(pars[p+2,]<=0)) {
      return(minval)
    } else if(M>=2 & sum(alphas[-M])>=1) {
      return(minval)
    } else if(any(alphas<=0)) {
      return(minval)
    } else if(!isStationary_int(p, M, params)) {
      return(minval)
    }
    if(StMAR==TRUE | GStMAR==TRUE) {
      if(any(dfs<=2)) {
        return(minval)
      }
    }
  }

  if(checks==TRUE) {
    data = checkAndCorrectData(data=data, p=p)
    parameterChecks(p=p, M=M_orig, params=params, pars=pars, alphas=alphas, StMAR=StMAR, GStMAR=GStMAR, constraints=constraints)
  }
  n_obs = length(data)

  #### Start evaluating the log-likelihood ####

  # Expected values mu_m (Kalliovirta 2015, s.250)
  mu = vapply(1:M, function(i1) pars[1, i1]/(1-sum(pars[2:(p+1), i1])), numeric(1))

  # Observed data: y_(-p+1),...,y_0,y_1,...,y_(n_obs-p). First row denotes vector y_0, i:th row vector y_[i-1] and last row denotes the vector y_T.
  Y = vapply(1:p, function(i1) data[(p-i1+1):(n_obs-i1+1)], numeric(n_obs-p+1) )

  # Calculate inverse Gamma_m and calculate the matrix products in mv normal and t-distribution (Galbraith and Galbraith 1974)
  matProd = matrix(nrow=n_obs-p+1, ncol=M)
  invG = array(dim=c(p, p, M))
  if(p==1) { # Form vectorized Gamma_m (Lutkepohl (2005), s.15-29)
    for(i1 in 1:M) {
      A = pars[p+1, i1]
      Sigma = as.matrix(sigmas[i1])
      VecGamma = solve(1-kronecker(A,A), Sigma)
      invG[,,i1] = as.matrix(1/VecGamma)
      matProd[,i1] = (Y-mu[i1])*invG[,,i1]*(Y-mu[i1])
    }
  } else { # Inverse formula by Galbraith, R., Galbraith, J., (1974)
    for(i1 in 1:M) {
      ARcoefs = pars[2:(p+1), i1]
      U = diag(1, p, p)
      V = diag(ARcoefs[p], p, p)
      for(i2 in 1:(p-1)) {
        U[(i2+1):p, i2] <- -ARcoefs[1:(p-i2)]
        V[(i2+1):p, i2] <- rev(ARcoefs[i2:(p-1)])
      }
      invG[,,i1] = (crossprod(U, U) - crossprod(V, V))/sigmas[i1]
      matProd[,i1] = rowSums((Y-mu[i1]*rep(1,p))%*%invG[,,i1]*(Y-mu[i1]*rep(1,p)))
    }
  }

  # Calculate the log multivariate normal or student's t values (Kalliovirta 2015, s250 eq.(7) or MPS 2018) for each vector y_t and for each m=1,..,M
  # First row for initial values y_0 (as denoted by Kalliovirta 2015) and i:th row for y_(i-1). First column for component m=1 and j:th column for m=j.
  logmv_values = matrix(nrow=(n_obs-p+1), ncol=M)
  if(StMAR==FALSE & GStMAR==FALSE) {
    for(i1 in 1:M) {
      detG = 1/det(as.matrix(invG[,,i1]))
      logmv_values[,i1] = -0.5*p*log(2*pi)-0.5*log(detG)-0.5*matProd[,i1]
    }
  } else if(StMAR==TRUE){
    for(i1 in 1:M) {
      detG = 1/det(as.matrix(invG[,,i1]))
      logC = lgamma(0.5*(p+dfs[i1]))-0.5*p*log(pi)-0.5*p*log(dfs[i1]-2)-lgamma(0.5*dfs[i1])
      logmv_values[,i1] = logC - 0.5*log(detG) - 0.5*(p+dfs[i1])*log(1 + matProd[,i1]/(dfs[i1]-2))
    }
  } else {  # If GStMAR==TRUE
    for(i1 in 1:M) {
      detG = 1/det(as.matrix(invG[,,i1]))
      if(i1 <= M1) { # Multinormals
        logmv_values[,i1] = -0.5*p*log(2*pi)-0.5*log(detG)-0.5*matProd[,i1]
      } else { # Multistudents
        logC = lgamma(0.5*(p+dfs[i1-M1]))-0.5*p*log(pi)-0.5*p*log(dfs[i1-M1]-2)-lgamma(0.5*dfs[i1-M1])
        logmv_values[,i1] = logC - 0.5*log(detG) - 0.5*(p+dfs[i1-M1])*log(1 + matProd[,i1]/(dfs[i1-M1]-2))
      }
    }
  }

  # Calculate the alpha_mt mixing weights (Kalliovirta 2015, s.250 eq.(8))
  # First row for t=1, second for t=2, and i:th for t=i. First column for m=1, second for m=2 and j:th column for m=j.
  logmv_values0 = logmv_values[1:(n_obs-p),] # The last row is not needed because alpha_mt uses vector Y_(t-1)
  if(!is.matrix(logmv_values0)) logmv_values0 = as.matrix(logmv_values0)

  l_0 = 0 # "The first term" of the exact log-likelihood
  if(M==1) {
    alpha_mt = as.matrix(rep(1, n_obs-p))
    if(conditional==FALSE) { # Calculate "the first term" of the log-likelihood (Kalliovirta ym 2015, s.254 eq.(12))
      l_0 = logmv_values[1]
    }
  } else if(any(logmv_values0 < epsilon)) { # Close to zero values handled with Brobdingnag if needed
    numerators = lapply(1:M, function(i1) alphas[i1]*Brobdingnag::as.brob(exp(1))^logmv_values0[,i1]) # alphas[i1]*exp( as.brob(mvn_values[,i1]) )
    denominator = Reduce("+", numerators) # For all t=0,...,T
    alpha_mt = vapply(1:M, function(i1) as.numeric(numerators[[i1]]/denominator), numeric(n_obs-p))

    if(conditional==FALSE) {
      l_0 = log(Reduce("+", lapply(1:M, function(i1) numerators[[i1]][1])))
    }
  } else {
    mv_values0 = exp(logmv_values0)
    denominator = colSums(alphas*t(mv_values0))
    alpha_mt = t(alphas*t(mv_values0/denominator))

    if(conditional==FALSE) {
      l_0 = log(sum(alphas*mv_values0[1,]))
    }
  }

  # Calculate the conditional means mu_mt (Kalliovirta 2015, s.249 eq.(2)). First row for t=1, second for t=2 etc. First column for m=1, second column for m=2 etc.
  if(p==1) {
    mu_mt = vapply(1:M, function(i1) rep(pars[1,i1], nrow(Y)-1) + Y[1:(nrow(Y)-1),]*pars[2,i1], numeric(n_obs-p) )
  } else {
    mu_mt = vapply(1:M, function(i1) rep(pars[1,i1], nrow(Y)-1) + colSums(pars[2:(p+1),i1]*t(Y[1:(nrow(Y)-1),])), numeric(n_obs-p) )
  }

  # Calculate "the second term" of the log-likelihood (Kalliovirta 2015, s.254 eq.(12)-(13) )
  Y2 = Y[2:nrow(Y),1] # Only first column and rows 2...T are needed
  if(StMAR==FALSE & GStMAR==FALSE) {
    invsqrt_sigmas = sigmas^(-1/2)
    if(M==1) {
      lt_tmp = invsqrt_sigmas*dnorm((Y2-mu_mt)*invsqrt_sigmas)
    } else {
      lt_tmp = alpha_mt*dnorm((Y2-mu_mt)%*%diag(invsqrt_sigmas))%*%diag(invsqrt_sigmas)
    }
  } else if(StMAR==TRUE) { # If StMAR=TRUE
    matProd0 = matProd[1:(n_obs-p),] # Last row is not needed because sigma_t uses y_{t-1}
    if(M==1) {
      sigma_mt = sigmas*(dfs - 2 + matProd0)/(dfs - 2 + p) # Conditional variances
      lt_tmp = ((exp(lgamma(0.5*(1+dfs+p))-lgamma(0.5*(dfs+p)))/sqrt(pi*(dfs+p-2)))/sqrt(sigma_mt))*(1 + ((Y2-mu_mt)^2)/((dfs+p-2)*sigma_mt))^(-0.5*(1+dfs+p))
    } else {
      sigma_mt = t(dfs - 2 + t(matProd0))%*%diag(1/(dfs - 2 + p))%*%diag(sigmas)
      lt_tmp = alpha_mt*t(exp(lgamma(0.5*(1+dfs+p))-lgamma(0.5*(dfs+p)))/sqrt(pi*(dfs+p-2))/t(sqrt(sigma_mt)))*t(t(1 + ((Y2-mu_mt)^2)/(sigma_mt%*%diag(dfs+p-2)))^(-0.5*(1+dfs+p)))
    }
  } else {  # If GStMAR=TRUE
    # GMAR-components
    invsqrt_sigmasM1 = sigmas[1:M1]^(-1/2)
    if(M1==1) {
      lt_tmpM1 = alpha_mt[,1]*invsqrt_sigmasM1*dnorm((Y2-mu_mt[,1])*invsqrt_sigmasM1)
    } else { # TESTAA!!
      lt_tmpM1 = alpha_mt[,1:M1]*dnorm((Y2-mu_mt[,1:M1])%*%diag(invsqrt_sigmasM1))%*%diag(invsqrt_sigmasM1)
    }
    # StMAR-components
    sigmasM2 = sigmas[(M1+1):M]
    matProd0 = matProd[1:(n_obs-p),(M1+1):M]
    if(M2==1) {
      sigma_mt = sigmasM2*(dfs - 2 + matProd0)/(dfs - 2 + p) # Conditional variances
      lt_tmpM2 = alpha_mt[,(M1+1):M]*((exp(lgamma(0.5*(1+dfs+p))-lgamma(0.5*(dfs+p)))/sqrt(pi*(dfs+p-2)))/sqrt(sigma_mt))*(1 + ((Y2-mu_mt[,(M1+1):M])^2)/((dfs+p-2)*sigma_mt))^(-0.5*(1+dfs+p))
    } else { # TESTAA!!
      sigma_mt = t(dfs - 2 + t(matProd0))%*%diag(1/(dfs - 2 + p))%*%diag(sigmasM2)
      lt_tmpM2 = alpha_mt[,(M1+1):M]*t(exp(lgamma(0.5*(1+dfs+p))-lgamma(0.5*(dfs+p)))/sqrt(pi*(dfs+p-2))/t(sqrt(sigma_mt)))*t(t(1 + ((Y2-mu_mt[,(M1+1):M])^2)/(sigma_mt%*%diag(dfs+p-2)))^(-0.5*(1+dfs+p)))
    }
    lt_tmp = cbind(lt_tmpM1, lt_tmpM2)
  }
  l_t = rowSums(lt_tmp)

  if(returnTerms==TRUE) {
    return(log(l_t))
  } else {
    return(l_0 + sum(log(l_t)))
  }
}


#' @import stats
#'
#' @title Compute the log-likelihood of Gaussian or/and Student's t Mixture Autoregressive model
#'
#' @description \code{loglikelihood} computes the log-likelihood value of the specified GMAR, StMAR or G-StMAR model for the given data.
#'
#' @inheritParams loglikelihood_int
#' @inherit loglikelihood_int return references
#' @examples
#' # GMAR model
#' params13 <- c(1.4, 0.88, 0.26, 2.46, 0.82, 0.74, 5.0, 0.68, 5.2, 0.72, 0.2)
#' loglikelihood(VIX, 1, 3, params13)
#'
#' # Restricted GMAR model
#' params12r <- c(1.4, 1.8, 0.88, 0.29, 3.18, 0.84)
#' loglikelihood(VIX, 1, 2, params12r, restricted=TRUE)
#'
#' # StMAR model
#' params12t <- c(1.9, 0.85, 1.16, 1.22, 0.89, 0.13, 0.64, 3.1, 7.0)
#' loglikelihood(VIX, 1, 2, params12t, StMAR=TRUE)
#'
#' # Non-mixture version of StMAR model
#' params11t <- c(0.76, 0.93, 1.35, 2.4)
#' loglikelihood(VIX, 1, 1, params11t, StMAR=TRUE)
#'
#' # G-StMAR model
#' params12gs <- c(1.2, 0.8, 0.6, 1.3, 0.6, 1.1, 0.6, 3)
#' loglikelihood(VIX, 1, c(1,1), params12gs, GStMAR=TRUE)
#'
#' # Restricted G-StMAR model
#' params13gsr <- c(1.3, 2.2, 1.4, 0.8, 2.4, 4.6, 0.4, 0.25, 0.15, 20)
#' loglikelihood(VIX, 1, c(2,1), params13gsr, GStMAR=TRUE, restricted=TRUE)
#'
#' # GMAR model as a mixture of AR(2) and AR(1) models
#' R <- list(diag(1, ncol=2, nrow=2), as.matrix(c(1, 0)))
#' params22c <- c(1.2, 0.85, 0.04, 0.3, 3.3, 0.77, 2.8, 0.77)
#' loglikelihood(VIX, 2, 2, params22c, constraints=TRUE, R=R)
#'
#' # Such StMAR(3,2) that the AR coefficients are restricted to be
#' # the same for both regimes and that the second AR coefficients are
#' # constrained to zero.
#' params32trc <- c(2.2, 1.8, 0.88, -0.03, 2.4, 0.27, 0.40, 3.9, 1000)
#' loglikelihood(VIX, 3, 2, params32trc, StMAR=TRUE, restricted=TRUE,
#'               constraints=TRUE, R=matrix(c(1, 0, 0, 0, 0, 1), ncol=2))
#' @export

loglikelihood <- function(data, p, M, params, StMAR=FALSE, GStMAR=FALSE, restricted=FALSE, constraints=FALSE, R, conditional=TRUE, returnTerms=FALSE) {
  checkLogicals(StMAR=StMAR, GStMAR=GStMAR)
  checkPM(p=p, M=M, GStMAR=GStMAR)
  if(length(params)!=nParams(p=p, M=M, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted, constraints=constraints, R=R)) {
    stop("The parameter vector has wrong dimension")
  }
  return(loglikelihood_int(data=data, p=p, M=M, params=params, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted, constraints=constraints, R=R, conditional=conditional, checks=TRUE, returnTerms=returnTerms))
}


#' @import stats
#'
#' @title Calculate mixing weights of GMAR, StMAR or G-StMAR model
#'
#' @description FOR INTERNAL USE. \code{mixingWeights_int} calculates the mixing weights of the specified GMAR, StMAR or G-StMAR model and returns them as a matrix.
#'
#' @inheritParams loglikelihood_int
#' @return Returns size \eqn{(TxM)} matrix containing the mixing weights of the specified GMAR, StMAR or G-StMAR model so that
#'  \eqn{i}:th column corresponds to \eqn{i}:th mixing component or regime.
#' @inherit loglikelihood_int references

mixingWeights_int <- function(data, p, M, params, StMAR=FALSE, GStMAR=FALSE, restricted=FALSE, constraints=FALSE, R, checks=TRUE, epsilon) {
  M_orig = M
  if(GStMAR==TRUE) {
    M1 = M[1]
    M2 = M[2]
    M = sum(M)
  }
  if(missing(epsilon)) {
    epsilon = round(log(.Machine$double.xmin)+10)
  }

  # Reform and collect parameters
  if(constraints==TRUE) {
    if(checks==TRUE) {
      checkConstraintMat(p=p, M=M_orig, R=R, restricted=restricted)
    }
    params = reformConstrainedPars(p=p, M=M_orig, params=params, R=R, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted)
  }
  if(M==1) {
    pars = matrix(params[1:(M*(p+2))], ncol=M)
    alphas = c(1)
    if(StMAR==TRUE) {
      dfs = params[p+3]
    }
  } else {
    if(restricted==FALSE) {
      pars = matrix(params[1:(M*(p+2))], ncol=M) # Component parameters by column (except alphas and dfs)
      alphas = params[(M*(p+2)+1):(M*(p+3)-1)]
      alphas = c(alphas, 1-sum(alphas))
      if(StMAR==TRUE) {
        dfs = params[(M*(p+3)):(M*(p+4)-1)] # degrees of freedom
      } else if(GStMAR==TRUE) {
        dfs = params[(M*(p+3)):(M*(p+3)+M2-1)]
      }
    } else {
      # If restricted TRUE: transform the restricted parameter vector into the standard form. Then do everything like for non-restricted.
      phi0 = params[1:M]
      arcoefs = matrix(rep(params[(M+1):(M+p)], M), ncol=M)
      variances = params[(M+p+1):(p+2*M)]
      pars = rbind(phi0, arcoefs, variances)
      alphas = params[(p+2*M+1):(3*M+p-1)]
      if(StMAR==TRUE) {
        dfs = params[(3*M+p):(4*M+p-1)] # degrees of freedom
        params = c(as.vector(pars), alphas, dfs)
      } else if(GStMAR==TRUE) {
        dfs = params[(3*M+p):(3*M+p-1+M2)]
        params = c(as.vector(pars), alphas, dfs)
      } else {
        params = c(as.vector(pars), alphas)
      }
      alphas = c(alphas, 1-sum(alphas))
    }
  }
  rownames(pars) = NULL
  sigmas = pars[p+2,]

  if(checks==TRUE) {
    data = checkAndCorrectData(data=data, p=p)
    parameterChecks(p=p, M=M_orig, params=params, pars=pars, alphas=alphas, StMAR=StMAR, GStMAR=GStMAR, constraints=constraints)
  }
  n_obs = length(data)

  #### Start evaluating the mixing weights ####

  # Expected values mu_m (Kalliovirta 2015, s.250)
  mu = vapply(1:M, function(i1) pars[1, i1]/(1-sum(pars[2:(p+1), i1])), numeric(1))

  # Observed data: y_(-p+1),...,y_0,y_1,...,y_(n_obs-p). First row denotes vector y_0, i:th row vector y_[i-1] and last row denotes the vector y_T.
  Y = vapply(1:p, function(i1) data[(p-i1+1):(n_obs-i1+1)], numeric(n_obs-p+1) )

  # Calculate inverse Gamma_m and calculate the matrix products in mv normal and t-distribution (Galbraith and Galbraith 1974)
  matProd = matrix(nrow=n_obs-p+1, ncol=M)
  invG = array(dim=c(p, p, M))
  if(p==1) { # Form vectorized Gamma_m (Lutkepohl (2005), s.15-29)
    for(i1 in 1:M) {
      A = pars[p+1, i1]
      Sigma = as.matrix(sigmas[i1])
      VecGamma = solve(1-kronecker(A,A), Sigma)
      invG[,,i1] = as.matrix(1/VecGamma)
      matProd[,i1] = (Y-mu[i1])*invG[,,i1]*(Y-mu[i1])
    }
  } else { # Inverse formula by Galbraith, R., Galbraith, J., (1974)
    for(i1 in 1:M) {
      ARcoefs = pars[2:(p+1), i1]
      U = diag(1, p, p)
      V = diag(ARcoefs[p], p, p)
      for(i2 in 1:(p-1)) {
        U[(i2+1):p, i2] <- -ARcoefs[1:(p-i2)]
        V[(i2+1):p, i2] <- rev(ARcoefs[i2:(p-1)])
      }
      invG[,,i1] = (crossprod(U, U) - crossprod(V, V))/sigmas[i1]
      matProd[,i1] = rowSums((Y-mu[i1]*rep(1,p))%*%invG[,,i1]*(Y-mu[i1]*rep(1,p)))
    }
  }

  # Calculate the log multivariate normal or student's t values (Kalliovirta 2015, s250 eq.(7) or MPS 2018) for each vector y_t and for each m=1,..,M
  # First row for initial values y_0 (as denoted by Kalliovirta 2015) and i:th row for y_(i-1). First column for component m=1 and j:th column for m=j.
  logmv_values = matrix(nrow=(n_obs-p+1), ncol=M)
  if(StMAR==FALSE & GStMAR==FALSE) {
    for(i1 in 1:M) {
      detG = 1/det(as.matrix(invG[,,i1]))
      logmv_values[,i1] = -0.5*p*log(2*pi)-0.5*log(detG)-0.5*matProd[,i1]
    }
  } else if(StMAR==TRUE){
    for(i1 in 1:M) {
      detG = 1/det(as.matrix(invG[,,i1]))
      logC = lgamma(0.5*(p+dfs[i1]))-0.5*p*log(pi)-0.5*p*log(dfs[i1]-2)-lgamma(0.5*dfs[i1])
      logmv_values[,i1] = logC - 0.5*log(detG) - 0.5*(p+dfs[i1])*log(1 + matProd[,i1]/(dfs[i1]-2))
    }
  } else {  # If GStMAR==TRUE
    for(i1 in 1:M) {
      detG = 1/det(as.matrix(invG[,,i1]))
      if(i1 <= M1) { # Multinormals
        logmv_values[,i1] = -0.5*p*log(2*pi)-0.5*log(detG)-0.5*matProd[,i1]
      } else { # Multistudents
        logC = lgamma(0.5*(p+dfs[i1-M1]))-0.5*p*log(pi)-0.5*p*log(dfs[i1-M1]-2)-lgamma(0.5*dfs[i1-M1])
        logmv_values[,i1] = logC - 0.5*log(detG) - 0.5*(p+dfs[i1-M1])*log(1 + matProd[,i1]/(dfs[i1-M1]-2))
      }
    }
  }

  # Calculate the alpha_mt mixing weights (Kalliovirta 2015, s.250 eq.(8))
  # First row for t=1, second for t=2, and i:th for t=i. First column for m=1, second for m=2 and j:th column for m=j.
  logmv_values0 = logmv_values[1:(n_obs-p),] # The last row is not needed because alpha_mt uses vector Y_(t-1)
  if(!is.matrix(logmv_values0)) logmv_values0 = as.matrix(logmv_values0)

  if(M==1) {
    alpha_mt = as.matrix(rep(1, n_obs-p))
  } else if(any(logmv_values0 < epsilon)) { # Close to zero values handled with Brobdingnag if needed
    numerators = lapply(1:M, function(i1) alphas[i1]*Brobdingnag::as.brob(exp(1))^logmv_values0[,i1]) # alphas[i1]*exp( as.brob(mvn_values[,i1]) )
    denominator = Reduce("+", numerators) # For all t=0,...,T
    alpha_mt = vapply(1:M, function(i1) as.numeric(numerators[[i1]]/denominator), numeric(n_obs-p))
  } else {
    mv_values0 = exp(logmv_values0)
    denominator = colSums(alphas*t(mv_values0))
    alpha_mt = t(alphas*t(mv_values0/denominator))
  }
  return(alpha_mt)
}


#' @import stats
#'
#' @title Calculate mixing weights of GMAR, StMAR or G-StMAR model
#'
#' @description \code{mixingWeights} calculates the mixing weights of the specified GMAR, StMAR or G-StMAR model and returns them as a matrix.
#'
#' @inheritParams mixingWeights_int
#' @inherit mixingWeights_int return references
#' @examples
#' # GMAR model
#' params12 <- c(1.1, 0.9, 0.3, 4.5, 0.7, 3.2, 0.8)
#' mw12 <- mixingWeights(VIX, 1, 2, params12)
#'
#' # Restricted GMAR model
#' params12r <- c(1.4, 1.8, 0.9, 0.3, 3.2, 0.8)
#' mw12r <- mixingWeights(VIX, 1, 2, params12r, restricted=TRUE)
#'
#' # StMAR model
#' params12t <- c(1.1, 0.9, 0.3, 4.5, 0.7, 3.2, 0.8, 5, 8)
#' mw12t <- mixingWeights(VIX, 1, 2, params12t, StMAR=TRUE)
#'
#' # Non-mixture version of StMAR model
#' params11t <- c(0.76, 0.93, 1.35, 2.4)
#' mw11t <- mixingWeights(VIX, 1, 1, params11t, StMAR=TRUE)
#'
#' # G-StMAR model
#' params12gs <- c(1.2, 0.8, 0.6, 1.3, 0.6, 1.1, 0.6, 3)
#' mw12gs <- mixingWeights(VIX, 1, c(1,1), params12gs, GStMAR=TRUE)
#'
#' # Restricted G-StMAR model
#' params13gsr <- c(1.3, 2.2, 1.4, 0.8, 2.4, 4.6, 0.4, 0.25, 0.15, 20)
#' mw13gsr <- mixingWeights(VIX, 1, c(2,1), params13gsr, GStMAR=TRUE, restricted=TRUE)
#'
#' # GMAR model as a mixture of AR(2) and AR(1) models
#' R <- list(diag(1, ncol=2, nrow=2), as.matrix(c(1, 0)))
#' params22c <- c(1.2, 0.8, 0.1, 0.3, 3.3, 0.8, 2.8, 0.8)
#' mw22c <- mixingWeights(VIX, 2, 2, params22c, constraints=TRUE, R=R)
#'
#' # Such StMAR(3,2) that the AR coefficients are restricted to be
#' # the same for both regimes and that the second AR coefficients are
#' # constrained to zero.
#' params32trc <- c(2.2, 1.8, 0.88, -0.03, 2.4, 0.27, 0.40, 3.9, 1000)
#' mw32trc <- mixingWeights(VIX, 3, 2, params32trc, StMAR=TRUE,
#'                          restricted=TRUE, constraints=TRUE,
#'                          R=matrix(c(1, 0, 0, 0, 0, 1), ncol=2))
#' @export

mixingWeights <- function(data, p, M, params, StMAR=FALSE, GStMAR=FALSE, restricted=FALSE, constraints=FALSE, R) {
  checkLogicals(StMAR=StMAR, GStMAR=GStMAR)
  checkPM(p=p, M=M, GStMAR=GStMAR)
  if(length(params)!=nParams(p=p, M=M, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted, constraints=constraints, R=R)) {
    stop("The parameter vector has wrong dimension")
  }
  return(mixingWeights_int(data=data, p=p, M=M, params=params, StMAR=StMAR, GStMAR=GStMAR, restricted=restricted, constraints=constraints, R=R, checks=TRUE))
}

#' @import stats
#'
#' @title Compute the log-likelihood of GMAR, StMAR or G-StMAR model
#'
#' @description \code{loglikelihood_int} computes the log-likelihood value of the specified GMAR, StMAR or G-StMAR model for the given data.
#'
#' @param data a numeric vector class \code{'ts'} object containing the data. \code{NA} values are not supported.
#' @param p a positive integer specifying the order of AR coefficients.
#' @param M \describe{
#'   \item{For \strong{GMAR} and \strong{StMAR} models:}{a positive integer specifying the number of mixture components.}
#'   \item{For \strong{G-StMAR} model:}{a size (2x1) vector specifying the number of \emph{GMAR-type} components \code{M1} in the
#'    first element and \emph{StMAR-type} components \code{M2} in the second. The total number of mixture components is \code{M=M1+M2}.}
#' }
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
#'          matrices \strong{C} that satisfy \strong{\eqn{\phi_{m}}}\eqn{=}\strong{\eqn{R_{m}\psi_{m}}} for all \eqn{m=1,...,M}, where
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
#'          \strong{\eqn{C}} that satisfies \strong{\eqn{\phi}}\eqn{=}\strong{\eqn{R\psi}}, where
#'          \strong{\eqn{\psi}}\eqn{=(\psi_{1},...,\psi_{q})}.}
#'      }
#'    }
#'  }
#'  Symbol \eqn{\phi} denotes an AR coefficient, \eqn{\sigma^2} a variance, \eqn{\alpha} a mixing weight and \eqn{\nu} a degrees of
#'  freedom parameter. If \code{parametrization=="mean"} just replace each intercept term \eqn{\phi_{m,0}} with regimewise mean
#'  \eqn{\mu_m = \phi_{m,0}/(1-\sum\phi_{i,m})}. In the \strong{G-StMAR} model the first \code{M1} components are \emph{GMAR-type}
#'  and the rest \code{M2} components are \emph{StMAR-type}.
#'  Note that in the case \strong{M=1} the parameter \eqn{\alpha} is dropped, and in the case of \strong{StMAR} or \strong{G-StMAR} model
#'  the degrees of freedom parameters \eqn{\nu_{m}} have to be larger than \eqn{2}.
#' @param restricted a logical argument stating whether the AR coefficients \eqn{\phi_{m,1},...,\phi_{m,p}} are restricted
#'  to be the same for all regimes.
#' @param model is "GMAR", "StMAR" or "G-StMAR" model considered? In G-StMAR model the first \code{M1} components
#'  are \emph{GMAR-type} and the rest \code{M2} components are \emph{StMAR-type}.
#' @param constraints specifies linear constraints applied to the autoregressive parameters.
#'   \describe{
#'   \item{For \strong{non-restricted} models:}{a list of size \eqn{(pxq_{m})} constraint matrices \strong{\eqn{C_{m}}} of full column rank
#'     satisfying \strong{\eqn{\phi_{m}}}\eqn{=}\strong{\eqn{C_{m}\psi_{m}}} for all \eqn{m=1,...,M}, where
#'     \strong{\eqn{\phi_{m}}}\eqn{=(\phi_{m,1},...,\phi_{m,p})} and \strong{\eqn{\psi_{m}}}\eqn{=(\psi_{m,1},...,\psi_{m,q_{m}})}.}
#'   \item{For \strong{restricted} models:}{a size \eqn{(pxq)} constraint matrix \strong{\eqn{C}} of full column rank satisfying
#'     \strong{\eqn{\phi}}\eqn{=}\strong{\eqn{C\psi}}, where \strong{\eqn{\phi}}\eqn{=(\phi_{1},...,\phi_{p})} and
#'     \strong{\eqn{\psi}}\eqn{=\psi_{1},...,\psi_{q}}.}
#'   }
#'   Symbol \eqn{\phi} denotes an AR coefficient. Note that regardless of any constraints, the nominal order of AR coefficients is
#'   alway \code{p} for all regimes.
#'   Ignore or set to \code{NULL} if applying linear constraints is \strong{not} desired.
#' @param conditional a logical argument specifying whether the conditional or exact log-likehood function should be used.
#' @param parametrization is the model parametrized with the "intercepts" \eqn{\phi_{m,0}} or
#'  "means" \eqn{\mu_m = \phi_{m,0}/(1-\sum\phi_{i,m})}?
#' @param boundaries a logical argument. If \code{TRUE} then \code{loglikelihood} returns \code{minval} if...
#' \itemize{
#'   \item any component variance is not larger than zero,
#'   \item any parametrized mixing weight \eqn{\alpha_{1},...,\alpha_{M-1}} is not larger than zero,
#'   \item sum of the parametrized mixing weights is not smaller than one,
#'   \item if the model is not stationary,
#'   \item or if \code{model=="StMAR"} or \code{model=="G-StMAR"} and any degrees of freedom parameter \eqn{\nu_{m}} is not larger than two.
#' }
#' Argument \code{minval} will be used only if \code{boundaries==TRUE}.
#' @param checks an (optional) logical argument defining whether argument checks are made. If \code{FALSE} then no argument checks
#' such as stationary checks etc are made. The default is \code{TRUE}.
#' @param to_return should the returned object be the log-likelihood value, mixing weights, mixing weights including
#'   value for \eqn{alpha_{m,T+1}}, a list containing log-likelihood value and mixing weights, the terms \eqn{l_{t}: t=1,..,T}
#'   in the log-likelihood function (see \emph{KMS 2015, eq.(13)}), regimewise conditional means, regimewise conditional variances,
#'   total conditional means or total conditional variances?
#'   Default is the log-likelihood value (\code{"loglik"}).
#' @param minval this will be returned when the parameter vector is outside the parameter space.
#' @return
#'  Note that the first p observations are taken as the initial values so
#'  mixing weights and conditional moments start from the p+1:th observation (interpreted as t=1).
#'  \describe{
#'   \item{By default:}{log-likelihood value of the specified model,}
#'   \item{If \code{to_return=="mw"}:}{a size ((n_obs-p)xM) matrix containing the mixing weights: for m:th component in m:th column.}
#'   \item{If \code{to_return=="mw_tplus1"}:}{a size ((n_obs-p+1)xM) matrix containing the mixing weights: for m:th component in m:th column.
#'     The last row is for \eqn{\alpha_{m,T+1}}.}
#'   \item{if \code{to_return=="loglik_and_mw"}:}{a list of two elements. The first element contains the log-likelihood value and the
#'     second element contains the mixing weights.}
#'   \item{If \code{to_return=="terms"}:}{a size ((n_obs-p)x1) numeric vector containing the terms \eqn{l_{t}}.}
#'   \item{if \code{to_return=="regime_cmeans"}:}{a size ((n_obs-p)xM) matrix containing the regime specific conditional means.}
#'   \item{if \code{to_return=="regime_cvars"}:}{a size ((n_obs-p)xM) matrix containing the regime specific conditional variances.}
#'   \item{if \code{to_return=="total_cmeans"}:}{a size ((n_obs-p)x1) vector containing the total conditional means.}
#'   \item{if \code{to_return=="total_cvars"}:}{a size ((n_obs-p)x1) vector containing the total conditional variances.}
#'  }
#' @references
#'  \itemize{
#'    \item Galbraith, R., Galbraith, J. 1974. On the inverses of some patterned matrices arising
#'            in the theory of stationary time series. \emph{Journal of Applied Probability} \strong{11}, 63-71.
#'    \item Kalliovirta L., Meitz M. and Saikkonen P. 2015. Gaussian Mixture Autoregressive model for univariate time series.
#'            \emph{Journal of Time Series Analysis}, \strong{36}, 247-266.
#'    \item Meitz M., Preve D., Saikkonen P. 2018. A mixture autoregressive model based on Student's t-distribution.
#'            arXiv:1805.04010 \strong{[econ.EM]}.
#'    \item There are currently no published references for the G-StMAR model, but it's a straightforward generalization with
#'            theoretical properties similar to the GMAR and StMAR models.
#'  }

loglikelihood_int <- function(data, p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL,
                              conditional=TRUE, parametrization=c("intercept", "mean"), boundaries=TRUE, checks=TRUE,
                              to_return=c("loglik", "mw", "mw_tplus1", "loglik_and_mw", "terms", "regime_cmeans", "regime_cvars",
                                          "total_cmeans", "total_cvars"), minval) {
  epsilon <- round(log(.Machine$double.xmin) + 10)
  model <- match.arg(model)
  parametrization <- match.arg(parametrization)
  to_return <- match.arg(to_return)
  M_orig <- M
  if(model == "G-StMAR") {
    M1 <- M[1]
    M2 <- M[2]
    M <- sum(M)
  }

  # Reform parameters to the "standard form" and collect them
  if(checks == TRUE) {
    checkConstraintMat(p=p, M=M_orig, restricted=restricted, constraints=constraints)
  }
  params <- removeAllConstraints(p=p, M=M_orig, params=params, model=model, restricted=restricted, constraints=constraints)
  pars <- pick_pars(p=p, M=M_orig, params=params, model=model, restricted=FALSE, constraints=NULL)
  alphas <- pick_alphas(p=p, M=M_orig, params=params, model=model, restricted=FALSE, constraints=NULL)
  dfs <- pick_dfs(p=p, M=M_orig, params=params, model=model)
  sigmas <- pars[p + 2,]

  # Return minval if parameters are out of their boundaries.
  if(boundaries == TRUE) {
    if(any(pars[p + 2,] <= 0)) {
      return(minval)
    } else if(M >= 2 & sum(alphas[-M]) >= 1) {
      return(minval)
    } else if(any(alphas <= 0)) {
      return(minval)
    } else if(!isStationary_int(p, M, params, restricted=FALSE)) {
      return(minval)
    }
    if(model == "StMAR" | model == "G-StMAR") {
      if(any(dfs <= 2 + 1e-8 | dfs > 1e+6)) return(minval)
    }
  }

  if(checks == TRUE) {
    data <- checkAndCorrectData(data=data, p=p)
    parameterChecks(p=p, M=M_orig, params=params, model=model, restricted=FALSE, constraints=NULL)
  }
  n_obs <- length(data)

  #### Start evaluating the log-likelihood ####

  # Expected values mu_m (Kalliovirta 2015, s.250)
  if(parametrization == "mean") {
    mu <- pars[1,]
    pars[1,] <- vapply(1:M, function(i1) mu[i1]*(1 - sum(pars[2:(p + 1), i1])), numeric(1))
  } else {
    mu <- vapply(1:M, function(i1) pars[1, i1]/(1 - sum(pars[2:(p + 1), i1])), numeric(1))
  }

  # Observed data: y_(-p+1),...,y_0,y_1,...,y_(n_obs-p). First row denotes vector y_0, i:th row vector y_[i-1] and last row denotes the vector y_T.
  Y <- vapply(1:p, function(i1) data[(p - i1 + 1):(n_obs - i1 + 1)], numeric(n_obs - p + 1) )

  # Calculate inverse Gamma_m and calculate the matrix products in mv normal and t-distribution (Galbraith and Galbraith 1974)
  matProd <- matrix(nrow=n_obs - p + 1, ncol=M)
  invG <- array(dim=c(p, p, M))
  if(p == 1) { # Inverse formula by Galbraith, R., Galbraith, J., (1974)
    for(i1 in 1:M) {
      invG[, , i1] <- (1 - pars[p + 1, i1]^2)/sigmas[i1]
      matProd[, i1] <- (Y - mu[i1])*invG[, , i1]*(Y - mu[i1])
    }
  } else {
    for(i1 in 1:M) {
      ARcoefs <- pars[2:(p + 1), i1]
      U <- diag(1, nrow=p, ncol=p)
      V <- diag(ARcoefs[p], nrow=p, ncol=p)
      for(i2 in 1:(p - 1)) {
        U[(i2 + 1):p, i2] <- -ARcoefs[1:(p - i2)]
        V[(i2 + 1):p, i2] <- rev(ARcoefs[i2:(p - 1)])
      }
      invG[, , i1] <- (crossprod(U, U) - crossprod(V, V))/sigmas[i1]
      matProd[, i1] <- rowSums((Y - mu[i1]*rep(1, p))%*%invG[, , i1]*(Y - mu[i1]*rep(1, p)))
    }
  }

  # Calculate the log multivariate normal or student's t values (Kalliovirta 2015, s250 eq.(7) or MPS 2018) for each vector y_t and for each m=1,..,M
  # First row for initial values y_0 (as denoted by Kalliovirta 2015) and i:th row for y_(i-1). First column for component m=1 and j:th column for m=j.
  logmv_values <- matrix(nrow=(n_obs - p + 1), ncol=M)
  if(model == "GMAR") {
    for(i1 in 1:M) {
      detG <- 1/det(as.matrix(invG[, , i1]))
      logmv_values[,i1] <- -0.5*p*log(2*pi) - 0.5*log(detG) - 0.5*matProd[,i1]
    }
  } else if(model == "StMAR"){
    for(i1 in 1:M) {
      detG <- 1/det(as.matrix(invG[, , i1]))
      logC <- lgamma(0.5*(p + dfs[i1])) - 0.5*p*log(pi) - 0.5*p*log(dfs[i1] - 2) - lgamma(0.5*dfs[i1])
      logmv_values[,i1] <- logC - 0.5*log(detG) - 0.5*(p + dfs[i1])*log(1 + matProd[,i1]/(dfs[i1] - 2))
    }
  } else {  # If model == "G-StMAR
    for(i1 in 1:M) {
      detG <- 1/det(as.matrix(invG[, , i1]))
      if(i1 <= M1) { # Multinormals
        logmv_values[,i1] <- -0.5*p*log(2*pi) - 0.5*log(detG) - 0.5*matProd[,i1]
      } else { # Multistudents
        logC <- lgamma(0.5*(p + dfs[i1 - M1])) - 0.5*p*log(pi) - 0.5*p*log(dfs[i1 - M1] - 2) - lgamma(0.5*dfs[i1 - M1])
        logmv_values[,i1] <- logC - 0.5*log(detG) - 0.5*(p + dfs[i1 - M1])*log(1 + matProd[,i1]/(dfs[i1 - M1] - 2))
      }
    }
  }

  # Calculate the alpha_mt mixing weights (Kalliovirta 2015, s.250 eq.(8))
  # First row for t=1, second for t=2, and i:th for t=i. First column for m=1, second for m=2 and j:th column for m=j.
  if(to_return != "mw_tplus1") {
    logmv_values0 <- logmv_values[1:(n_obs - p),] # The last row is not needed because alpha_mt uses vector Y_(t-1)
  } else {
    logmv_values0 <- logmv_values # The last row is needed for alpha_{m,t+1}
  }
  if(!is.matrix(logmv_values0)) logmv_values0 <- as.matrix(logmv_values0)

  l_0 <- 0 # "The first term" of the exact log-likelihood
  if(M == 1) {
    alpha_mt <- as.matrix(rep(1, nrow(logmv_values0)))
    if(conditional == FALSE) { # Calculate "the first term" of the log-likelihood (Kalliovirta ym 2015, s.254 eq.(12))
      l_0 <- logmv_values[1]
    }
  } else if(any(logmv_values0 < epsilon)) { # Close to zero values handled with Brobdingnag if needed
    numerators <- lapply(1:M, function(i1) alphas[i1]*exp(Brobdingnag::as.brob(logmv_values0[,i1]))) # alphas[i1]*Brobdingnag::as.brob(exp(1))^logmv_values0[,i1]
    denominator <- Reduce("+", numerators) # For all t=0,...,T
    alpha_mt <- vapply(1:M, function(i1) as.numeric(numerators[[i1]]/denominator), numeric(nrow(logmv_values0)))

    if(conditional == FALSE) {
      l_0 <- log(Reduce("+", lapply(1:M, function(i1) numerators[[i1]][1])))
    }
  } else {
    mv_values0 <- exp(logmv_values0)
    denominator <- as.vector(mv_values0%*%alphas)
    alpha_mt <- (mv_values0/denominator)%*%diag(alphas)

    if(conditional == FALSE) {
      l_0 <- log(sum(alphas*mv_values0[1,]))
    }
  }
  if(to_return == "mw" | to_return == "mw_tplus1") {
    return(alpha_mt)
  }

  # Calculate the conditional means mu_mt (Kalliovirta 2015, s.249 eq.(2)). First row for t=1, second for t=2 etc. First column for m=1, second column for m=2 etc.
  if(p == 1) {
    mu_mt <- vapply(1:M, function(i1) rep(pars[1, i1], nrow(Y) - 1) + Y[1:(nrow(Y) - 1),]*pars[2, i1], numeric(n_obs - p))
  } else {
    mu_mt <- vapply(1:M, function(i1) rep(pars[1, i1], nrow(Y) - 1) + colSums(pars[2:(p + 1), i1]*t(Y[1:(nrow(Y) - 1),])), numeric(n_obs - p))
  }

  # Calculate/return conditional means
  if(to_return == "regime_cmeans") {
    return(mu_mt)
  } else if(to_return == "total_cmeans") {
    return(rowSums(alpha_mt*mu_mt))
  }

  # Calculate "the second term" of the log-likelihood (Kalliovirta 2015, s.254 eq.(12)-(13) )
  Y2 <- Y[2:nrow(Y), 1] # Only first column and rows 2...T are needed
  if(model == "GMAR") {
    invsqrt_sigmas <- sigmas^(-1/2)
    if(M == 1) {
      lt_tmp <- invsqrt_sigmas*dnorm((Y2 - mu_mt)*invsqrt_sigmas)
    } else {
      lt_tmp <- alpha_mt*dnorm((Y2 - mu_mt)%*%diag(invsqrt_sigmas))%*%diag(invsqrt_sigmas)
    }
  } else if(model == "StMAR") {
    matProd0 <- matProd[1:(n_obs - p),] # Last row is not needed because sigma_t uses y_{t-1}
    if(M == 1) {
      sigma_mt <- sigmas*(dfs - 2 + matProd0)/(dfs - 2 + p) # Conditional variances
      lt_tmp <- ((exp(lgamma(0.5*(1 + dfs + p)) - lgamma(0.5*(dfs + p)))/sqrt(pi*(dfs + p - 2)))/sqrt(sigma_mt))*(1 + ((Y2 - mu_mt)^2)/((dfs + p - 2)*sigma_mt))^(-0.5*(1 + dfs + p))
    } else {
      sigma_mt <- t(dfs - 2 + t(matProd0))%*%diag(1/(dfs - 2 + p))%*%diag(sigmas)
      lt_tmp <- alpha_mt*t(exp(lgamma(0.5*(1 + dfs + p)) - lgamma(0.5*(dfs + p)))/sqrt(pi*(dfs + p - 2))/t(sqrt(sigma_mt)))*t(t(1 + ((Y2 - mu_mt)^2)/(sigma_mt%*%diag(dfs + p - 2)))^(-0.5*(1 + dfs + p)))
    }
  } else {  # If model == "G-StMAR"
    # GMAR-components
    invsqrt_sigmasM1 <- sigmas[1:M1]^(-1/2)
    if(M1 == 1) {
      lt_tmpM1 <- alpha_mt[,1]*invsqrt_sigmasM1*dnorm((Y2 - mu_mt[,1])*invsqrt_sigmasM1)
    } else {
      lt_tmpM1 <- alpha_mt[,1:M1]*dnorm((Y2 - mu_mt[,1:M1])%*%diag(invsqrt_sigmasM1))%*%diag(invsqrt_sigmasM1)
    }
    # StMAR-components
    sigmasM2 <- sigmas[(M1 + 1):M]
    matProd0 <- matProd[1:(n_obs - p), (M1 + 1):M]
    if(M2 == 1) {
      sigma_mt <- sigmasM2*(dfs - 2 + matProd0)/(dfs - 2 + p) # Conditional variances
      lt_tmpM2 <- alpha_mt[,(M1 + 1):M]*((exp(lgamma(0.5*(1 + dfs + p)) - lgamma(0.5*(dfs + p)))/sqrt(pi*(dfs + p - 2)))/sqrt(sigma_mt))*(1 + ((Y2-mu_mt[,(M1 + 1):M])^2)/((dfs + p - 2)*sigma_mt))^(-0.5*(1 + dfs + p))
    } else {
      sigma_mt <- t(dfs - 2 + t(matProd0))%*%diag(1/(dfs - 2 + p))%*%diag(sigmasM2)
      lt_tmpM2 <- alpha_mt[,(M1 + 1):M]*t(exp(lgamma(0.5*(1 + dfs + p))-lgamma(0.5*(dfs + p)))/sqrt(pi*(dfs + p - 2))/t(sqrt(sigma_mt)))*t(t(1 + ((Y2 - mu_mt[,(M1 + 1):M])^2)/(sigma_mt%*%diag(dfs + p - 2)))^(-0.5*(1 + dfs + p)))
    }
    lt_tmp <- cbind(lt_tmpM1, lt_tmpM2)
  }
  l_t <- rowSums(lt_tmp)

  # Calculate/return conditional variances
  if(to_return == "regime_cvars" | to_return == "total_cvars") {
    if(model == "GMAR") {
      sigma_mt <- matrix(rep(sigmas, n_obs - p), ncol=M, byrow=TRUE)
    } else if(model == "StMAR") {
      sigma_mt <- as.matrix(sigma_mt)
    } else if(model == "G-StMAR") {
      sigma_mt1 <- matrix(rep(sigmas[1:M1], n_obs - p), ncol=M1, byrow=TRUE)
      sigma_mt <- cbind(sigma_mt1, sigma_mt)
      colnames(sigma_mt) <- NULL
    }
    if(to_return == "regime_cvars") {
      return(sigma_mt)
    } else { # Calculate and return total conditional variances
      return(rowSums(alpha_mt*sigma_mt) - rowSums(alpha_mt*(mu_mt - rowSums(alpha_mt*mu_mt))^2))
    }
  }

  if(to_return == "terms") {
    ret <- log(l_t)
  } else if(to_return == "loglik_and_mw") {
    ret <- list(loglik=l_0 + sum(log(l_t)), mw=alpha_mt)
  } else {
    ret <- l_0 + sum(log(l_t))
  }
  ret
}


#' @import stats
#'
#' @title Compute the log-likelihood of GMAR, StMAR or G-StMAR model
#'
#' @description \code{loglikelihood} computes the log-likelihood value of the specified GMAR, StMAR or G-StMAR model.
#'   Exists for convenience if one wants to for example plot profile log-likelihoods or employ other estimation algorithms.
#'   Use \code{minval} to control what happens when the parameter vector is outside the parameter space.
#'
#' @inheritParams loglikelihood_int
#' @param returnTerms should the terms \eqn{l_{t}: t=1,..,T} in the log-likelihood function (see \emph{KMS 2015, eq.(13)})
#'   be returned instead of the log-likelihood value?
#' @return Returns the log-likelihood value or the terms described in \code{returnTerms}.
#' @inherit loglikelihood_int references
#' @seealso \code{\link{fitGSMAR}}, \code{\link{GSMAR}}, \code{\link{quantileResiduals}},
#'  \code{\link{mixingWeights}}, \code{\link{calc_gradient}}
#' @examples
#' # GMAR model
#' params12 <- c(0.18, 0.93, 0.01, 0.86, 0.68, 0.02, 0.88)
#' loglikelihood(logVIX, 1, 2, params12)
#'
#' # Restricted GMAR model, outside parameter space
#' params12r <- c(0.21, 0.23, 0.92, 0.01, 0.02, 0.86)
#' loglikelihood(logVIX, 1, 2, params12r, restricted=TRUE)
#'
#' # Non-mixture version of StMAR model, outside parameter space
#' params11t <- c(0.16, 0.93, 0.00, 3.01)
#' loglikelihood(logVIX, 1, 1, params11t, model="StMAR", minval="Hello")
#'
#' # G-StMAR model
#' params12gs <- c(0.86, 0.68, 0.02, 0.18, 0.93, 0.01, 0.11, 44.36)
#' loglikelihood(logVIX, 1, c(1, 1), params12gs, model="G-StMAR")
#'
#' # Restricted G-StMAR model
#' params12gsr <- c(0.31, 0.33, 0.88, 0.01, 0.02, 0.77, 2.72)
#' loglikelihood(logVIX, 1, c(1, 1), params12gsr, model="G-StMAR",
#'  restricted=TRUE)
#'
#' # GMAR model as a mixture of AR(2) and AR(1) models
#' constraints <- list(diag(1, ncol=2, nrow=2), as.matrix(c(1, 0)))
#' params22c <- c(0.61, 0.83, -0.06, 0.02, 0.21, 0.91, 0.01, 0.16)
#' loglikelihood(logVIX, 2, 2, params22c, constraints=constraints)
#'
#' # Such StMAR(3,2) that the AR coefficients are restricted to be
#' # the same for both regimes and that the second AR coefficients are
#' # constrained to zero.
#' params32trc <- c(0.35, 0.33, 0.88, -0.02, 0.01, 0.01, 0.36, 4.53, 1000)
#' loglikelihood(logVIX, 3, 2, params32trc, model="StMAR", restricted=TRUE,
#'               constraints=matrix(c(1, 0, 0, 0, 0, 1), ncol=2))
#' @export

loglikelihood <- function(data, p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL,
                          conditional=TRUE, parametrization=c("intercept", "mean"), returnTerms=FALSE, minval=NA) {
  model <- match.arg(model)
  check_model(model)
  parametrization <- match.arg(parametrization)
  stopifnot(parametrization %in% c("intercept", "mean"))
  checkPM(p=p, M=M, model=model)
  check_params_length(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints)
  to_ret <- ifelse(returnTerms, "terms", "loglik")
  loglikelihood_int(data=data, p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints,
                    conditional=conditional, parametrization=parametrization, boundaries=TRUE, checks=FALSE,
                    to_return=to_ret, minval=minval)
}


#' @import stats
#'
#' @title Calculate mixing weights of GMAR, StMAR or G-StMAR model
#'
#' @description \code{mixingWeights_int} calculates the mixing weights of the specified GMAR, StMAR or G-StMAR model and returns them as a matrix.
#'
#' @inheritParams loglikelihood_int
#' @param to_return should the returned object the mixing weights or mixing weights (\code{"mw"}) including
#'   value for \eqn{alpha_{m,T+1}} (\code{"mw_tplus1"})?
#' @return
#'  \describe{
#'   \item{If \code{to_return=="mw"}:}{a size ((n_obs-p)xM) matrix containing the mixing weights: for m:th component in m:th column.}
#'   \item{If \code{to_return=="mw_tplus1"}:}{a size ((n_obs-p+1)xM) matrix containing the mixing weights: for m:th component in m:th column.
#'     The last row is for \eqn{\alpha_{m,T+1}}}.
#'  }
#' @inherit loglikelihood_int references

mixingWeights_int <- function(data, p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL,
                              parametrization=c("intercept", "mean"), checks=TRUE, to_return=c("mw", "mw_tplus1")) {
  to_ret <- match.arg(to_return)
  parametrization <- match.arg(parametrization)
  loglikelihood_int(data=data, p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints,
                    parametrization=parametrization, boundaries=FALSE, checks=checks, to_return=to_ret, minval=NA)
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
#' params12 <- c(0.18, 0.93, 0.01, 0.86, 0.68, 0.02, 0.88)
#' mixingWeights(logVIX, 1, 2, params12)
#'
#' # Restricted GMAR model, outside parameter space
#' params12r <- c(0.21, 0.23, 0.92, 0.01, 0.02, 0.86)
#' mixingWeights(logVIX, 1, 2, params12r, restricted=TRUE)
#'
#' # Non-mixture version of StMAR model, outside parameter space
#' params11t <- c(0.16, 0.93, 0.01, 3.01)
#' mixingWeights(logVIX, 1, 1, params11t, model="StMAR")
#'
#' # G-StMAR model
#' params12gs <- c(0.86, 0.68, 0.02, 0.18, 0.93, 0.01, 0.11, 44.36)
#' mixingWeights(logVIX, 1, c(1, 1), params12gs, model="G-StMAR")
#'
#' # Restricted G-StMAR model
#' params12gsr <- c(0.31, 0.33, 0.88, 0.01, 0.02, 0.77, 2.72)
#' mixingWeights(logVIX, 1, c(1, 1), params12gsr, model="G-StMAR", restricted=TRUE)
#'
#' # GMAR model as a mixture of AR(2) and AR(1) models
#' constraints <- list(diag(1, ncol=2, nrow=2), as.matrix(c(1, 0)))
#' params22c <- c(0.61, 0.83, -0.06, 0.02, 0.21, 0.91, 0.01, 0.16)
#' mixingWeights(logVIX, 2, 2, params22c, constraints=constraints)
#'
#' # Such StMAR(3,2) that the AR coefficients are restricted to be
#' # the same for both regimes and that the second AR coefficients are
#' # constrained to zero.
#' params32trc <- c(0.35, 0.33, 0.88, -0.02, 0.01, 0.01, 0.36, 4.53, 1000)
#' mixingWeights(logVIX, 3, 2, params32trc, model="StMAR", restricted=TRUE,
#'               constraints=matrix(c(1, 0, 0, 0, 0, 1), ncol=2))
#' @export

mixingWeights <- function(data, p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL,
                          parametrization=c("intercept", "mean")) {
  model <- match.arg(model)
  check_model(model)
  parametrization <- match.arg(parametrization)
  stopifnot(parametrization %in% c("intercept", "mean"))
  checkPM(p=p, M=M, model=model)
  check_params_length(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints)
  mixingWeights_int(data=data, p=p, M=M, params=params, model=model, restricted=restricted,
                    constraints=constraints, parametrization=parametrization, checks=TRUE, to_return="mw")
}


#' @import stats
#'
#' @title Calculate conditional moments of GMAR, StMAR or G-StMAR model
#'
#' @description \code{condMoments} calculates the regime specific conditional means and variances and total
#'  conditional means and variances of the specified GMAR, StMAR or G-StMAR model.
#'
#' @inheritParams loglikelihood_int
#' @param to_return calculate regimewise conditional means (\code{regime_cmeans}), regimewise conditional variances
#'  (\code{regime_cvars}), total conditional means (\code{total_cmeans}), or total conditional variances (\code{total_cvars})?
#' @inherit loglikelihood_int references
#' @family moment functions
#' @return
#'  Note that the first p observations are taken as the initial values so the conditional moments
#'  start form the p+1:th observation (interpreted as t=1).
#'  \describe{
#'   \item{if \code{to_return=="regime_cmeans"}:}{a size ((n_obs-p)xM) matrix containing the regime specific conditional means.}
#'   \item{if \code{to_return=="regime_cvars"}:}{a size ((n_obs-p)xM) matrix containing the regime specific conditional variances.}
#'   \item{if \code{to_return=="total_cmeans"}:}{a size ((n_obs-p)x1) vector containing the total conditional means.}
#'   \item{if \code{to_return=="total_cvars"}:}{a size ((n_obs-p)x1) vector containing the total conditional variances.}
#'  }
#' @examples
#' # GMAR model
#' params12 <- c(0.18, 0.93, 0.01, 0.86, 0.68, 0.02, 0.88)
#' rcm12 <- condMoments(logVIX, 1, 2, params12, to_return="regime_cmeans")
#' rcv12 <- condMoments(logVIX, 1, 2, params12, to_return="regime_cvars")
#' tcm12 <- condMoments(logVIX, 1, 2, params12, to_return="total_cmeans")
#' tcv12 <- condMoments(logVIX, 1, 2, params12, to_return="total_cvars")
#'
#' # StMAR model
#' params12t <- c(0.17,  0.93, 0.01, 4.87, -0.90, 0.01, 0.98, 4.22, 1000)
#' rcm12t <- condMoments(logVIX, 1, 2, params12t, model="StMAR",
#'  to_return="regime_cmeans")
#' rcv12t <- condMoments(logVIX, 1, 2, params12t, model="StMAR",
#'  to_return="regime_cvars")
#'
#' # G-StMAR model
#' params12gs <- c(0.86, 0.68, 0.02, 0.18, 0.93, 0.01, 0.11, 44)
#' rcv12gs <- condMoments(logVIX, 1, c(1,1), params12gs, model="G-StMAR",
#'  to_return="regime_cvars")
#' tcv12gs <- condMoments(logVIX, 1, c(1,1), params12gs, model="G-StMAR",
#'  to_return="total_cvars")
#' @export
condMoments <- function(data, p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL,
                        parametrization=c("intercept", "mean"), to_return=c("regime_cmeans", "regime_cvars", "total_cmeans", "total_cvars")) {
  model <- match.arg(model)
  check_model(model)
  parametrization <- match.arg(parametrization)
  stopifnot(parametrization %in% c("intercept", "mean"))
  checkPM(p=p, M=M, model=model)
  check_params_length(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints)
  loglikelihood_int(data=data, p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints,
                    conditional=TRUE, parametrization=parametrization, boundaries=FALSE, checks=TRUE, to_return=to_return)
}


#' @title Calculate AIC, HQIC and BIC
#'
#' @description \code{get_IC} calculates AIC, HQIC and BIC
#'
#' @param loglik log-likelihood value
#' @param npars number of (freely estimated) parameters in the model.
#' @param obs numbers of observations with starting values excluded for conditional models.
#' @return Returns a data frame containing the information criteria values.

get_IC <- function(loglik, npars, obs) {
  AIC <- -2*loglik + 2*npars
  HQIC <- -2*loglik + 2*npars*log(log(obs))
  BIC <- -2*loglik + npars*log(obs)
  data.frame(AIC=AIC, HQIC=HQIC, BIC=BIC)
}

#' @import stats
#'
#' @title Compute the log-likelihood of GMAR, StMAR, or G-StMAR model
#'
#' @description \code{loglikelihood_int} computes the log-likelihood of the specified GMAR, StMAR, or G-StMAR model.
#'
#' @param data a numeric vector or class \code{'ts'} object containing the data. \code{NA} values are not supported.
#' @param p a positive integer specifying the autoregressive order of the model.
#' @param M \describe{
#'   \item{For \strong{GMAR} and \strong{StMAR} models:}{a positive integer specifying the number of mixture components.}
#'   \item{For \strong{G-StMAR} models:}{a size (2x1) integer vector specifying the number of \emph{GMAR type} components \code{M1} in the
#'    first element and \emph{StMAR type} components \code{M2} in the second element. The total number of mixture components is \code{M=M1+M2}.}
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
#'          \sigma_{1}^2,...,\sigma_{M}^2,\alpha_{1},...,\alpha_{M-1})}, where \strong{\eqn{\phi}}=\eqn{(\phi_{1},...,\phi_{p})}.}
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
#'  Symbol \eqn{\phi} denotes an AR coefficient, \eqn{\sigma^2} a variance, \eqn{\alpha} a mixing weight, and \eqn{\nu} a degrees of
#'  freedom parameter. If \code{parametrization=="mean"}, just replace each intercept term \eqn{\phi_{m,0}} with regimewise mean
#'  \eqn{\mu_m = \phi_{m,0}/(1-\sum\phi_{i,m})}. In the \strong{G-StMAR} model, the first \code{M1} components are \emph{GMAR type}
#'  and the rest \code{M2} components are \emph{StMAR type}.
#'  Note that in the case \strong{M=1}, the parameter \eqn{\alpha} is dropped, and in the case of \strong{StMAR} or \strong{G-StMAR} model,
#'  the degrees of freedom parameters \eqn{\nu_{m}} have to be larger than \eqn{2}.
#' @param restricted a logical argument stating whether the AR coefficients \eqn{\phi_{m,1},...,\phi_{m,p}} are restricted
#'  to be the same for all regimes.
#' @param model is "GMAR", "StMAR", or "G-StMAR" model considered? In the G-StMAR model, the first \code{M1} components
#'  are \emph{GMAR type} and the rest \code{M2} components are \emph{StMAR type}.
#' @param constraints specifies linear constraints applied to the autoregressive parameters.
#'   \describe{
#'   \item{For \strong{non-restricted} models:}{a list of size \eqn{(pxq_{m})} constraint matrices \strong{\eqn{C_{m}}} of full column rank
#'     satisfying \strong{\eqn{\phi_{m}}}\eqn{=}\strong{\eqn{C_{m}\psi_{m}}} for all \eqn{m=1,...,M}, where
#'     \strong{\eqn{\phi_{m}}}\eqn{=(\phi_{m,1},...,\phi_{m,p})} and \strong{\eqn{\psi_{m}}}\eqn{=(\psi_{m,1},...,\psi_{m,q_{m}})}.}
#'   \item{For \strong{restricted} models:}{a size \eqn{(pxq)} constraint matrix \strong{\eqn{C}} of full column rank satisfying
#'     \strong{\eqn{\phi}}\eqn{=}\strong{\eqn{C\psi}}, where \strong{\eqn{\phi}}\eqn{=(\phi_{1},...,\phi_{p})} and
#'     \strong{\eqn{\psi}}\eqn{=\psi_{1},...,\psi_{q}}.}
#'   }
#'   Symbol \eqn{\phi} denotes an AR coefficient. Note that regardless of any constraints, the nominal autoregressive order
#'   is always \code{p} for all regimes.
#'   Ignore or set to \code{NULL} if applying linear constraints is \strong{not} desired.
#' @param conditional a logical argument specifying whether the conditional or exact log-likelihood function should be used.
#' @param parametrization is the model parametrized with the "intercepts" \eqn{\phi_{m,0}} or
#'  "means" \eqn{\mu_m = \phi_{m,0}/(1-\sum\phi_{i,m})}?
#' @param boundaries a logical argument. If \code{TRUE}, then \code{loglikelihood} returns \code{minval} if...
#' \itemize{
#'   \item some component variance is not larger than zero,
#'   \item some parametrized mixing weight \eqn{\alpha_{1},...,\alpha_{M-1}} is not larger than zero,
#'   \item sum of the parametrized mixing weights is not smaller than one,
#'   \item if the model is not stationary,
#'   \item or if \code{model=="StMAR"} or \code{model=="G-StMAR"} and some degrees of freedom parameter \eqn{\nu_{m}} is not larger than two.
#' }
#' Argument \code{minval} will be used only if \code{boundaries==TRUE}.
#' @param checks \code{TRUE} or \code{FALSE} specifying whether argument checks, such as stationarity checks, should be done.
#' @param to_return should the returned object be the log-likelihood value, mixing weights, mixing weights including
#'   value for \eqn{alpha_{m,T+1}}, a list containing log-likelihood value and mixing weights, the terms \eqn{l_{t}: t=1,..,T}
#'   in the log-likelihood function (see \emph{KMS 2015, eq.(13)}), regimewise conditional means, regimewise conditional variances,
#'   total conditional means, total conditional variances, or quantile residuals?
#'   Default is the log-likelihood value (\code{"loglik"}).
#' @param minval this will be returned when the parameter vector is outside the parameter space and \code{boundaries==TRUE}.
#' @return
#'  Note that the first p observations are taken as the initial values so the mixing weights and conditional moments start
#'  from the p+1:th observation (interpreted as t=1).
#'  \describe{
#'   \item{By default:}{log-likelihood value of the specified model,}
#'   \item{If \code{to_return=="mw"}:}{a size ((n_obs-p)xM) matrix containing the mixing weights: for m:th component in the m:th column.}
#'   \item{If \code{to_return=="mw_tplus1"}:}{a size ((n_obs-p+1)xM) matrix containing the mixing weights: for m:th component in the m:th column.
#'     The last row is for \eqn{\alpha_{m,T+1}}.}
#'   \item{If \code{to_return=="loglik_and_mw"}:}{a list of two elements. The first element contains the log-likelihood value and the
#'     second element contains the mixing weights.}
#'   \item{If \code{to_return=="terms"}:}{a size ((n_obs-p)x1) numeric vector containing the terms \eqn{l_{t}}.}
#'   \item{If \code{to_return=="regime_cmeans"}:}{a size ((n_obs-p)xM) matrix containing the regime specific conditional means.}
#'   \item{If \code{to_return=="regime_cvars"}:}{a size ((n_obs-p)xM) matrix containing the regime specific conditional variances.}
#'   \item{If \code{to_return=="total_cmeans"}:}{a size ((n_obs-p)x1) vector containing the total conditional means.}
#'   \item{If \code{to_return=="total_cvars"}:}{a size ((n_obs-p)x1) vector containing the total conditional variances.}
#'   \item{If \code{to_return=="qresiduals"}:}{a size ((n_obs-p)x1) vector containing the quantile residuals.}
#'  }
#' @references
#'  \itemize{
#'    \item Galbraith, R., Galbraith, J. 1974. On the inverses of some patterned matrices arising
#'            in the theory of stationary time series. \emph{Journal of Applied Probability} \strong{11}, 63-71.
#'    \item Kalliovirta L. (2012) Misspecification tests based on quantile residuals.
#'            \emph{The Econometrics Journal}, \strong{15}, 358-393.
#'    \item Kalliovirta L., Meitz M. and Saikkonen P. 2015. Gaussian Mixture Autoregressive model for univariate time series.
#'            \emph{Journal of Time Series Analysis}, \strong{36}, 247-266.
#'    \item Meitz M., Preve D., Saikkonen P. 2018. A mixture autoregressive model based on Student's t-distribution.
#'            arXiv:1805.04010 \strong{[econ.EM]}.
#'    \item Virolainen S. 2020. A mixture autoregressive model based on Gaussian and Student's t-distribution.	arXiv:2003.05221 [econ.EM].
#'  }

loglikelihood_int <- function(data, p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL,
                              conditional=TRUE, parametrization=c("intercept", "mean"), boundaries=TRUE, checks=TRUE,
                              to_return=c("loglik", "mw", "mw_tplus1", "loglik_and_mw", "terms", "regime_cmeans", "regime_cvars",
                                          "total_cmeans", "total_cvars", "qresiduals"), minval) {
  epsilon <- round(log(.Machine$double.xmin) + 10)
  model <- match.arg(model)
  parametrization <- match.arg(parametrization)
  to_return <- match.arg(to_return)
  M_orig <- M
  if(model == "G-StMAR") {
    M1 <- M[1]
    M2 <- M[2]
    M <- sum(M)
  } else if(model == "GMAR") {
    M1 <- M
    M2 <- 0
  } else { # model == "StMAR
    M1 <- 0
    M2 <- M
  }

  # Reform parameters to the "standard form" and collect them
  if(checks) check_constraint_mat(p=p, M=M_orig, restricted=restricted, constraints=constraints)
  params <- remove_all_constraints(p=p, M=M_orig, params=params, model=model, restricted=restricted, constraints=constraints)
  pars <- pick_pars(p=p, M=M_orig, params=params, model=model, restricted=FALSE, constraints=NULL)
  alphas <- pick_alphas(p=p, M=M_orig, params=params, model=model, restricted=FALSE, constraints=NULL)
  dfs <- pick_dfs(p=p, M=M_orig, params=params, model=model)
  sigmas <- pars[p + 2,] # sigma^2

  # Return minval if the parameter is outside the parameters space
  if(boundaries) {
    if(any(pars[p + 2,] <= 0)) {
      return(minval)
    } else if(M >= 2 & sum(alphas[-M]) >= 1) {
      return(minval)
    } else if(any(alphas <= 0)) {
      return(minval)
    } else if(!is_stationary_int(p, M, params, restricted=FALSE)) {
      return(minval)
    }
    if(model == "StMAR" | model == "G-StMAR") {
      if(any(dfs <= 2 + 1e-8 | dfs > 1e+5)) return(minval)
    }
  }

  if(checks) {
    data <- check_and_correct_data(data=data, p=p)
    parameter_checks(p=p, M=M_orig, params=params, model=model, restricted=FALSE, constraints=NULL)
  }
  n_obs <- length(data)


  ## Start evaluating the log-likelihood ##

  # Unconditional regimewise means, mu_m (KMS 2015, s.250, and MDP 2018, eq.(4))
  if(parametrization == "mean") {
    mu <- pars[1,]
    pars[1,] <- mu*(1 - colSums(pars[2:(p + 1), , drop=FALSE]))
  } else {
    mu <- pars[1, ]/(1 - colSums(pars[2:(p + 1), , drop=FALSE]))
  }

  # Observed data: y_(-p+1),...,y_0,y_1,...,y_(n_obs-p). First row denotes vector y_0, i:th row vector y_[i-1] and last row denotes the vector y_T.
  Y <- vapply(1:p, function(i1) data[(p - i1 + 1):(n_obs - i1 + 1)], numeric(n_obs - p + 1))

  # Calculate inverse Gamma_m (see the covariance matrix Gamma_p in MPS 2018, p.3 - we calculate this for all mixture components using
  # the inverse formula in Galbraith and Galbraith 1974). Also, calculate the matrix products in multivariate normal and t-distribution
  # densities.
  matProd <- matrix(nrow=n_obs - p + 1, ncol=M)
  invG <- array(dim=c(p, p, M))
  if(p == 1) {
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

  # Calculate the multivariate normal or student's t values (KMS 2015, eq.(7) and MPS 2018, Theorem 1) in log for each vector y_t and for each m=1,..,M.
  # First row for initial values \bm{y}_0 (as denoted by KMS 2015) and i:th row for \bm{y}_(i-1). First column for component m=1 and j:th column for m=j.
  logmv_values <- matrix(nrow=(n_obs - p + 1), ncol=M)
  if(model == "GMAR" | model == "G-StMAR") { # Multinormals
    for(i1 in 1:M1) {
      detG <- 1/det(as.matrix(invG[, , i1]))
      logmv_values[,i1] <- -0.5*p*log(2*base::pi) - 0.5*log(detG) - 0.5*matProd[,i1]
    }
  }
  if(model == "StMAR" | model == "G-StMAR") { # Multistudents
    for(i1 in (M1 + 1):M) {
      detG <- 1/det(as.matrix(invG[, , i1]))
      logC <- lgamma(0.5*(p + dfs[i1 - M1])) - 0.5*p*log(base::pi) - 0.5*p*log(dfs[i1 - M1] - 2) - lgamma(0.5*dfs[i1 - M1])
      logmv_values[,i1] <- logC - 0.5*log(detG) - 0.5*(p + dfs[i1 - M1])*log(1 + matProd[,i1]/(dfs[i1 - M1] - 2))
    }
  }

  # Calculate the mixing weights alpha_mt (KMS 2015, eq.(8) and MPS 2018, eq.(11)).
  # First row for t=1, second for t=2, and i:th for t=i. First column for m=1, second for m=2 and j:th column for m=j.
  if(to_return != "mw_tplus1") {
    logmv_values0 <- logmv_values[1:(n_obs - p),] # The last row is not needed because alpha_mt uses vector Y_(t-1)
  } else {
    logmv_values0 <- logmv_values # The last row is needed for alpha_{m,t+1}
  }
  if(!is.matrix(logmv_values0)) logmv_values0 <- as.matrix(logmv_values0)

  alpha_mt_and_l_0 <- get_alpha_mt(M=M, log_mvnvalues=logmv_values0, alphas=alphas,
                                   epsilon=epsilon, conditional=conditional, to_return=to_return,
                                   also_l_0=TRUE)
  alpha_mt <- alpha_mt_and_l_0$alpha_mt
  l_0 <- alpha_mt_and_l_0$l_0 # The first term in the exact log-likelihood function (=0 for conditional)

  if(to_return == "mw" | to_return == "mw_tplus1") {
    return(alpha_mt)
  }

  # Calculate the conditional means mu_mt (KMS 2015, eq.(2), MPS 2018, eq.(5)). First row for t=1, second for t=2 etc. First column for m=1,
  # second column for m=2 etc.
  mu_mt <- t(pars[1,] + t(Y[-nrow(Y),]%*%pars[2:(p + 1), , drop=FALSE]))

  # Calculate/return conditional means
  if(to_return == "regime_cmeans") {
    return(mu_mt)
  } else if(to_return == "total_cmeans") { # KMS 2015, eq.(4), MPS 2018, eq.(13)
    return(rowSums(alpha_mt*mu_mt))
  }

  # Calculate "the second term" of the log-likelihood (KMS 2015, eq.(12)-(13), MPS 2018, eq.(14)-(15)) or quantile residuals
  Y2 <- Y[2:nrow(Y), 1] # Only the first column and rows 2...T are needed

  # GMAR type components
  if(model == "GMAR" | model == "G-StMAR") {
    invsqrt_sigmasM1 <- sigmas[1:M1]^(-1/2) # M1 = M for the GMAR model
    smat <- diag(x=invsqrt_sigmasM1, nrow=length(invsqrt_sigmasM1), ncol=length(invsqrt_sigmasM1))

    if(to_return == "qresiduals") { # Calculate quantile residuals; see Kalliovirta 2012 for the general formulas and framework
      resM1 <- alpha_mt[,1:M1]*pnorm((Y2 - mu_mt[,1:M1])%*%smat)
      lt_tmpM1 <- resM1 # We exploit the same names
    } else {
      lt_tmpM1 <- alpha_mt[,1:M1]*dnorm((Y2 - mu_mt[,1:M1])%*%smat)%*%smat
    }
  }

  # StMAR type components
  if(model == "StMAR" | model == "G-StMAR") {
    sigmasM2 <- sigmas[(M1 + 1):M] # M1 = 0 and M2 = M for the StMAR model
    matProd0 <- matProd[1:(n_obs - p), (M1 + 1):M] # The last row is not needed because sigma_t uses y_{t-1}
    smat <- diag(x=sigmasM2, nrow=length(sigmasM2), ncol=length(sigmasM2))
    dfmat1 <- diag(x=1/(dfs - 2 + p), nrow=length(dfs), ncol=length(dfs))
    dfmat2 <- diag(x=dfs + p - 2, nrow=length(dfs), ncol=length(dfs))
    sigma_mt <- crossprod(dfs - 2 + t(matProd0), dfmat1)%*%diag(x=sigmasM2, nrow=length(sigmasM2), ncol=length(sigmasM2))

    if(to_return == "qresiduals") { # Calculate the integrals for the quantile residuals
      resM2 <- matrix(ncol=M2, nrow=n_obs - p)

      # Function for numerical integration of the pdf
      my_integral <- function(i1, i2) { # Takes in the regime index i1 and the observation index i2 for the upper bound
        f_mt <- function(y_t) { # The conditional density function to be integrated numerically
          alpha_mt[i2, M1 + i1]*exp(lgamma(0.5*(1 + dfs[i1] + p)) - lgamma(0.5*(dfs[i1] + p)))/sqrt(sigma_mt[i2, i1]*base::pi*(dfs[i1] + p - 2))*
            (1 + ((y_t - mu_mt[i2, M1 + i1])^2)/((dfs[i1] + p - 2)*sigma_mt[i2, i1]))^(-0.5*(1 + dfs[i1] + p))
        }
        tryCatch(integrate(f_mt, lower=-Inf, upper=Y2[i2])$value, # Integrate PDF numerically
                 error=function(e) {
                   warning("Couldn't analytically nor numerically integrate all quantile residuals:")
                   warning(e)
                   return(NA)
                 })
      }

#      is_gsl <- requireNamespace("gsl", quietly = TRUE) # If 'gsl' available, calculate with hypergeometric function what can be calculated
      for(i1 in 1:M2) { # Go through StMAR type regimes
      #  if(is_gsl) {
        whichDef <- which(abs(mu_mt[, M1 + i1] - Y2) < sqrt(sigma_mt[,i1]*(dfs[i1] + p - 2))) # Which ones can be calculated with hypergeometric function
        whichNotDef <- (1:length(Y2))[-whichDef]
       # } #else {
        #  whichDef <- integer(0)
        #  whichNotDef <- 1:length(Y2)
        #}

        if(length(whichDef) > 0) { # Calculate the CDF values at y_t using hypergeometric function whenever it's defined
          Y0 <- Y2[whichDef]
          alpha_mt0 <- alpha_mt[whichDef, M1 + i1]
          mu_mt0 <- mu_mt[whichDef, M1 + i1]
          sigma_mt0 <- sigma_mt[whichDef, i1]
          a0 <- exp(lgamma(0.5*(1 + dfs[i1] + p)) - lgamma(0.5*(dfs[i1] + p)))/sqrt(sigma_mt0*base::pi*(dfs[i1] + p - 2))
          resM2[whichDef, i1] <- alpha_mt0*(0.5 - a0*(mu_mt0 - Y0)*gsl::hyperg_2F1(0.5, 0.5*(1 + dfs[i1] + p), 1.5,
                                                                                   -((mu_mt0 - Y0)^2)/(sigma_mt0*(dfs[i1] + p - 2)),
                                                                                   give=FALSE, strict=TRUE))
        }
        # Calculate the CDF values at y_t that can't be calculated with the hypergeometric function (from the package 'gsl')
        if(length(whichNotDef) > 0) {
          for(i2 in whichNotDef) {
            resM2[i2, i1] <- my_integral(i1, i2)
          }
        }
      }
      lt_tmpM2 <- resM2 # We exploit the same names
    } else { # Calculate l_t in the log-likelihood function
      lt_tmpM2 <- alpha_mt[,(M1 + 1):M]*t(exp(lgamma(0.5*(1 + dfs + p)) - lgamma(0.5*(dfs + p)))/sqrt(base::pi*(dfs + p - 2))/t(sqrt(sigma_mt)))*
                  t(t(1 + ((Y2 - mu_mt[,(M1 + 1):M])^2)/(sigma_mt%*%dfmat2))^(-0.5*(1 + dfs + p)))
    }
  }

  if(model == "GMAR") {
    lt_tmp <- as.matrix(lt_tmpM1)
  } else if(model == "StMAR") {
    lt_tmp <- as.matrix(lt_tmpM2)
  } else { # model == "G-StMAR
    lt_tmp <- as.matrix(cbind(lt_tmpM1, lt_tmpM2))
  }
  l_t <- rowSums(lt_tmp)

  # Return quantile residuals (note that l_t is different if qresiduals are not to be returned)
  if(to_return == "qresiduals") {
    res <- l_t

    # To prevent problems with numerical approximations
    res[which(res >= 1)] <- 1 - 2e-16
    res[which(res <= 0)] <- 2e-16
    return(qnorm(res))
  }

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
    } else { # Calculate and return the total conditional variances (KMS 2015, eq.(5), MPS 2018, eq.(13))
      return(rowSums(alpha_mt*sigma_mt) - rowSums(alpha_mt*(mu_mt - rowSums(alpha_mt*mu_mt))^2))
    }
  }

  if(to_return == "terms") {
    ret <- log(l_t)
  } else if(to_return == "loglik_and_mw") {
    ret <- list(loglik=l_0 + sum(log(l_t)), mw=alpha_mt)
  } else {
    ret <- l_0 + sum(log(l_t)) # KMS 2015, eq.(12)-(13), MPS 2018, eq.(14)-(15)
  }
  ret
}



#' @title Get mixing weights alpha_mt (this function is for internal use)
#'
#' @description \code{get_alpha_mt} computes the mixing weights based on
#'   the logarithm of the multivariate normal densities in the definition of
#'   the mixing weights.
#'
#' @inheritParams loglikelihood_int
#' @param log_mvnvalues \eqn{T x M} matrix containing the log multivariate normal densities.
#' @param alphas \eqn{M x 1} vector containing the mixing weight pa
#' @param epsilon the smallest number such that its exponent is wont classified as numerically zero
#'   (around \code{-698} is used).
#' @param also_l_0 return also l_0 (the first term in the exact log-likelihood function)?
#' @details Note that we index the time series as \eqn{-p+1,...,0,1,...,T} as in Kalliovirta et al. (2015).
#' @return Returns the mixing weights a matrix of the same dimension as \code{log_mvnvalues} so
#'   that the t:th row is for the time point t and m:th column is for the regime m.
#' @inherit loglikelihood_int references
#' @seealso \code{\link{loglikelihood_int}}

get_alpha_mt <- function(M, log_mvnvalues, alphas, epsilon, conditional, to_return, also_l_0=FALSE) {
  if(M == 1) {
    if(!is.matrix(log_mvnvalues)) log_mvnvalues <- as.matrix(log_mvnvalues) # Possibly many time points but only one regime
    alpha_mt <- as.matrix(rep(1, nrow(log_mvnvalues)))
  } else {
    if(!is.matrix(log_mvnvalues)) log_mvnvalues <- t(as.matrix(log_mvnvalues)) # Only one time point but multiple regimes

    log_mvnvalues_orig <- log_mvnvalues
    small_logmvns <- log_mvnvalues < epsilon
    if(any(small_logmvns)) {
      # If too small or large non-log-density values are present (i.e., that would yield -Inf or Inf),
      # we replace them with ones that are not too small or large but imply the same mixing weights
      # up to negligible numerical tolerance.
      which_change <- rowSums(small_logmvns) > 0 # Which rows contain too small  values
      to_change <- log_mvnvalues[which_change, , drop=FALSE]
      largest_vals <- do.call(pmax, split(to_change, f=rep(1:ncol(to_change), each=nrow(to_change)))) # The largest values of those rows
      diff_to_largest <- to_change - largest_vals # Differences to the largest value of the row

      # For each element in each row, check the (negative) distance from the largest value of the row. If the difference
      # is smaller than epsilon, replace the with epsilon. The results are then the new log_mvn values.
      diff_to_largest[diff_to_largest < epsilon] <- epsilon

      # Replace the old log_mvnvalues with the new ones
      log_mvnvalues[which_change,] <- diff_to_largest
    }

    mvnvalues <- exp(log_mvnvalues)
    denominator <- as.vector(mvnvalues%*%alphas)
    alpha_mt <- (mvnvalues/denominator)%*%diag(alphas)
  }

  if(!also_l_0) {
    return(alpha_mt)
  } else {
    # First term of the exact log-likelihood (Kalliovirta et al. 2016, eq.(9))
    l_0 <- 0
    if(M == 1 && conditional == FALSE && (to_return == "loglik" | to_return == "loglik_and_mw")) {
      l_0 <- log_mvnvalues[1]
    } else if(M > 1 && conditional == FALSE && (to_return == "loglik" | to_return == "loglik_and_mw")) {
      if(any(log_mvnvalues_orig[1,] < epsilon)) { # Need to use Brobdingnag
        l_0 <- log(Reduce("+", lapply(1:M, function(i1) alphas[i1]*exp(Brobdingnag::as.brob(log_mvnvalues_orig[1, i1])))))
      } else {
        l_0 <- log(sum(alphas*mvnvalues[1,]))
      }
    }
    return(list(alpha_mt=alpha_mt,
                l_0=l_0))
  }
}


#' @import stats
#'
#' @title Compute the log-likelihood of GMAR, StMAR, or G-StMAR model
#'
#' @description \code{loglikelihood} computes the log-likelihood of the specified GMAR, StMAR, or G-StMAR model.
#'   Exists for convenience if one wants to for example plot profile log-likelihoods or employ other estimation algorithms.
#'   Use \code{minval} to control what happens when the parameter vector is outside the parameter space.
#'
#' @inheritParams loglikelihood_int
#' @param returnTerms should the terms \eqn{l_{t}: t=1,..,T} in the log-likelihood function (see \emph{KMS 2015, eq.(13)}
#'   or MPS 2018, eq.(15)) be returned instead of the log-likelihood value?
#' @return Returns the log-likelihood value or the terms described in \code{returnTerms}.
#' @inherit loglikelihood_int references
#' @seealso \code{\link{fitGSMAR}}, \code{\link{GSMAR}}, \code{\link{quantile_residuals}},
#'  \code{\link{mixing_weights}}, \code{\link{calc_gradient}}
#' @examples
#' # StMAR model
#' params43 <- c(0.09, 1.31, -0.46, 0.33, -0.23, 0.04, 0.01, 1.15,
#'  -0.3, -0.03, 0.03, 1.54, 0.06, 1.19, -0.3, 0.42, -0.4, 0.01,
#'   0.57, 0.22, 8.05, 2.02, 10000)
#' loglikelihood(T10Y1Y, p=4, M=3, params=params43, model="StMAR")
#'
#' # Restricted G-StMAR-model
#' params42gsr <- c(0.11, 0.03, 1.27, -0.39, 0.24, -0.17, 0.03, 1.01, 0.3, 2.03)
#' loglikelihood(T10Y1Y, p=4, M=c(1, 1), params=params42gsr, model="G-StMAR",
#'   restricted=TRUE)
#'
#' # Two-regime GMAR p=2 model with the second AR coeffiecient of
#' # of the second regime contrained to zero.
#' constraints <- list(diag(1, ncol=2, nrow=2), as.matrix(c(1, 0)))
#' params22c <- c(0.03, 1.27, -0.29, 0.03, -0.01, 0.91, 0.34, 0.88)
#' loglikelihood(T10Y1Y, p=2, M=2, params=params22c, model="GMAR",
#'  constraints=constraints)
#' @export

loglikelihood <- function(data, p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL,
                          conditional=TRUE, parametrization=c("intercept", "mean"), returnTerms=FALSE, minval=NA) {
  model <- match.arg(model)
  check_model(model)
  parametrization <- match.arg(parametrization)
  check_pM(p=p, M=M, model=model)
  check_params_length(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints)
  to_ret <- ifelse(returnTerms, "terms", "loglik")
  loglikelihood_int(data=data, p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints,
                    conditional=conditional, parametrization=parametrization, boundaries=TRUE, checks=FALSE,
                    to_return=to_ret, minval=minval)
}


#' @import stats
#'
#' @title Calculate mixing weights of a GMAR, StMAR, or G-StMAR model
#'
#' @description \code{mixing_weights_int} calculates the mixing weights of the specified GMAR, StMAR, or G-StMAR model
#'  and returns them as a matrix.
#'
#' @inheritParams loglikelihood_int
#' @param to_return should the returned object contain mixing weights for t=1,..,T (\code{"mw"}) or for t=1,..,T+1 (\code{"mw_tplus1"})?
#' @details The first p observations are taken to be the initial values.
#' @return
#'  \describe{
#'   \item{If \code{to_return=="mw"}:}{a size ((n_obs-p)xM) matrix containing the mixing weights: for m:th component in m:th column.}
#'   \item{If \code{to_return=="mw_tplus1"}:}{a size ((n_obs-p+1)xM) matrix containing the mixing weights: for m:th component in m:th column.
#'     The last row is for \eqn{\alpha_{m,T+1}}}.
#'  }
#' @inherit loglikelihood_int references

mixing_weights_int <- function(data, p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL,
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
#' @description \code{mixing_weights} calculates the mixing weights of the specified GMAR, StMAR or G-StMAR model and returns them as a matrix.
#'
#' @inheritParams mixing_weights_int
#' @inherit mixing_weights_int return details references
#' @examples
#' # StMAR model
#' params43 <- c(0.09, 1.31, -0.46, 0.33, -0.23, 0.04, 0.01, 1.15,
#'  -0.3, -0.03, 0.03, 1.54, 0.06, 1.19, -0.3, 0.42, -0.4, 0.01,
#'   0.57, 0.22, 8.05, 2.02, 10000)
#' mixing_weights(T10Y1Y, p=4, M=3, params=params43, model="StMAR")
#'
#' # Restricted G-StMAR-model
#' params42gsr <- c(0.11, 0.03, 1.27, -0.39, 0.24, -0.17, 0.03, 1.01, 0.3, 2.03)
#' mixing_weights(T10Y1Y, p=4, M=c(1, 1), params=params42gsr, model="G-StMAR",
#'   restricted=TRUE)
#'
#' # Two-regime GMAR p=2 model with the second AR coeffiecient of
#' # of the second regime contrained to zero.
#' constraints <- list(diag(1, ncol=2, nrow=2), as.matrix(c(1, 0)))
#' params22c <- c(0.03, 1.27, -0.29, 0.03, -0.01, 0.91, 0.34, 0.88)
#' mixing_weights(T10Y1Y, p=2, M=2, params=params22c, model="GMAR",
#'  constraints=constraints)
#' @export

mixing_weights <- function(data, p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL,
                          parametrization=c("intercept", "mean")) {
  model <- match.arg(model)
  check_model(model)
  parametrization <- match.arg(parametrization)
  check_pM(p=p, M=M, model=model)
  check_params_length(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints)
  mixing_weights_int(data=data, p=p, M=M, params=params, model=model, restricted=restricted,
                    constraints=constraints, parametrization=parametrization, checks=TRUE, to_return="mw")
}


#' @import stats
#'
#' @title Calculate conditional moments of GMAR, StMAR, or G-StMAR model
#'
#' @description \code{cond_moments} calculates the regime specific conditional means and variances and total
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
#' params12 <- c(0.01, 0.99, 0.02, 0.03, 0.91, 0.32, 0.86)
#' rcm12 <- cond_moments(T10Y1Y, 1, 2, params12, to_return="regime_cmeans")
#' rcv12 <- cond_moments(T10Y1Y, 1, 2, params12, to_return="regime_cvars")
#' tcm12 <- cond_moments(T10Y1Y, 1, 2, params12, to_return="total_cmeans")
#' tcv12 <- cond_moments(T10Y1Y, 1, 2, params12, to_return="total_cvars")
#'
#' # StMAR model
#' params43 <- c(0.09, 1.31, -0.46, 0.33, -0.23, 0.04, 0.01, 1.15,
#'  -0.3, -0.03, 0.03, 1.54, 0.06, 1.19, -0.3, 0.42, -0.4, 0.01,
#'   0.57, 0.22, 8.05, 2.02, 10000)
#' rcm43t <- cond_moments(T10Y1Y, 4, 3, params43, model="StMAR",
#'  to_return="regime_cmeans")
#' rcv43t <- cond_moments(T10Y1Y, 4, 3, params43, model="StMAR",
#'  to_return="regime_cvars")
#'
#' # G-StMAR model
#' params42gsr <- c(0.11, 0.03, 1.27, -0.39, 0.24, -0.17, 0.03, 1.01, 0.3, 2.03)
#' rcv42gsr <- cond_moments(T10Y1Y, 4, c(1, 1), params42gsr, model="G-StMAR",
#'  restricted=TRUE, to_return="regime_cvars")
#' tcv42gs <- cond_moments(T10Y1Y, 4, c(1, 1), params42gsr, model="G-StMAR",
#'  restricted=TRUE, to_return="total_cvars")
#' @export

cond_moments <- function(data, p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL,
                        parametrization=c("intercept", "mean"), to_return=c("regime_cmeans", "regime_cvars", "total_cmeans", "total_cvars")) {
  model <- match.arg(model)
  check_model(model)
  parametrization <- match.arg(parametrization)
  check_pM(p=p, M=M, model=model)
  check_params_length(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints)
  loglikelihood_int(data=data, p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints,
                    conditional=TRUE, parametrization=parametrization, boundaries=FALSE, checks=TRUE, to_return=to_return)
}


#' @title Calculate AIC, HQIC and BIC
#'
#' @description \code{get_IC} calculates AIC, HQIC and BIC
#'
#' @param loglik log-likelihood value
#' @param npars the number of (freely estimated) parameters in the model.
#' @param obs the number of observations with initial values excluded for conditional models.
#' @return Returns a data frame containing the information criteria values.

get_IC <- function(loglik, npars, obs) {
  AIC <- -2*loglik + 2*npars
  HQIC <- -2*loglik + 2*npars*log(log(obs))
  BIC <- -2*loglik + npars*log(obs)
  data.frame(AIC=AIC, HQIC=HQIC, BIC=BIC)
}

#' @import stats
#'
#' @title Compute quantile residuals of GMAR, StMAR or G-StMAR model
#'
#' @description \code{quantileResiduals_int} computes the quantile residuals of the specified GMAR, StMAR or G-StMAR model.
#'
#' @inheritParams loglikelihood_int
#' @return Returns a \eqn{(Tx1)} numeric vector containing the quantile residuals associated with
#'  the specified GMAR, StMAR or G-StMAR model.
#' @section Suggested packages:
#'   Install the suggested package "gsl" for faster evaluation of the quantile residuals of StMAR and G-StMAR models.
#' @references
#'  \itemize{
#'    \item Galbraith, R., Galbraith, J. 1974. On the inverses of some patterned matrices arising
#'            in the theory of stationary time series. \emph{Journal of Applied Probability} \strong{11}, 63-71.
#'    \item Kalliovirta L. (2012) Misspecification tests based on quantile residuals.
#'            \emph{The Econometrics Journal}, \strong{15}, 358-393.
#'    \item Kalliovirta L., Meitz M. and Saikkonen P. 2015. Gaussian Mixture Autoregressive model for univariate time series.
#'            \emph{Journal of Time Series Analysis}, \strong{36}, 247-266.
#'    \item Lutkepohl H. 2005. New Introduction to Multiple Time Series Analysis.
#'            \emph{Springer}.
#'    \item Meitz M., Preve D., Saikkonen P. 2018. A mixture autoregressive model based on Student's t-distribution.
#'            arXiv:1805.04010 \strong{[econ.EM]}.
#'    \item There are currently no published references for G-StMAR model, but it's a straight forward generalization with
#'            theoretical properties similar to GMAR and StMAR models.
#'  }

quantileResiduals_int <- function(data, p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE,
                                  constraints=NULL, parametrization=c("intercept", "mean")) {
  epsilon <- round(log(.Machine$double.xmin) + 10)
  model <- match.arg(model)
  check_model(model)
  parametrization <- match.arg(parametrization)
  M_orig <- M
  if(model == "G-StMAR") {
    M1 <- M[1]
    M2 <- M[2]
    M <- sum(M)
  }

  # Reform and collect parameters
  # Reform parameters to the "standard form" and collect them
  checkConstraintMat(p=p, M=M_orig, restricted=restricted, constraints=constraints)
  params <- removeAllConstraints(p=p, M=M_orig, params=params, model=model, restricted=restricted, constraints=constraints)
  pars <- pick_pars(p=p, M=M_orig, params=params, model=model, restricted=FALSE, constraints=NULL)
  alphas <- pick_alphas(p=p, M=M_orig, params=params, model=model, restricted=FALSE, constraints=NULL)
  dfs <- pick_dfs(p=p, M=M_orig, params=params, model=model)
  sigmas <- pars[p + 2,]

  data <- checkAndCorrectData(data, p)
  parameterChecks(p=p, M=M_orig, params=params, model=model, restricted=FALSE, constraints=NULL)
  n_obs <- length(data)

  #### Start calculating the quantile residuals ####

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
  matProd <- matrix(nrow=n_obs - p + 1, ncol = M)
  invG <- array(dim=c(p, p, M))
  if(p == 1) { # Lutkepohl (2005), s.15-29
    for(i1 in 1:M) {
      A <- pars[p + 1, i1]
      Sigma <- as.matrix(sigmas[i1])
      VecGamma <- solve(1 - kronecker(A, A), Sigma)
      invG[, , i1] <- as.matrix(1/VecGamma)
      matProd[,i1] <- (Y - mu[i1])*invG[, , i1]*(Y - mu[i1])
    }
  } else { # Galbraith, R., Galbraith, J., (1974)
    for(i1 in 1:M) {
      ARcoefs <- pars[2:(p + 1), i1]
      U <- diag(1, nrow=p, ncol=p)
      V <- diag(ARcoefs[p], nrow=p, ncol=p)
      for(i2 in 1:(p-1)) {
        U[(i2 + 1):p, i2] <- -ARcoefs[1:(p - i2)]
        V[(i2 + 1):p, i2] <- rev(ARcoefs[i2:(p - 1)])
      }
      invG[, , i1] <- (crossprod(U, U) - crossprod(V, V))/sigmas[i1]
      matProd[,i1] <- rowSums((Y - mu[i1]*rep(1, p))%*%invG[, , i1]*(Y - mu[i1]*rep(1, p)))
    }
  }

  # Calculate the log multivariate normal or student's t values (Kalliovirta 2015, s250 eq.(7) or MPS 2007) for each vector y_t and for each m=1,..,M
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
  } else {  # If model == "G-StMAR"
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

  # Calculate the alpha_mt mixing weights (Kalliovirta 2015, s.250 eq.(8)). First row for t=1, second for t=2, and i:th for t=i. First column for m=1, second for m=2 and j:th column for m=j.
  logmv_values0 <- logmv_values[1:(n_obs - p),] # The last row is not needed because alpha_mt uses vector Y_(t-1)
  if(!is.matrix(logmv_values0)) logmv_values0 <- as.matrix(logmv_values0)
  if(M == 1) {
    alpha_mt <- as.matrix(rep(1, n_obs - p))
  } else if(any(logmv_values0 < epsilon)) { # Close to zero values handled with Brobdingnag if needed
    numerators <- lapply(1:M, function(i1) alphas[i1]*exp(Brobdingnag::as.brob(logmv_values0[,i1])))
    denominator <- Reduce("+", numerators) # For all t=0,...,T
    alpha_mt <- vapply(1:M, function(i1) as.numeric(numerators[[i1]]/denominator), numeric(n_obs - p))
  } else {
    mv_values0 <- exp(logmv_values0)
    denominator <- as.vector(mv_values0%*%alphas)
    alpha_mt <- (mv_values0/denominator)%*%diag(alphas)
  }

  # First the mu_mt values (Kalliovirta 2015, s.249 eq.(2)). First row for t=1, second for t=2 etc. First column for m=1, second column for m=2 etc.
  if(p == 1) {
    mu_mt <- vapply(1:M, function(i1) rep(pars[1, i1], nrow(Y) - 1) + Y[1:(nrow(Y) - 1),]*pars[2, i1], numeric(n_obs - p))
  } else {
    mu_mt <- vapply(1:M, function(i1) rep(pars[1, i1], nrow(Y) - 1) + colSums(pars[2:(p + 1), i1]*t(Y[1:(nrow(Y) - 1),])), numeric(n_obs - p) )
  }

  # Calculate the quantile residuals
  Y2 <- Y[2:nrow(Y),1] # Only first column and rows 2...T are needed
  if(model == "GMAR") {
    invsqrt_sigmas <- sigmas^(-1/2)
    res0 <- rowSums(vapply(1:M, function(i1) alpha_mt[,i1]*pnorm((Y2 - mu_mt[,i1])*invsqrt_sigmas[i1]), numeric(n_obs - p)))
  } else if(model == "StMAR") {
    matProd0 <- matProd[1:(n_obs - p),] # Last row is not needed because sigma_t uses y_{t-1}
    if(M == 1) {
      sigma_mt <- as.matrix(sigmas*(dfs - 2 + matProd0)/(dfs - 2 + p))
    } else {
      sigma_mt <- t(dfs - 2 + t(matProd0))%*%diag(1/(dfs - 2 + p))%*%diag(sigmas) # Calculate conditional variances
    }
    res0 <- matrix(ncol=M, nrow=n_obs - p)
    if(requireNamespace("gsl", quietly = TRUE)) { # Calculate with hypergeometric function what can be calculated
      for(i1 in 1:M) {
        whichDef <- which(abs(mu_mt[,i1] - Y2) < sqrt(sigma_mt[,i1]*(dfs[i1] + p - 2))) # Which ones can be calculated with hyper geometric function
        whichNotDef <- (1:length(Y2))[-whichDef]

        if(length(whichDef) > 0) { # Calculate CDF for values y_t that are bigger than mu_mt (or equal)
          Y0 = Y2[whichDef]
          alpha_mt0 <- alpha_mt[whichDef, i1]
          mu_mt0 <- mu_mt[whichDef, i1]
          sigma_mt0 <- sigma_mt[whichDef, i1]
          a0 <- exp(lgamma(0.5*(1 + dfs[i1] + p)) - lgamma(0.5*(dfs[i1] + p)))/sqrt(sigma_mt0*pi*(dfs[i1] + p - 2))
          res0[whichDef, i1] <- alpha_mt0*(0.5 - a0*(mu_mt0 - Y0)*gsl::hyperg_2F1(0.5, 0.5*(1 + dfs[i1] + p), 1.5, -((mu_mt0 - Y0)^2)/(sigma_mt0*(dfs[i1] + p - 2)), give=FALSE, strict=TRUE))
        }
        # Calculate CDF for the values that can't be calculated with hypergeometric function
        if(length(whichNotDef) > 0) {
          for(i2 in whichNotDef) {
            # Conditional density function to integrate numerically
            f_mt <- function(y_t) {
              alpha_mt[i2, i1]*exp(lgamma(0.5*(1 + dfs[i1] + p)) - lgamma(0.5*(dfs[i1] + p)))/sqrt(sigma_mt[i2, i1]*pi*(dfs[i1] + p - 2))*
                (1 + ((y_t-mu_mt[i2, i1])^2)/((dfs[i1] + p - 2)*sigma_mt[i2, i1]))^(-0.5*(1 + dfs[i1] + p))
            }
            res0[i2, i1] <- integrate(f_mt, lower=-Inf, upper=Y2[i2])$value # Integrate PDF numerically
          }
        }
      }
    } else { # Numerically integrate everything if package "gsl" is not available: slow but works always with all platforms
      for(i1 in 1:M) {
        for(i2 in 1:length(Y2)) {
          # Conditional density function to integrate numerically
          f_mt <- function(y_t) {
            alpha_mt[i2, i1]*exp(lgamma(0.5*(1 + dfs[i1] + p)) - lgamma(0.5*(dfs[i1] + p)))/sqrt(sigma_mt[i2, i1]*pi*(dfs[i1] + p - 2))*
              (1 + ((y_t - mu_mt[i2, i1])^2)/((dfs[i1] + p - 2)*sigma_mt[i2, i1]))^(-0.5*(1 + dfs[i1] + p))
          }
          res0[i2, i1] <- integrate(f_mt, lower=-Inf, upper=Y2[i2])$value # Integrate PDF numerically
        }
      }
    }
    res0 <- rowSums(res0)
  } else { # If model == "G-StMAR"
    # GMAR components
    invsqrt_sigmas <- sigmas[1:M1]^(-1/2)
    resM1 <- rowSums(vapply(1:M1, function(i1) alpha_mt[,i1]*pnorm((Y2 - mu_mt[,i1])*invsqrt_sigmas[i1]), numeric(n_obs - p)))

    # StMAR components
    matProd0 <- matProd[1:(n_obs - p),(M1 + 1):M] # Last row is not needed because sigma_t uses y_{t-1}
    sigmasM2 <- sigmas[(M1 + 1):M]
    if(M2==1) {
      sigma_mt <- as.matrix(sigmasM2*(dfs - 2 + matProd0)/(dfs - 2 + p))
    } else {
      sigma_mt <- t(dfs - 2 + t(matProd0))%*%diag(1/(dfs - 2 + p))%*%diag(sigmasM2) # Calculate conditional variances
    }
    resM2 <- matrix(ncol=M2, nrow=n_obs - p)
    if(requireNamespace("gsl", quietly=TRUE)) { # Calculate with hypergeometric function what can be calculated
      for(i1 in 1:M2) {
        whichDef <- which(abs(mu_mt[,M1 + i1] - Y2) < sqrt(sigma_mt[,i1]*(dfs[i1] + p - 2))) # Which ones can be calculated with hyper geometric function
        whichNotDef <- (1:length(Y2))[-whichDef]

        if(length(whichDef) > 0) { # Calculate CDF for values y_t that are bigger than mu_mt (or equal)
          Y0 <- Y2[whichDef]
          alpha_mt0 <- alpha_mt[whichDef, M1 + i1]
          mu_mt0 <- mu_mt[whichDef, M1 + i1]
          sigma_mt0 <- sigma_mt[whichDef, i1]
          a0 <- exp(lgamma(0.5*(1 + dfs[i1] + p)) - lgamma(0.5*(dfs[i1] + p)))/sqrt(sigma_mt0*pi*(dfs[i1] + p - 2))
          resM2[whichDef, i1] <- alpha_mt0*(0.5 - a0*(mu_mt0 - Y0)*gsl::hyperg_2F1(0.5, 0.5*(1 + dfs[i1] + p), 1.5, -((mu_mt0 - Y0)^2)/(sigma_mt0*(dfs[i1] + p - 2)), give=FALSE, strict=TRUE))
        }
        # Calculate CDF for the values that can't be calculated with hypergeometric function
        if(length(whichNotDef) > 0) {
          for(i2 in whichNotDef) {
            # Conditional density function to integrate numerically
            f_mt <- function(y_t) {
              alpha_mt[i2, M1 + i1]*exp(lgamma(0.5*(1 + dfs[i1] + p)) - lgamma(0.5*(dfs[i1] + p)))/sqrt(sigma_mt[i2, i1]*pi*(dfs[i1] + p - 2))*
                (1 + ((y_t - mu_mt[i2, M1 + i1])^2)/((dfs[i1] + p - 2)*sigma_mt[i2, i1]))^(-0.5*(1 + dfs[i1] + p))
            }
            resM2[i2, i1] <- integrate(f_mt, lower=-Inf, upper=Y2[i2])$value # Integrate PDF numerically
          }
        }
      }
    } else { # Numerically integrate everything if package "gsl" is not available
      for(i1 in 1:M2) {
        for(i2 in 1:length(Y2)) {
          # Conditional density function to integrate numerically
          f_mt <- function(y_t) {
            alpha_mt[i2, M1 + i1]*exp(lgamma(0.5*(1 + dfs[i1] + p)) - lgamma(0.5*(dfs[i1] + p)))/sqrt(sigma_mt[i2, i1]*pi*(dfs[i1] + p - 2))*
              (1 + ((y_t - mu_mt[i2, M1 + i1])^2)/((dfs[i1] + p - 2)*sigma_mt[i2, i1]))^(-0.5*(1 + dfs[i1] + p))
          }
          resM2[i2, i1] <- integrate(f_mt, lower=-Inf, upper=Y2[i2])$value # Integrate PDF numerically
        }
      }
    }
    res0 = rowSums(cbind(resM1, resM2))
  }
  # To prevent problems with numerical approximations
  res0[which(res0 >= 1)] <- 1 - 2e-16
  res0[which(res0 <= 0)] <- 2e-16
  qnorm(res0)
}


#' @import stats
#'
#' @title Compute quantile residuals of GMAR, StMAR or G-StMAR model
#'
#' @description \code{quantileResiduals} computes the quantile residuals of the specified GMAR, StMAR or G-StMAR model.
#'
#' @inheritParams quantileResiduals_int
#' @inherit quantileResiduals_int return references
#' @section Suggested packages:
#'   Install the suggested package "gsl" for faster evaluation of the quantile residuals of StMAR and G-StMAR models.
#' @examples
#' # GMAR model
#' params12 <- c(1.12, 0.9, 0.29, 4.53, 0.7, 3.22, 0.84)
#' qr12 <- quantileResiduals(VIX, 1, 2, params12)
#'
#' # Restricted GMAR model
#' params12r <- c(1.4, 1.8, 0.88, 0.29, 3.18, 0.84)
#' qr12r <- quantileResiduals(VIX, 1, 2, params12r, restricted=TRUE)
#'
#' # StMAR model
#' params12t <- c(1.1, 0.9, 0.3, 4.5, 0.7, 3.2, 0.8, 5, 8)
#' qr12t <- quantileResiduals(VIX, 1, 2, params12t, model="StMAR")
#'
#' # Non-mixture version of StMAR model
#' params11t <- c(0.76, 0.93, 1.35, 2.4)
#' qr11t <- quantileResiduals(VIX, 1, 1, params11t, model="StMAR")
#'
#' # G-StMAR model
#' params12gs <- c(1.5, 0.8, 1.5, 2.9, 0.8, 1.1, 0.6, 3)
#' qr12gs <- quantileResiduals(VIX, 1, c(1, 1), params12gs, model="G-StMAR")
#'
#' # Restricted G-StMAR model
#' params13gsr <- c(1.3, 1, 1.4, 0.8, 0.4, 2, 0.2, 0.25, 0.15, 20)
#' qr12gsr <- quantileResiduals(VIX, 1, c(2, 1), params13gsr, model="G-StMAR", restricted=TRUE)
#'
#' # GMAR model as a mixture of AR(2) and AR(1) models
#' constraints <- list(diag(1, ncol=2, nrow=2), as.matrix(c(1, 0)))
#' params22c <- c(1.2, 0.85, 0.04, 0.3, 3.3, 0.77, 2.8, 0.77)
#' qr22c <- quantileResiduals(VIX, 2, 2, params22c, constraints=constraints)
#'
#' # Such StMAR(3,2) that the AR coefficients are restricted to be
#' # the same for both regimes and that the second AR coefficients are
#' # constrained to zero.
#' params32trc <- c(2.2, 1.8, 0.88, -0.03, 2.4, 0.27, 0.40, 3.9, 1000)
#' qr32trc <- quantileResiduals(VIX, 3, 2, params32trc, model="StMAR",
#'  restricted=TRUE, constraints=matrix(c(1, 0, 0, 0, 0, 1), ncol=2))
#' @export

quantileResiduals <- function(data, p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE,
                              constraints=NULL, parametrization=c("intercept", "mean")) {
  model <- match.arg(model)
  check_model(model)
  parametrization <- match.arg(parametrization)
  stopifnot(parametrization %in% c("intercept", "mean"))
  checkPM(p, M, model=model)
  check_params_length(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints)
  quantileResiduals_int(data=data, p=p, M=M, params=params, model=model, restricted=restricted,
                        constraints=constraints, parametrization=parametrization)
}

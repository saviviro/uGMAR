#' @import stats
#'
#' @title Generate covariance matrix Omega for quantile residual tests
#'
#' @description \code{getOmega} generates the covariance matrix Omega used in the quantile residual tests.
#'
#' @inheritParams loglikelihood_int
#' @param g a function specifying the transformation.
#' @param dim_g output dimension of the transformation \code{g}.
#' @details This function is used for quantile residuals tests in \code{quantileResidualTests}.
#' @return Returns size (\code{dim_g}x\code{dim_g}) covariance matrix Omega.
#' @inherit quantileResiduals_int references

getOmega <- function(data, p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL,
                     parametrization=c("intercept", "mean"), g, dim_g) {

  # Function for calculating gradient of g
  f <- function(params) {
    g(quantileResiduals_int(data=data, p=p, M=M, params=params, model=model, restricted=restricted,
                            constraints=constraints, parametrization=parametrization))
  }
  # Function for calculating gradient of the log-likelihood
  l <- function(params) {
    loglikelihood_int(data=data, p=p, M=M, params=params, model=model, restricted=restricted,
                      constraints=constraints, boundaries=FALSE, conditional=FALSE,
                      parametrization=parametrization, checks=FALSE, to_return="terms")
  }

  diff <- 6e-06 # Difference for numerical derivates
  npars <- length(params) # Dimension of the parameter vector

  # Compute the gradient of g: (v)x(npars)x(T)
  gres <- f(params)
  T0 <- nrow(gres)

  I <- diag(rep(1, npars))
  dg <- array(dim=c(dim_g, npars, T0))  # row per g_i, column per derivative and slice per t=1,..,T.
  for(i1 in 1:npars) {
    dg[,i1,] <- t((f(params + I[i1,]*diff) - f(params - I[i1,]*diff))/(2*diff))
  }

  # Compute gradient of the log-likelihood: (T)x(npars)
  dl <- vapply(1:npars, function(i1) (l(params + I[,i1]*diff) - l(params - I[,i1]*diff))/(2*diff), numeric(length(data) - p)) # NOTE: "returnTerms" in loglik is TRUE

  # Estimate Fisher's information matrix
  FisInf <- crossprod(dl, dl)/nrow(dl)
  invFisInf <- solve(FisInf) # Can cause error which needs to be handed in the call function

  # Calculate G (Kalliovirta 2012 eq.(2.4))
  G <- rowMeans(dg, dims=2)

  # Calculate PSI
  dif0 <- nrow(dl) - T0
  PSI <- crossprod(gres, dl[(dif0+1):nrow(dl),])/T0

  # Calculate H
  H <- crossprod(gres, gres)/T0

  # Calculate Omega
  Omega <- G%*%tcrossprod(invFisInf, G) + PSI%*%tcrossprod(invFisInf, G) + G%*%tcrossprod(invFisInf, PSI) + H
  Omega
}

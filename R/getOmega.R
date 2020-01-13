#' @import stats
#'
#' @title Generate the covariance matrix Omega for quantile residual tests
#'
#' @description \code{getOmega} generates the covariance matrix Omega used in the quantile residual tests.
#'
#' @inheritParams loglikelihood_int
#' @param g a function specifying the transformation.
#' @param dim_g output dimension of the transformation \code{g}.
#' @details This function is used for quantile residuals tests in \code{quantileResidualTests}.
#' @seealso \code{\link{quantileResidualTests}}
#' @return Returns size (\code{dim_g}x\code{dim_g}) covariance matrix Omega.
#' @inherit quantileResidualTests references

getOmega <- function(data, p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE, constraints=NULL,
                     parametrization=c("intercept", "mean"), g, dim_g) {
  model <- match.arg(model)
  parametrization <- match.arg(parametrization)

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

  npars <- length(params) # Dimension of the parameter vector
  h <- get_varying_h(p=p, M=M, params=params, model=model) # Differences for numerical derivates: adjust for overly large dfs parameters to avoid numerical problems

  # Compute the gradient of g: (v)x(npars)x(T)
  gres <- f(params)
  T0 <- nrow(gres)

  I <- diag(rep(1, npars))
  dg <- array(dim=c(dim_g, npars, T0))  # Row per g_i, column per derivative, and slice per t=1,..,T.
  for(i1 in 1:npars) {
    dg[,i1,] <- t((f(params + I[i1,]*h[i1]) - f(params - I[i1,]*h[i1]))/(2*h[i1]))
  }

  # Compute gradient of the log-likelihood: (T)x(npars)
  dl <- vapply(1:npars, function(i1) (l(params + I[,i1]*h[i1]) - l(params - I[,i1]*h[i1]))/(2*h[i1]), numeric(length(data) - p)) # NOTE: "returnTerms" in loglik in 'l' is TRUE

  # Estimate Fisher's information matrix
  FisInf <- crossprod(dl, dl)/nrow(dl)
  invFisInf <- solve(FisInf) # Can cause error which needs to be handled in the function which calls getOmega

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

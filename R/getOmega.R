#' @import stats
#'
#' @title Generate covariance matrix omega for quantile residual tests
#'
#' @description \code{getOmega} generates the covariance matrix omega for quantile residual tests.
#'
#' @inheritParams loglikelihood
#' @param g a function specifying the transformation.
#' @param dim_g output dimension of the transformation \code{g}.
#' @param qresiduals optionally provide the quantile residuals of the model to save time when calculating multiple omegas for the same model.
#' @details This function is used for quantile residuals tests in \code{quantileResidualTests}.
#' @return Returns size (\code{dim_g}x\code{dim_g}) covariance matrix Omega.
#' @references
#'  \itemize{
#'    \item Kalliovirta L. (2012) Misspecification tests based on quantile residuals.
#'          \emph{The Econometrics Journal}, \strong{15}, 358-393.
#'  }

getOmega <- function(data, p, M, params, StMAR=FALSE, restricted=FALSE, constraints=FALSE, R, g, dim_g, qresiduals) {

  # Function for calculating gradient of g
  f <- function(params) {
    g(quantileResiduals_int(data, p, M, params, StMAR=StMAR, restricted=restricted, constraints=constraints, R=R))
  }
  # Function for calculating gradient of the log-likelihood
  l <- function(params) {
    loglikelihood_int(data, p, M, params, StMAR=StMAR, restricted=restricted, constraints=constraints, R=R, boundaries=FALSE,
                  conditional=FALSE, checks=FALSE, returnTerms=TRUE)
  }

  diff = 1e-08 # Difference for numerical derivates
  d = length(params) # Dimension of the parameter vector
  dataSize = length(data) # Size of the given data
  if(missing(qresiduals)) {
    qresiduals = quantileResiduals_int(data, p, M, params, StMAR=StMAR, restricted=restricted, constraints=constraints, R=R) # Quantile residuals
  }

  # Compute the gradient of g: (v)x(d)x(T)
  T0 = nrow(f(params))

  I = diag(rep(1, d))
  dg = array(dim=c(dim_g, d, T0))  # row per g_i, column per derivative and slice per t=1,..,T.
  for(i1 in 1:d) {
    dg[,i1,] = t((f(params+I[i1,]*diff)-f(params-I[i1,]*diff))/(2*diff))
  }

  # Compute gradient of the log-likelihood: (T)x(d)
  dl = sapply(1:d, function(i1) (l(params+I[,i1]*diff)-l(params-I[,i1]*diff))/(2*diff)) # NOTE: "returnTerms" in loglik is TRUE

  # Estimate Fisher's information matrix
  FisInf <- matrix(0, ncol=d, nrow=d)
  for(i1 in 1:T0) {
    FisInf = FisInf + tcrossprod(dl[i1,], dl[i1,])
  }
  FisInf = FisInf/(dataSize-p)

  invFisInf <- tryCatch(solve(FisInf), error=function(cond) {
    message("The covariance matrix Omega cannot be solved")
    return(NA)
  })
  if(is.na(invFisInf[1])) { # If cannot be solved
    return(matrix(ncol=dim_g, nrow=dim_g))
  }

  # Calculate G (Kalliovirta 2012 eq.(2.4))
  G = apply(dg, c(1,2), mean)

  # Calculate PSI
  gres = t(g(qresiduals))
  dif0 = nrow(dl)-T0
  PSI = gres%*%dl[(dif0+1):nrow(dl),]/T0

  # Calculate H
  H = gres%*%t(gres)/T0

  # Calculate Omega
  Omega = G%*%invFisInf%*%t(G) + PSI%*%invFisInf%*%t(G) + G%*%invFisInf%*%t(PSI) + H
  return(Omega)
}

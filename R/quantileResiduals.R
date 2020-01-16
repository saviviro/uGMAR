#' @import stats
#'
#' @title Compute quantile residuals of GMAR, StMAR, or G-StMAR model
#'
#' @description \code{quantileResiduals_int} computes the quantile residuals of the specified GMAR, StMAR, or G-StMAR model.
#'
#' @inheritParams loglikelihood_int
#' @return Returns a \eqn{(Tx1)} numeric vector containing the quantile residuals of the specified GMAR, StMAR or G-StMAR model.
#' @details Numerical integration is employed if the quantile residuals cannot be obtained analytically with the
#'  hypergeometric function using the package 'gsl'.
#' @section Suggested packages:
#'   Install the suggested package "gsl" for faster evaluation of the quantile residuals for the StMAR and G-StMAR models.
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
#'    \item There are currently no published references for G-StMAR model, but it's a straight forward generalization with
#'            theoretical properties similar to GMAR and StMAR models.
#'  }

quantileResiduals_int <- function(data, p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE,
                                  constraints=NULL, parametrization=c("intercept", "mean")) {
  loglikelihood_int(data=data, p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints,
                    parametrization=parametrization, checks=TRUE, boundaries=FALSE, to_return="qresiduals", minval=NA)
}


#' @import stats
#'
#' @title Compute quantile residuals of GMAR, StMAR, or G-StMAR model
#'
#' @description \code{quantileResiduals} computes the quantile residuals of the specified GMAR, StMAR, or G-StMAR model.
#'
#' @inheritParams quantileResiduals_int
#' @inherit quantileResiduals_int return details references
#' @section Suggested packages:
#'   Install the suggested package "gsl" for faster evaluation of the quantile residuals of StMAR and G-StMAR models.
#' @examples
#' # GMAR model
#' params12 <- c(0.18, 0.93, 0.01, 0.86, 0.68, 0.02, 0.88)
#' quantileResiduals(logVIX, 1, 2, params12)
#'
#' # Restricted GMAR model, outside parameter space
#' params12r <- c(0.21, 0.23, 0.92, 0.01, 0.02, 0.86)
#' quantileResiduals(logVIX, 1, 2, params12r, restricted=TRUE)
#'
#' # Non-mixture version of StMAR model, outside parameter space
#' params11t <- c(0.16, 0.93, 0.01, 3.01)
#' quantileResiduals(logVIX, 1, 1, params11t, model="StMAR")
#'
#' # G-StMAR model
#' params12gs <- c(0.86, 0.68, 0.02, 0.18, 0.93, 0.01, 0.11, 44.36)
#' quantileResiduals(logVIX, 1, c(1, 1), params12gs, model="G-StMAR")
#'
#' # Restricted G-StMAR model
#' params12gsr <- c(0.31, 0.33, 0.88, 0.01, 0.02, 0.77, 2.72)
#' quantileResiduals(logVIX, 1, c(1, 1), params12gsr, model="G-StMAR",
#'  restricted=TRUE)
#'
#' # GMAR model as a mixture of AR(2) and AR(1) models
#' constraints <- list(diag(1, ncol=2, nrow=2), as.matrix(c(1, 0)))
#' params22c <- c(0.61, 0.83, -0.06, 0.02, 0.21, 0.91, 0.01, 0.16)
#' quantileResiduals(logVIX, 2, 2, params22c, constraints=constraints)
#'
#' # Such StMAR(3,2) that the AR coefficients are restricted to be
#' # the same for both regimes and that the second AR coefficients are
#' # constrained to zero.
#' params32trc <- c(0.35, 0.33, 0.88, -0.02, 0.01, 0.01, 0.36, 4.53, 1000)
#' quantileResiduals(logVIX, 3, 2, params32trc, model="StMAR",
#'  restricted=TRUE, constraints=matrix(c(1, 0, 0, 0, 0, 1), ncol=2))
#' @export

quantileResiduals <- function(data, p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE,
                              constraints=NULL, parametrization=c("intercept", "mean")) {
  model <- match.arg(model)
  check_model(model)
  parametrization <- match.arg(parametrization)
  checkPM(p, M, model=model)
  check_params_length(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints)
  loglikelihood_int(data=data, p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints,
                    parametrization=parametrization, checks=TRUE, boundaries=TRUE, to_return="qresiduals", minval=NA)
}

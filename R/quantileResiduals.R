#' @import stats
#'
#' @title Compute quantile residuals of GMAR, StMAR, or G-StMAR model
#'
#' @description \code{quantile_residuals_int} computes the quantile residuals of the specified GMAR, StMAR, or G-StMAR model.
#'
#' @inheritParams loglikelihood_int
#' @return Returns a \eqn{(Tx1)} numeric vector containing the quantile residuals of the specified GMAR, StMAR or G-StMAR model.
#'  Note that there are no quantile residuals for the first \code{p} observations as they are the initial values.
#' @details Numerical integration is employed if the quantile residuals cannot be obtained analytically with the
#'  hypergeometric function using the package 'gsl'.
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

quantile_residuals_int <- function(data, p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE,
                                  constraints=NULL, parametrization=c("intercept", "mean")) {
  loglikelihood_int(data=data, p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints,
                    parametrization=parametrization, checks=TRUE, boundaries=FALSE, to_return="qresiduals", minval=NA)
}


#' @import stats
#'
#' @title Compute quantile residuals of GMAR, StMAR, or G-StMAR model
#'
#' @description \code{quantile_residuals} computes the quantile residuals of the specified GMAR, StMAR, or G-StMAR model.
#'
#' @inheritParams quantile_residuals_int
#' @inherit quantile_residuals_int return details references
#' @examples
#' # GMAR model
#' params12 <- c(1.70, 0.85, 0.30, 4.12, 0.73, 1.98, 0.63)
#' quantile_residuals(simudata, p=1, M=2, params=params12, model="GMAR")
#'
#' # G-StMAR-model
#' params42gs <- c(0.04, 1.34, -0.59, 0.54, -0.36, 0.01, 0.06, 1.28, -0.36,
#'                 0.2, -0.15, 0.04, 0.19, 9.75)
#' quantile_residuals(M10Y1Y, p=4, M=c(1, 1), params=params42gs, model="G-StMAR")
#' @export

quantile_residuals <- function(data, p, M, params, model=c("GMAR", "StMAR", "G-StMAR"), restricted=FALSE,
                              constraints=NULL, parametrization=c("intercept", "mean")) {
  model <- match.arg(model)
  check_model(model)
  parametrization <- match.arg(parametrization)
  check_pM(p, M, model=model)
  check_params_length(p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints)
  loglikelihood_int(data=data, p=p, M=M, params=params, model=model, restricted=restricted, constraints=constraints,
                    parametrization=parametrization, checks=TRUE, boundaries=TRUE, to_return="qresiduals", minval=NA)
}

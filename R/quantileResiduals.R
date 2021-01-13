#' @import stats
#'
#' @title Compute quantile residuals of GMAR, StMAR, or G-StMAR model
#'
#' @description \code{quantile_residuals_int} computes the quantile residuals of the specified GMAR, StMAR, or G-StMAR model.
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
#' @section Suggested packages:
#'   Install the suggested package "gsl" for faster evaluation of the quantile residuals of StMAR and G-StMAR models.
#' @examples
#' # StMAR model
#' params43 <- c(0.09, 1.31, -0.46, 0.33, -0.23, 0.04, 0.01, 1.15,
#'  -0.3, -0.03, 0.03, 1.54, 0.06, 1.19, -0.3, 0.42, -0.4, 0.01,
#'   0.57, 0.22, 8.05, 2.02, 10000)
#' quantile_residuals(T10Y1Y, p=4, M=3, params=params43, model="StMAR")
#'
#' # Restricted G-StMAR-model
#' params42gsr <- c(0.11, 0.03, 1.27, -0.39, 0.24, -0.17, 0.03, 1.01, 0.3, 2.03)
#' quantile_residuals(T10Y1Y, p=4, M=c(1, 1), params=params42gsr, model="G-StMAR",
#'   restricted=TRUE)
#'
#' # Two-regime GMAR p=2 model with the second AR coeffiecient of
#' # of the second regime contrained to zero.
#' constraints <- list(diag(1, ncol=2, nrow=2), as.matrix(c(1, 0)))
#' params22c <- c(0.03, 1.27, -0.29, 0.03, -0.01, 0.91, 0.34, 0.88)
#' quantile_residuals(T10Y1Y, p=2, M=2, params=params22c, model="GMAR",
#'  constraints=constraints)
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

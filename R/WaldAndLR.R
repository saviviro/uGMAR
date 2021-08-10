#' @title Perform Wald test
#'
#' @description \code{Wald_test} performs a Wald test for a GMAR, StMAR, or G-StMAR model.
#'
#' @inheritParams simulate.gsmar
#' @inheritParams calc_gradient
#' @param A a size \eqn{(k x n_params)} matrix with full row rank specifying a part of the null hypothesis,
#'   where \eqn{n_params} is the number of parameters in the (unconstrained) model.
#'   See details for more information.
#' @param c a length \eqn{k} vector specifying a part of the null hypothesis. See details for more information.
#' @details Denoting the true parameter value by \eqn{\theta_{0}}, we test the null hypothesis \eqn{A\theta_{0}=c}.
#'   Under the null, the test statistic is asymptotically \eqn{\chi^2}-distributed with \eqn{k}
#'   (\code{=nrow(A)}) degrees of freedom. The parameter \eqn{\theta_{0}} is assumed to have the same form as in
#'   the model supplied in the argument \code{gsmar} and it is presented in the documentation of the argument
#'   \code{params} in the function \code{GSMAR} (see \code{?GSMAR}).
#'
#'   Note that this function does \strong{not} check whether the specified constraints are feasible (e.g., whether
#'   the implied constrained model would be stationary or have positive definite error term covariance matrices).
#' @return A list with class "htest" containing the following components:
#'   \item{statistic}{the value of the Wald statistics.}
#'   \item{parameter}{the degrees of freedom of the Wald statistic.}
#'   \item{p.value}{the p-value of the test.}
#'   \item{alternative}{a character string describing the alternative hypothesis.}
#'   \item{method}{a character string indicating the type of the test (Wald test).}
#'   \item{data.name}{a character string giving the names of the supplied model, constraint matrix A, and vector c.}
#'   \item{gsmar}{the supplied argument gsmar.}
#'   \item{A}{the supplied argument A.}
#'   \item{c}{the supplied argument c.}
#'   \item{h}{the supplied argument h.}
#' @seealso \code{\link{LR_test}}, \code{\link{fitGSMAR}}, \code{\link{GSMAR}}, \code{\link{diagnostic_plot}},
#'  \code{\link{profile_logliks}}, \code{\link{quantile_residual_tests}}, \code{\link{cond_moment_plot}}
#' @inherit is_stationary references
#' @examples
#' \donttest{
#' # GMAR p=1, M=2 model:
#' fit12 <- fitGSMAR(simudata, p=1, M=2, model="GMAR", ncalls=1, seeds=1)
#'
#' # Test with Wald test whether the AR coefficients are the same in both
#' # regimes:
#' # There are 7 parameters in the model and the AR coefficient of the
#' # first regime is the 2nd element, whereas the AR coefficient of the second
#' # regime is in the 5th element.
#' A <- matrix(c(0, 1, 0, 0, -1, 0, 0), nrow=1, ncol=7)
#' c <- 0
#' Wald_test(fit12, A=A, c=c)
#' }
#' @export

Wald_test <- function(gsmar, A, c, h=6e-6) {
  # Checks
  check_gsmar(gsmar)
  params <- gsmar$params
  stopifnot(is.matrix(A) && ncol(A) == length(params) && nrow(A) <= ncol(A))
  stopifnot(length(c) == nrow(A))
  stopifnot(!is.null(gsmar$data))
  if(qr(A)$rank != nrow(A)) stop("The constraint matrix 'A' should have full row rank")

  # Calculate Hessian matrix at the estimate
  minval <- get_minval(gsmar$data)
  loglik_fn <- function(pars) {
    tryCatch(loglikelihood_int(data=gsmar$data, p=gsmar$model$p, M=gsmar$model$M, params=pars, model=gsmar$model$model,
                               conditional=gsmar$model$conditional, parametrization=gsmar$model$parametrization,
                               restricted=gsmar$model$restricted, constraints=gsmar$model$constraints, boundaries=TRUE,
                               to_return="loglik", minval=minval),
             error=function(e) {
               print(paste("Failed to evualuate log-likelihood function in the approximation of Hessian matrix:", e))
               return(NA)
             })
  }
  Hess <- calc_hessian(x=params, fn=loglik_fn, h=h)
  if(anyNA(Hess)) stop("Unable to fully calculate Hessian matrix of the log-likelihood function using central difference numerical approximation. Check whether there is something funny in the estimates or maybe try another difference 'h'?")

  # Invert the Hessian matrix
  inv_Hess <- tryCatch(solve(Hess), error=function(e) {
    print(paste("Failed to invert Hessian matrix:", e))
    return(NA)
  })
  if(anyNA(inv_Hess)) stop("Couldn't invert Hessian matrix of the log-likelihood function. This might happen when the mixing weights are very close to zero for some regime (if so, reduce the redundant regime from the model).")

  # Calculate the Wald test statistic
  test_stat <- as.numeric(crossprod(A%*%params - c, solve(-tcrossprod(A%*%inv_Hess, A), A%*%params - c))) # t(A%*%params - c)%*%solve(-A%*%inv_Hess%*%t(A))%*%(A%*%params - c)

  # Calculate the p-value
  df <- nrow(A)
  p_value <- pchisq(test_stat, df=df, lower.tail=FALSE)

  dname <- paste0(deparse(substitute(gsmar)),", ", deparse(substitute(A)), ", ", deparse(substitute(c)))


  structure(list(statistic=c("W"=test_stat),
                 parameter=c("df"=df),
                 p.value=p_value,
                 alternative="the true parameter theta does not satisfy A%*%theta = c",
                 data.name=dname,
                 method="Wald test",
                 gsmar=gsmar,
                 A=A,
                 c=c,
                 h=h),
            class="htest")
}


#' @title Perform likelihood ratio test
#'
#' @description \code{LR_test} performs a likelihood ratio test for a GMAR, StMAR, or G-StMAR model.
#'
#' @param gsmar1 an object of class \code{'gsmar'} generated by \code{fitGSMAR} or \code{GSMAR}, containing
#'   the \strong{unconstrained} model.
#' @param gsmar2 an object of class \code{'gsmar'} generated by \code{fitGSMAR} or \code{GSMAR}, containing
#'   the \strong{constrained} model.
#' @details Performs a likelihood ratio test, testing the null hypothesis that the true parameter value lies
#'   in the constrained parameter space specified by constraints imposed to the model \code{gsmar2}.
#'   Under the null, the test statistic is asymptotically \eqn{\chi^2}-distributed with \eqn{k} degrees of freedom,
#'   \eqn{k} being the difference in the dimensions of the unconstrained and constrained parameter spaces.
#'
#'   Note that this function does \strong{not} verify that the two models are actually nested. Notably, GSMAR models
#'   with different autoregressive orders are not nested, whereas testing models with different numbers of regimes
#'   induce an identification problem and thereby unreliable test results (see the discussion related to Theorem 2
#'   in Virolainen, 2021).
#' @return A list with class "htest" containing the following components:
#'   \item{statistic}{the value of the likelihood ratio statistics.}
#'   \item{parameter}{the degrees of freedom of the likelihood ratio statistic.}
#'   \item{p.value}{the p-value of the test.}
#'   \item{alternative}{a character string describing the alternative hypothesis.}
#'   \item{method}{a character string indicating the type of the test (likelihood ratio test).}
#'   \item{data.name}{a character string giving the names of the supplied models, gsmar1 and gsmar2.}
#'   \item{gsmar1}{the supplied argument gsmar1}
#'   \item{gsmar2}{the supplied argument gsmar2}
#' @seealso \code{\link{Wald_test}}, \code{\link{fitGSMAR}}, \code{\link{GSMAR}}, \code{\link{diagnostic_plot}},
#'  \code{\link{profile_logliks}}, \code{\link{quantile_residual_tests}}, \code{\link{cond_moment_plot}}
#' @inherit Wald_test references
#' @examples
#' \donttest{
#' # GMAR p=1, M=2 model:
#' fit12 <- fitGSMAR(simudata, p=1, M=2, model="GMAR", ncalls=1, seeds=1)
#'
#' # GMAR p=1, M=2 model with AR parameters restricted to be the same in both
#' # regimes:
#' fit12r <- fitGSMAR(simudata, p=1, M=2, model="GMAR", restricted=TRUE,
#'                    ncalls=1, seeds=1)
#'
#' # Test with likelihood ratio test whether the AR parameters are the same in
#' # both regimes:
#' LR_test(fit12, fit12r)
#' }
#' @export

LR_test <- function(gsmar1, gsmar2) {
  # Checks
  check_gsmar(gsmar1, object_name="gsmar1")
  check_gsmar(gsmar2, object_name="gsmar2")
  stopifnot(length(gsmar1$params) > length(gsmar2$params))
  stopifnot(gsmar1$loglik >= gsmar2$loglik)

  # Calculate the likelihood ratio test
  test_stat <- as.numeric(2*(gsmar1$loglik - gsmar2$loglik))
  df <- length(gsmar1$params) - length(gsmar2$params)
  p_value <- pchisq(test_stat, df=df, lower.tail=FALSE)

#  structure(list(gsmar1=gsmar1,
#                 gsmar2=gsmar2,
#                 test_stat=test_stat,
#                 df=df,
#                 p_value=p_value),
#            class="lr")

  dname <- paste(deparse(substitute(gsmar1)), "and", deparse(substitute(gsmar2)))

  structure(list(statistic=c("LR"=test_stat),
                 parameter=c("df"=df),
                 p.value=p_value,
                 alternative=paste("the true parameter does not satisfy the constraints in", deparse(substitute(gsmar2))),
                 data.name=dname,
                 method="Likelihood ratio test",
                 gsmar1=gsmar1,
                 gsmar2=gsmar2),
            class="htest")
}

#' @describeIn GSMAR Log-likelihood method
#' @inheritParams plot.gsmar
#' @param object object of class \code{'gsmar'} created with \code{fitGSMAR} or \code{GSMAR}.
#' @export
logLik.gsmar <- function(object, ...) object$loglik


#' @describeIn GSMAR residuals method to extract multivariate quantile residuals
#' @inheritParams logLik.gsmar
#' @export
residuals.gsmar <- function(object, ...) object$quantile_residuals


#' @describeIn GSMAR summary method, standard errors in brackets
#' @inheritParams logLik.gsmar
#' @export
summary.gsmar <- function(object, ..., digits=2) {
  gsmar <- object
  check_data(gsmar)
  structure(list(gsmar=gsmar,
                 regime_means=get_regime_means(gsmar),
                 regime_vars=get_regime_vars(gsmar),
                 abs_ar_roots=get_ar_roots(gsmar),
                 digits=digits),
            class="gsmarsum")
}

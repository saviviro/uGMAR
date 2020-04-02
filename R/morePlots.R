#' @title Condinional mean plot for GMAR, StMAR, and G-StMAR models
#'
#' @description \code{condmeanPlot} plots the in-sample conditional mean and regimewise conditional means of the model
#'   along with the time series contained in the model (e.g. the time series the model was fitted to).
#'
#' @inheritParams simulateGSMAR
#' @param ... Grapchical parameters passed to FILL IN
#' @details FILL IN
#' @return \code{condmeanPlot} only plots to a graphical device and does not return anything. Numerical values
#'  of the conditional means can be extracted from the model with the dollar sign.
#' @inherit simulateGSMAR references
#' @seealso \code{\link{profile_logliks}}, \code{\link{diagnosticPlot}}, \code{\link{fitGSMAR}}, \code{\link{GSMAR}}, \code{\link{quantileResidualTests}},
#'  \code{\link{quantileResidualPlot}}, \code{\link{simulateGSMAR}}
#' @examples
#' \donttest{
#' # Examples to be written.
#' }
#' @export

condmeanPlot <- function(gsmar, ...) {
  if(is.null(gsmar$data)) stop("The model needs to contain data!")

  invisible(gsmar)
}

#' @title fitGMAR is deprecated
#' @description This function is deprecated! Use \code{fitGSMAR()} instead!
#' @param ... deprecated
#' @export

fitGMAR <- function(...) {
  .Deprecated(msg="'fitGMAR' is deprecated.\nUse 'fitGSMAR' instead.\nApologies for missing backward compatibility!")
  invisible(NULL)
}


#' @title simulateGMAR is deprecated
#' @description This function is deprecated! Use \code{simulateGSMAR()} instead!
#' @param ... deprecated
#' @export

simulateGMAR <- function(...) {
  .Deprecated(msg="'simulateGMAR' is deprecated.\nUse 'simulateGSMAR' instead.\nApologies for the missing backward compatibility!")
  invisible(NULL)
}


#' @title forecastGMAR is deprecated
#' @description This function is deprecated! Use \code{predict.gsmar()} instead!
#' @param ... deprecated
#' @export

forecastGMAR <- function(...) {
  .Deprecated(msg="'forecastGMAR' is deprecated.\nUse 'predict.gsmar' instead.\nApologies for the missing backward compatibility!")
  invisible(NULL)
}


#' @title plotGMAR is deprecated
#' @description This function is deprecated! Use \code{diagnosticPlot()} instead!
#' @param ... deprecated
#' @export

plotGMAR <- function(...) {
  .Deprecated(msg="'plotGMAR' is deprecated.\nUse 'diagnosticPlot' instead.\nApologies for the missing backward compatibility!")
  invisible(NULL)
}

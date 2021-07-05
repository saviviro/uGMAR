#' Spread between 10-Year and 1-Year Treasury rates: T10Y1Y
#'
#' A dataset containing monthly U.S. interest rate spread between the 10-Year Treasury constant
#' maturity rate and 1-Year Treasury constant maturity rate from 1953IV to 2020II.
#'
#' @format A class 'ts' time series object containing 803 observations.
#'
#' @source \url{https://fred.stlouisfed.org/series/GS10} \url{https://fred.stlouisfed.org/series/GS1}
"T10Y1Y"


#' Spread between 10-Year and 1-Year Treasury rates: M10Y1Y
#'
#' A dataset containing monthly U.S. interest rate spread between the 10-Year Treasury constant
#' maturity rate and 1-Year Treasury constant maturity rate from 1982 January to 2020 December.
#'
#' @format A class 'ts' time series object containing 468 observations.
#'
#' @source \url{https://fred.stlouisfed.org/series/GS10} \url{https://fred.stlouisfed.org/series/GS1}
"M10Y1Y"


#' Simulated data
#'
#' A dataset containing 200 observations simulated from a GMAR p=1, M=2 process.
#'
#' @format A numeric vector of length 200.
#'
#' @source Simulated
"simudata"


#' Spread between the 3-month Treasury bill rate and the effective federal funds rate: TBFF
#'
#' A dataset containing the monthly U.S. interest rate spread between the 3-month Treasury bill secondary
#' market rate and the effective federal funds rate from 1954 July to 2019 July (781 observations).
#' This series was studied in the empirical application of Virolainen (forthcoming) introducing the
#' G-StMAR model.
#'
#' @format A class 'ts' time series object containing 781 observations.
#'
#' @references
#'  \itemize{
#'    \item Virolainen S. forthcoming. A mixture autoregressive model based on Gaussian and Student's t-distributions.
#'          Studies in Nonlinear Dynamics & Econometrics, (preprint available as arXiv:2003.05221).
#'  }
#' @source \url{https://fred.stlouisfed.org/series/TB3SMFFM}
"TBFF"



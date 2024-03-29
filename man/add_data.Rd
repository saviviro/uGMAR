% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GSMARconstructor.R
\name{add_data}
\alias{add_data}
\title{Add data to object of class 'gsmar' defining a GMAR, StMAR, or G-StMAR model}
\usage{
add_data(
  data,
  gsmar,
  calc_qresiduals = TRUE,
  calc_cond_moments = TRUE,
  calc_std_errors = FALSE,
  custom_h = NULL
)
}
\arguments{
\item{data}{a numeric vector or class \code{'ts'} object containing the data. \code{NA} values are not supported.}

\item{gsmar}{a class 'gsmar' object, typically generated by \code{fitGSMAR} or \code{GSMAR}.}

\item{calc_qresiduals}{should quantile residuals be calculated? Default is \code{TRUE} iff the model contains data.}

\item{calc_cond_moments}{should conditional means and variances be calculated? Default is \code{TRUE} iff the model contains data.}

\item{calc_std_errors}{should approximate standard errors be calculated?}

\item{custom_h}{A numeric vector with same the length as the parameter vector: i:th element of custom_h is the difference
used in central difference approximation for partial differentials of the log-likelihood function for the i:th parameter.
If \code{NULL} (default), then the difference used for differentiating overly large degrees of freedom parameters
is adjusted to avoid numerical problems, and the difference is \code{6e-6} for the other parameters.}
}
\value{
Returns an object of class 'gsmar' defining the GMAR, StMAR, or G-StMAR model with the data added to the model.
  If the object already contained data, the data will be updated. Does not modify the 'gsmar' object given as argument!
}
\description{
\code{add_data} adds or updates data to object of class '\code{gsmar}' that defines a GMAR, StMAR,
 or G-StMAR model. Also calculates empirical mixing weights, conditional moments, and quantile residuals accordingly.
}
\examples{
# G-StMAR model without data
params42gs <- c(0.04, 1.34, -0.59, 0.54, -0.36, 0.01, 0.06, 1.28, -0.36,
                0.2, -0.15, 0.04, 0.19, 9.75)
gstmar42 <- GSMAR(p=4, M=c(1, 1), params=params42gs, model="G-StMAR")
gstmar42

# Add data to the model
gstmar42 <- add_data(data=M10Y1Y, gsmar=gstmar42)
gstmar42
}
\references{
\itemize{
   \item Kalliovirta L., Meitz M. and Saikkonen P. 2015. Gaussian Mixture Autoregressive model for univariate time series.
           \emph{Journal of Time Series Analysis}, \strong{36}(2), 247-266.
   \item Meitz M., Preve D., Saikkonen P. 2023. A mixture autoregressive model based on Student's t-distribution.
         \emph{Communications in Statistics - Theory and Methods}, \strong{52}(2), 499-515.
   \item Virolainen S. 2022. A mixture autoregressive model based on Gaussian and Student's t-distributions.
         Studies in Nonlinear Dynamics & Econometrics, \strong{26}(4) 559-580.
 }
}
\seealso{
\code{\link{fitGSMAR}}, \code{\link{GSMAR}}, \code{\link{iterate_more}}, \code{\link{get_gradient}},
 \code{\link{get_regime_means}}, \code{\link{swap_parametrization}}, \code{\link{stmar_to_gstmar}}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecatedFunctions.R
\name{condmomentPlot}
\alias{condmomentPlot}
\title{DEPRECATED, USE \code{cond_moment_plot} INSTEAD! Conditional mean or variance plot for GMAR, StMAR, and G-StMAR models}
\usage{
condmomentPlot(gsmar, which_moment = c("mean", "variance"))
}
\arguments{
\item{gsmar}{a class 'gsmar' object, typically generated by \code{fitGSMAR} or \code{GSMAR}.}

\item{which_moment}{should conditional means or variances be plotted?}
}
\value{
\code{cond_moment_plot} only plots to a graphical device and does not return anything. Numerical values
 of the conditional means/variances can be extracted from the model with the dollar sign.
}
\description{
\code{condmomentPlot} plots the one-step in-sample conditional means/variances of the model along with
 the time series contained in the model (e.g. the time series the model was fitted to). Also plots
 the regimewise conditional means/variances multiplied with the mixing weights.
 DEPRECATED, USE \code{cond_moment_plot} INSTEAD!
}
\details{
DEPRECATED, USE \code{cond_moment_plot} INSTEAD!
}
\references{
\itemize{
   \item Galbraith, R., Galbraith, J. 1974. On the inverses of some patterned matrices arising
           in the theory of stationary time series. \emph{Journal of Applied Probability} \strong{11}, 63-71.
   \item Kalliovirta L. (2012) Misspecification tests based on quantile residuals.
           \emph{The Econometrics Journal}, \strong{15}, 358-393.
   \item Kalliovirta L., Meitz M. and Saikkonen P. 2015. Gaussian Mixture Autoregressive model for univariate time series.
           \emph{Journal of Time Series Analysis}, \strong{36}(2), 247-266.
   \item Meitz M., Preve D., Saikkonen P. 2023. A mixture autoregressive model based on Student's t-distribution.
         \emph{Communications in Statistics - Theory and Methods}, \strong{52}(2), 499-515.
   \item Virolainen S. 2022. A mixture autoregressive model based on Gaussian and Student's t-distributions.
         Studies in Nonlinear Dynamics & Econometrics, \strong{26}(4) 559-580.
 }
}
\seealso{
\code{\link{profile_logliks}}, \code{\link{diagnostic_plot}}, \code{\link{fitGSMAR}}, \code{\link{GSMAR}}, \code{\link{quantile_residual_tests}},
 \code{\link{quantileResidualPlot}}
}

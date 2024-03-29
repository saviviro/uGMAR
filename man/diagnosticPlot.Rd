% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecatedFunctions.R
\name{diagnosticPlot}
\alias{diagnosticPlot}
\title{DEPRECATED, USE \code{diagnostic_plot} INSTEAD! Quantile residual based diagnostic plots for GMAR, StMAR, and G-StMAR models}
\usage{
diagnosticPlot(gsmar, nlags = 20, nsimu = 1, plot_indstats = FALSE)
}
\arguments{
\item{gsmar}{a class 'gsmar' object, typically generated by \code{fitGSMAR} or \code{GSMAR}.}

\item{nlags}{a positive integer specifying how many lags should be calculated for the autocorrelation and
conditional heteroscedasticity statistics.}

\item{nsimu}{a positive integer specifying to how many simulated values from the process the covariance
matrix "Omega" (used to compute the tests) should be based on. Larger number of simulations may result
more reliable tests but takes longer to compute. If smaller than data size, then "Omega" will be based
on the given data. Ignored if \code{plot_indstats==FALSE}.}

\item{plot_indstats}{set \code{TRUE} if the individual statistics discussed in Kalliovirta (2012) should be
plotted with their approximate 95\% critical bounds (this may take some time).}
}
\value{
\code{diagnostic_plot} only plots to a graphical device and does not return anything. Use the
 function \code{quantile_residual_tests} in order to obtain the individual statistics.
}
\description{
\code{diagnosticPlot} plots quantile residual time series, normal QQ-plot, autocorrelation function,
 and squared quantile residual autocorrelation function. There is an option to also plot the individual statistics
 associated with the quantile residual tests (for autocorrelation and conditional heteroskedasticity) divided by
 their approximate standard errors with their approximate 95\% critical bounds (see Kalliovirta 2012, Section 3).
 DEPRECATED, USE \code{diagnostic_plot} INSTEAD!
}
\details{
DEPRECATED, USE \code{diagnostic_plot} INSTEAD!

 Sometimes the individual statistics are not plotted because it was not (numerically) possible
 to calculate all the required statistics. This may suggest that the model is misspecified.

 The dashed lines plotted with autocorrelation functions (for quantile residuals and their squares) are
 plus-minus \eqn{1.96*T^{-1/2}} where \eqn{T} is the sample size (minus the \eqn{p} initial values for
 conditional models).
}
\section{Suggested packages}{

  Install the suggested package "gsl" for faster evaluations in the cases of StMAR and G-StMAR models.
  For large StMAR and G-StMAR models with large data the calculations to obtain the individual statistics
  may take a significantly long time without the package "gsl".
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
\code{\link{profile_logliks}}, \code{\link{get_foc}}, \code{\link{fitGSMAR}}, \code{\link{cond_moment_plot}}, \code{\link{quantile_residual_tests}},
 \code{\link{quantile_residual_plot}}, \code{\link{simulate.gsmar}}, \code{\link{LR_test}}, \code{\link{Wald_test}}
}

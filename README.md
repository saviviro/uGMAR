
<!-- README.md is generated from README.Rmd. Please edit that file -->
uGMAR
=====

<!-- badges: start -->
<!-- badges: end -->
The goal of uGMAR is to provide tools for analysing Gaussian mixture autoregressive (GMAR), Student's t mixture Autoregressive (StMAR) and Gaussian and Student's t mixture autoregressive (G-StMAR) models. uGMAR provides functions for unconstrained and constrained maximum likelihood estimation of the model parameters, quantile residual based model diagnostics, simulation from the processes, and forecasting.

Installation
------------

You can install the released version of uGMAR from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("uGMAR")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("saviviro/uGMAR")
```

Example
-------

This is a basic example how to estimate a GMAR model to data. For details about the example data "logVIX" see ?logVIX. The estimation process is computationally heavy and takes advantage of parallel computing. After estimating the model, it's by simple examples how to conduct some further analysis.

``` r
## Estimate a GMAR(1, 2) model to logarithmized VIX data
data(logVIX, package="uGMAR")
fit <- fitGSMAR(data=logVIX, p=1, M=2, model="GMAR")
fit
summary(fit) # Approximate standard errors in brackets
plot(fit)

get_gradient(fit) # The first order condition
get_soc(fit) # The second order condition (eigen values of approximated Hessian)

## Quantile residual diagnostics
quantileResidualPlot(fit)
diagnosticPlot(fit)
qrt <- quantileResidualTests(fit)

## Simulate a sample path from the estimated process
sim <- simulateGSMAR(fit, nsimu=10)

## Forecast future values of the process
predict(fit, n_ahead=10)
```

References
----------

-   Kalliovirta L., Meitz M. and Saikkonen P. 2015. Gaussian Mixture Autoregressive model for univariate time series. *Journal of Time Series Analysis*, **36**, 247-266.
-   Meitz M., Preve D., Saikkonen P. 2018. A mixture autoregressive model based on Student's t-distribution. arXiv:1805.04010 **
    *e**c**o**n*.*E**M*
    **
-   There are currently no published references for the G-StMAR model, but it's a straightforward generalization with theoretical properties similar to the GMAR and StMAR models.

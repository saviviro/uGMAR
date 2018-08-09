<!-- README.md is generated from README.Rmd. Please edit that file -->
uGMAR
=====

The goal of uGMAR is to provide tools to work with Gaussian Mixture Autoregressive (GMAR), Student's t Mixture Autoregressive (StMAR) and Gaussian and Student's t Mixture Autoregressive (G-StMAR) models. G-StMAR is such model that some of its mixture components are similar to the ones that GMAR model uses and some similar to the ones that StMAR model uses. Most importantly uGMAR provides function `fitGSMAR` for two phase maximum likelihood estimation, but it also contains tools for quantile residual based model diagnostics, forecasting and simulations for example. With uGMAR it's easy to apply general linear constraints to the autoregressive parameters or to restrict them to be the same for regimes.

Example
-------

This is a basic example how to estimate a GMAR, StMAR or G-StMAR model to data. The data "VIX", that is used in this example, comes with the package (for details see ?VIX). The estimation process is computationally heavy and uses parallel computing.

``` r
## Estimate GMAR(1, 2) model to VIX data
fit12 <- fitGSMAR(data=VIX, p=1, M=2)
fit12

## Estimate StMAR(1, 1) model to VIX data
fit11t <- fitGSMAR(data=VIX, p=1, M=1, model="StMAR")
fit11t

## Estimate G-StMAR(1, 1, 1) model to VIX data
fit12gs <- fitGSMAR(data=VIX, p=1, M=c(1, 1), model="G-StMAR")
fit12gs
```

References
----------

-   Kalliovirta L., Meitz M. and Saikkonen P. 2015. Gaussian Mixture Autoregressive model for univariate time series. *Journal of Time Series Analysis*, **36**, 247-266.
-   Meitz M., Preve D., Saikkonen P. 2018. A mixture autoregressive model based on Student's t-distribution. arXiv:1805.04010 **\[econ.EM\]**.
-   There are currently no published references for G-StMAR model, but it's a straight forward generalization with theoretical properties similar to GMAR and StMAR models.

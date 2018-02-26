<!-- README.md is generated from README.Rmd. Please edit that file -->
uGMAR
=====

The goal of uGMAR is to provide tools to work with Gaussian Mixture Autoregressive (GMAR), Student's t Mixture Autoregressive (StMAR) and Gaussian and Student's t Mixture Autoregressive (G-StMAR) models. G-StMAR is such model that some of its mixture components are similar to the ones that GMAR model uses and some similar to the ones that StMAR model uses. Most importantly uGMAR provides function `fitGMAR` for two phase maximum likelihood estimation, but it also contains tools for quantile residual based model diagnostics, forecasting and simulation for example. With uGMAR it's easy to apply general linear constraints to the autoregressive parameters or to restrict them to be the same for regimes.

Simple example
--------------

This is a basic example how to estimate a GMAR, StMAR or G-StMAR model to data. I'll use example data "VIX" which comes with the package (for details see ?VIX). The estimation process is computationally heavy and uses parallel computing.

``` r
## Estimate GMAR(1, 2) model to VIX data
fit12 <- fitGMAR(data=VIX, p=1, M=2)

## Estimate StMAR(1, 1) model to VIX data
fit11t <- fitGMAR(data=VIX, p=1, M=1, StMAR=TRUE)

## Estimate G-StMAR(1, 1, 1) model to VIX data
fit12gs <- fitGMAR(data=VIX, p=1, M=c(1, 1), GStMAR=TRUE)
```

Vignette and about references
-----------------------------

See vignette for more detailed introduction to the package and see the references for information about the models. Unfortunately there are not yet articles published considering the Student's t Mixture Autoregressive model. This package is based on working papers considering the model and the references will be updated after the papers have been published.

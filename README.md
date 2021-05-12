
<!-- README.md is generated from README.Rmd. Please edit that file -->

# uGMAR

<!-- badges: start -->

<!-- badges: end -->

uGMAR provides tools for estimating and analyzing Gaussian mixture
autoregressive (GMAR), Student’s t mixture Autoregressive (StMAR) and
Gaussian and Student’s t mixture autoregressive (G-StMAR) models,
including functions for unconstrained and constrained maximum likelihood
estimation of the model parameters, quantile residual based model
diagnostics, simulation from the processes, and forecasting.

## Installation

You can install the released version of uGMAR from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("uGMAR")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("saviviro/uGMAR")
```

## Example

This is a basic example how to estimate GSMAR model and further analyze
it.. The example data is simulated from a GMAR p=1, M=2 process. The
estimation process is computationally demanding and takes advantage of
parallel computing.

``` r
## Estimate a GMAR(1, 2) model and examine the estimates
data(simudata, package="uGMAR")
fit <- fitGSMAR(data=simudata, p=1, M=2, model="GMAR", ncalls=10, seeds=1:10)
fit
summary(fit) # Approximate standard errors in brackets
plot(fit)

get_gradient(fit) # The first order condition
get_soc(fit) # The second order condition (eigenvalues of approximated Hessian)
profile_logliks(fit) # Plot the profile log-likelihood functions

## Quantile residual diagnostics
quantile_residual_plot(fit)
diagnostic_plot(fit)
qrt <- quantile_residual_tests(fit)

## Simulate a sample path from the estimated process
sim <- simulateGSMAR(fit, nsimu=100)
plot.ts(sim$sample)

## Forecast future values of the process
predict(fit, n_ahead=10, pi=c(0.95, 0.8))

# Estimate a GMAR(1, 2) model with the autoregressive coefficients restricted
# to be the same in both regimes:
fitr <- fitGSMAR(data=simudata, p=1, M=2, model="GMAR", restricted=TRUE,
                 ncalls=10, seeds=1:10)

# Test with likelihood ratio tests whether the AR parameters are the same in
# both regimes (see also the function 'Wald_test'):
LR_test(fit, fitr)

# Conditional mean and variance plots:
cond_moment_plot(fit, which_moment="mean")
cond_moment_plot(fit, which_moment="variance")
```

## References

  - Kalliovirta L., Meitz M. and Saikkonen P. 2015. Gaussian Mixture
    Autoregressive model for univariate time series. *Journal of Time
    Series Analysis*, **36**, 247-266.
  - Meitz M., Preve D., Saikkonen P. forthcoming. A mixture
    autoregressive model based on Student’s t-distribution.
    *Communications in Statistics - Theory and Methods*, doi:
    10.1080/03610926.2021.1916531
  - Virolainen S. 2020. A mixture autoregressive model based on Gaussian
    and Student’s t-distributions. arXiv:2003.05221 \[econ.EM\].

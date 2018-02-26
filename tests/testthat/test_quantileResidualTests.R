library(uGMAR)
context("quantileResidualTests")

# Covariance matrices omegas are very sensitive to even very small changes in the quantile residuals. The quantile residuals are
# calculated only with numerical integration if package "gsl" is not available, which will cause these very
# small differences (in 1e-6 or 1e-8 tolerance for eample) that will then accumulate to big differences in Omega.
# And machine precision might also affect it in some cases.
# And a big difference in Omega will cause a "big difference" in the test results.
# This is why these tests are commented out (and used only for development).
test_that("quantileResiduals works", {
  expect_equal(quantileResiduals_int(VIX, 1, 1, c(-3, 0.9, 2), StMAR=FALSE)[13], 2.513058, tolerance=1e-3)
})

# quantileResidualTests calls other functions to deal with constraints so they are not considered here
# params11t <- c(-2, 0.8, 1, 12) # StMAR
# params12 <- c(1.0, 0.9, 0.25, 4.5, 0.7, 3.0, 0.8)
# params23 <- c(2.7, 0.8, -0.06, 0.3, 3.5, 0.8, -0.07, 2.6, 7.2, 0.3, -0.01, 0.1, 0.6, 0.25)
# params23t <- c(1.8, 0.9, -0.06, 0.4, 7.2, 0.3, -0.009, 0.1, 3.1, 0.83, -0.05, 3.7, 0.7, 0.2, 11, 339, 198) # StMAR
# params12r <- c(1.4, 1.8, 0.9, 0.3, 3.3, 0.8)
# params12gs <- c(1.5, 0.8, 1.5, 2.9, 0.8, 1.1, 0.6, 3)
# params13gsr <- c(1.3, 1, 1.4, 0.8, 0.4, 2, 0.2, 0.25, 0.15, 20) # M1=2, M2=1
#
# set.seed(1)
# qrt11t <- quantileResidualTests(VIX[20:50], 1, 1, params11t, StMAR=TRUE, printRes=FALSE, nsimu=20)
# qrt12 <- quantileResidualTests(VIX[1:100], 1, 2, params12, StMAR=FALSE, printRes=FALSE, nsimu=1)
# qrt23 <- quantileResidualTests(VIX[120:220], 2, 3, params23, StMAR=FALSE, lagsAC=c(1), lagsCH=c(3), printRes=FALSE, nsimu=1)
# qrt23t <- quantileResidualTests(VIX[120:220], 2, 3, params23t, StMAR=TRUE, lagsAC=c(2), lagsCH=c(2), printRes=FALSE, nsimu=1)
# qrt12r <- quantileResidualTests(VIX[100:133], 1, 2, params12r, StMAR=FALSE, lagsAC=c(1, 3), lagsCH=c(1, 3),
#                                 restricted=TRUE, printRes=FALSE, nsimu=1)
# qrt12gs <- quantileResidualTests(VIX[1:50], p=1, M=c(1, 1), params=params12gs, GStMAR=TRUE, lagsAC=c(2), lagsCH=c(1), printRes=FALSE, nsimu=1)
# qrt13gsr <- quantileResidualTests(VIX[1:50], p=1, M=c(2, 1), params=params13gsr, GStMAR=TRUE, restricted=TRUE, lagsAC=c(1), lagsCH=c(1), printRes=FALSE, nsimu=1)
#
# test_that("quantile residual test for normality works", {
#   expect_equal(qrt12gs$normality$testStat, 6.040615, tolerance=1e-3)
#   expect_equal(qrt13gsr$normality$testStat, 26.49644, tolerance=1e-3)
#   expect_equal(qrt11t$normality$testStat, 42.0166, tolerance=1e-3)
#   expect_equal(qrt12$normality$testStat, 32.0942, tolerance=1e-3)
#   expect_equal(qrt23$normality$testStat, 63.9019, tolerance=1e-3)
#   expect_equal(qrt23t$normality$testStat, 12.36796, tolerance=1e-3)
#   expect_equal(qrt12r$normality$testStat, 1.574353, tolerance=1e-3)
# })
#
# test_that("quantile residuals tests for autocorrelation work", {
#   expect_equal(qrt12gs$autocorrelation$testStat, 5.136115, tolerance=1e-3)
#   expect_equal(qrt13gsr$autocorrelation$testStat, 1.674596, tolerance=1e-3)
#   expect_equal(qrt11t$autocorrelation$testStat[4], 1.031244, tolerance=1e-3)
#   expect_equal(qrt11t$autocorrelation$indStat[1], 1.295102, tolerance=1e-3)
#   expect_equal(qrt12$autocorrelation$testStat[5], 4.22105, tolerance=1e-3)
#   expect_equal(qrt12$autocorrelation$indStat[1], 0.05724395, tolerance=1e-3)
#   expect_equal(qrt12$autocorrelation$pvalue[6], 0.8732391, tolerance=1e-3)
#   expect_equal(qrt12$autocorrelation$stdError[3], 0.07095191, tolerance=1e-3)
#   expect_equal(qrt23$autocorrelation$testStat, 2.102267, tolerance=1e-3)
#   expect_equal(qrt23$autocorrelation$indStat, 0.1210913, tolerance=1e-3)
#   expect_equal(qrt23$autocorrelation$stdError, 0.08309301, tolerance=1e-3)
#   expect_equal(qrt23t$autocorrelation$testStat, 0.1728601, tolerance=1e-3)
#   expect_equal(qrt23t$autocorrelation$indStat, 0.05101753, tolerance=1e-3)
#   expect_equal(qrt23t$autocorrelation$stdError, 0.3162522, tolerance=1e-3)
#   expect_equal(qrt12r$autocorrelation$testStat[1], 0.4788158, tolerance=1e-3)
#   expect_equal(qrt12r$autocorrelation$indStat[2], -0.06335654, tolerance=1e-3)
#   expect_equal(qrt12r$autocorrelation$testStat[2], 1.296524, tolerance=1e-3)
#   expect_equal(qrt12$autocorrelation$stdError[1], 0.06846993, tolerance=1e-3)
# })
#
# test_that("quantile residual tests for conditional heteroskedasticity work", {
#   expect_equal(qrt12gs$cond.heteroscedasticity$testStat, 0.02008846, tolerance=1e-3)
#   expect_equal(qrt13gsr$cond.heteroscedasticity$indStat, 0.5920743, tolerance=1e-3)
#   expect_equal(qrt11t$cond.heteroscedasticity$testStat[1], 0.2450658, tolerance=1e-3)
#   expect_equal(qrt11t$cond.heteroscedasticity$stdError[3], 6.227949, tolerance=1e-3)
#   expect_equal(qrt12$cond.heteroscedasticity$testStat[1], 0.4655979, tolerance=1e-3)
#   expect_equal(qrt12$cond.heteroscedasticity$pvalue[2], 0.7655539, tolerance=1e-3)
#   expect_equal(qrt12$cond.heteroscedasticity$indStat[3], -0.05662181, tolerance=1e-3)
#   expect_equal(qrt12$cond.heteroscedasticity$stdError[4], 0.1169547, tolerance=1e-3)
#   expect_equal(qrt23$cond.heteroscedasticity$testStat, 13.90113, tolerance=1e-3)
#   expect_equal(qrt23$cond.heteroscedasticity$indStat, -0.3046883, tolerance=1e-3)
#   expect_equal(qrt23$cond.heteroscedasticity$stdError, 0.1015375, tolerance=1e-3)
#   expect_equal(qrt23t$cond.heteroscedasticity$testStat, 5.7698, tolerance=1e-3)
#   expect_equal(qrt23t$cond.heteroscedasticity$indStat, -0.3604339, tolerance=1e-3)
#   expect_equal(qrt23t$cond.heteroscedasticity$stdError, 0.4287505, tolerance=1e-3)
#   expect_equal(qrt12r$cond.heteroscedasticity$testStat[1], 0.07087807, tolerance=1e-3)
#   expect_equal(qrt12r$cond.heteroscedasticity$testStat[2], 0.4333143, tolerance=1e-3)
#   expect_equal(qrt12r$cond.heteroscedasticity$indStat[1], 0.1396549, tolerance=1e-3)
#   expect_equal(qrt12r$cond.heteroscedasticity$stdError[2], 0.4099274, tolerance=1e-3)
# })

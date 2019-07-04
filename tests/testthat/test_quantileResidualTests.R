library(uGMAR)
context("quantileResidualTests")

# The results slightly differ depending on whether numerical integration (without the package "gsl")
# or hypergeometric function is used to calculate the quantile residuals. Also very small differences
# in machine accuracy etc. might accumulate to differences notable enough. As the tests are moreover
# computationally demanding, they are commented out and used for development only.
test_that("quantileResiduals works", {
  expect_equal(quantileResiduals_int(VIX, 1, 1, c(-3, 0.9, 2), model="GMAR")[13], 2.513058, tolerance=1e-3)
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
# stmar11 <- GSMAR(data=VIX[20:50], p=1, M=1, params=params11t, model="StMAR")
# gmar12 <- GSMAR(data=VIX[1:100], p=1, M=2, params=params12, model="GMAR")
# gmar23 <- GSMAR(data=VIX[120:220], p=2, M=3, params=params23, model="GMAR")
# stmar23 <- GSMAR(data=VIX[120:220], p=2, M=3, params=params23t, model="StMAR")
# gmar12r <- GSMAR(data=VIX[100:133], p=1, M=2, params=params12r, model="GMAR", restricted=TRUE)
# gstmar12 <- GSMAR(data=VIX[1:50], p=1, M=c(1, 1), params=params12gs, model="G-StMAR")
# gstmar13r <- GSMAR(data=VIX[1:50], p=1, M=c(2, 1), params=params13gsr, model="G-StMAR", restricted=TRUE)
#
# set.seed(1)
# qrt11t <- quantileResidualTests(stmar11, printRes=FALSE, nsimu=20)
# qrt12 <- quantileResidualTests(gmar12, printRes=FALSE, nsimu=1)
# qrt23 <- quantileResidualTests(gmar23, lagsAC=c(1), lagsCH=c(3), printRes=FALSE, nsimu=1)
# qrt23t <- quantileResidualTests(stmar23, lagsAC=c(3), lagsCH=c(2), printRes=FALSE, nsimu=1)
# qrt12r <- quantileResidualTests(gmar12r, lagsAC=c(1, 3), lagsCH=c(1, 3), printRes=FALSE, nsimu=1)
# qrt12gs <- quantileResidualTests(gstmar12, lagsAC=c(2), lagsCH=c(1), printRes=FALSE, nsimu=1)
# qrt13gsr <- quantileResidualTests(gstmar13r, lagsAC=c(1), lagsCH=c(1), printRes=FALSE, nsimu=1)
#
# #eq_normstat <- function(qrt, val) excpect_equal(qrt$norm_res$testStat, val, tolerance=1e-3) # Return normality test statistic
#
# test_that("quantile residual test for normality works", {
#   expect_equal(qrt12gs$norm_res$testStat, 6.041794, tolerance=1e-3)
#   expect_equal(qrt13gsr$norm_res$testStat, 35.57558, tolerance=1e-3)
#   expect_equal(qrt11t$norm_res$testStat, 42.01667, tolerance=1e-3)
#   expect_equal(qrt12$norm_res$testStat, 32.09423, tolerance=1e-3)
#   expect_equal(qrt23$norm_res$testStat, 63.89753, tolerance=1e-3)
#  # expect_equal(qrt23t$norm_res$testStat, 4.525679, tolerance=1e-3)
#   expect_equal(qrt12r$norm_res$testStat, 1.574353, tolerance=1e-3)
# })
#
# test_that("quantile residuals tests for autocorrelation work", {
#   expect_equal(qrt12gs$ac_res$testStat, 1.227147, tolerance=1e-3)
#   expect_equal(qrt13gsr$ac_res$testStat, 12.54948, tolerance=1e-3)
#   expect_equal(qrt11t$ac_res$testStat, c(8.417805, 8.703875, 12.971856, 45.147962), tolerance=1e-3)
#   expect_equal(qrt11t$ac_res$indStat, c(1.295102, 1.355809, 1.319642, 1.339617), tolerance=1e-3)
#   expect_equal(qrt12$ac_res$testStat, c(0.6927903, 1.0223100, 1.4459217, 5.1346440), tolerance=1e-3)
#   expect_equal(qrt12$ac_res$indStat, c(0.05724395, 0.03412973, 0.01891017, -0.00966325), tolerance=1e-3)
#   expect_equal(qrt12$ac_res$pvalue, c(0.4052169, 0.5998024, 0.9192248, 0.8820040), tolerance=1e-3)
#   expect_equal(qrt12$ac_res$stdError, c(0.06842647, 0.07189492, 0.05256620, 0.06739349), tolerance=1e-3)
#   expect_equal(qrt23$ac_res$testStat, 1.18069, tolerance=1e-3)
#   expect_equal(qrt23$ac_res$indStat, 0.1210913, tolerance=1e-3)
#   expect_equal(qrt23$ac_res$stdError, 0.1108767, tolerance=1e-3)
#   expect_equal(qrt23t$ac_res$testStat, 0.6423269, tolerance=1e-3)
#   expect_equal(qrt23t$ac_res$indStat, -0.02966591, tolerance=1e-3)
#   expect_equal(qrt23t$ac_res$stdError, 0.09152835, tolerance=1e-3)
#   expect_equal(qrt12r$ac_res$testStat, c(0.4761958, 1.4598831), tolerance=1e-3)
#   expect_equal(qrt12r$ac_res$indStat, c(0.26191738, -0.06335654), tolerance=1e-3)
#   expect_equal(qrt12r$ac_res$testStat[2], 1.459883, tolerance=1e-3)
#   expect_equal(qrt12$ac_res$stdError[1], 0.06842647, tolerance=1e-3)
# })
#
# test_that("quantile residual tests for conditional heteroskedasticity work", {
#   expect_equal(qrt12gs$ch_res$testStat, 0.01772478, tolerance=1e-3)
#   expect_equal(qrt13gsr$ch_res$indStat, 0.5920743, tolerance=1e-3)
#   expect_equal(qrt11t$ch_res$testStat, c(0.3393064, 0.4802340, 2.5123071, 1.1701140), tolerance=1e-3)
#   expect_equal(qrt11t$ch_res$stdError, c(0.7286874, 1.0237356, 0.8122795, 0.8237535), tolerance=1e-3)
#   expect_equal(qrt12$ch_res$testStat, c(0.4647988, 0.5310061, 10.1436096, 20.4697598), tolerance=1e-3)
#   expect_equal(qrt12$ch_res$pvalue, c(0.49538973, 0.76682012, 0.07126690, 0.02510988), tolerance=1e-3)
#   expect_equal(qrt12$ch_res$indStat[3], -0.4115146, tolerance=1e-3)
#   expect_equal(qrt12$ch_res$stdError[4], 0.122487, tolerance=1e-3)
#   expect_equal(qrt23$ch_res$testStat, 18.6777, tolerance=1e-3)
#   expect_equal(qrt23$ch_res$indStat, -0.3046883, tolerance=1e-3)
#   expect_equal(qrt23$ch_res$stdError, 0.09125149, tolerance=1e-3)
#   expect_equal(qrt23t$ch_res$testStat, 6.521116, tolerance=1e-3)
#   expect_equal(qrt23t$ch_res$indStat, -0.3604339, tolerance=1e-3)
#   expect_equal(qrt23t$ch_res$stdError, 0.1398849, tolerance=1e-3)
#   expect_equal(qrt12r$ch_res$testStat[1], 0.07133745, tolerance=1e-3)
#   expect_equal(qrt12r$ch_res$testStat[2], 0.43759, tolerance=1e-3)
#   expect_equal(qrt12r$ch_res$indStat, c(0.1396549, -0.2607705), tolerance=1e-3)
#   expect_equal(qrt12r$ch_res$stdError[2], 0.4147108, tolerance=1e-3)
# })

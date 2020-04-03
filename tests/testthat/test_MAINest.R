library(uGMAR)
context("functions in MAINest")

# Note that fitGSMAR is not covered at all
test_that("get_minval works correctly", {
  expect_equal(get_minval(simudata), -9999)
  expect_equal(get_minval(rep(0, 1000)), -9999)
  expect_equal(get_minval(rep(0, 1001)), -99999)
})

params12 <- c(1.7, 0.85, 0.3, 4.12, 0.73, 1.98, 0.63)
gmar12 <- GSMAR(data=simudata, p=1, M=2, params=params12, model="GMAR")

params11t <- c(0.9, 0.92, 1.01, 2.89)
stmar11 <- GSMAR(data=simudata, p=1, M=1, params=params11t, model="StMAR")

params12gs <- c(4.13, 0.73, 1.98, 1.7, 0.85, 0.3, 0.37, 9) # M1=1, M2=1
gstmar12 <- GSMAR(data=simudata, p=1, M=c(1, 1), params=params12gs, model="G-StMAR")


test_that("iterate_more works", {
  gmar12it <- iterate_more(gmar12, calc_std_errors=FALSE, maxit=1)
  stmar11it <- iterate_more(stmar11, calc_std_errors=FALSE, maxit=1)
  gstmar12it <- iterate_more(gstmar12, calc_std_errors=FALSE, maxit=1)
  # We only test for no errors
})


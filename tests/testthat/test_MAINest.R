library(uGMAR)
context("functions in MAINest")

test_that("fitGSMAR does not throw errors", {
  fitGSMAR0 <- function(M, model) {
    suppressMessages(fitGSMAR(simudata, p=1, M=M, model=model, ncalls=1, ncores=1, maxit=1, seeds=1, print_res=FALSE, ngen=2))
  }
  tmp <- fitGSMAR0(M=2, model="GMAR")
  tmp <- fitGSMAR0(M=1, model="StMAR")
  tmp <- fitGSMAR0(M=c(1, 1), model="G-StMAR")
  expect_true(TRUE)
})

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


test_that("iterate_more does not throw errors", {
  iterate_more0 <- function(gsmar) suppressMessages(iterate_more(gsmar, calc_std_errors=FALSE, maxit=1))
  gmar12it <- iterate_more0(gmar12)
  stmar11it <- iterate_more0(stmar11)
  gstmar12it <- iterate_more0(gstmar12)
  expect_true(TRUE)
})


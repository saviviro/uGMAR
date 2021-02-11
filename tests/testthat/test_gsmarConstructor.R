library(uGMAR)
context("GSMAR etc")

params12 <- c(0.18, 0.93, 0.01, 0.86, 0.68, 0.02, 0.88)
params12t <- c(1.38, 0.88, 0.27, 3.8, 0.74, 3.15, 0.8, 300, 3.6)
params42gsr <- c(0.11, 0.03, 1.27, -0.39, 0.24, -0.17, 0.03, 1.01, 0.3, 2.03)

constraints <- list(diag(1, ncol=2, nrow=2), as.matrix(c(1, 0)))
params22c <- c(0.03, 1.27, -0.29, 0.03, -0.01, 0.91, 0.34, 0.88)

test_that("GSMAR does not throw errors", {
  tmp <- GSMAR(p=1, M=2, params=params12, model="GMAR")
  tmp <- GSMAR(p=1, M=2, params=params12t, model="StMAR")
  tmp <- GSMAR(data=T10Y1Y, p=4, M=c(1, 1), params=params42gsr, model="G-StMAR", restricted=TRUE)
  tmp <- GSMAR(data=T10Y1Y, p=2, M=2, params=params22c, model="GMAR", constraints=constraints)
  expect_true(TRUE)
})


test_that("add_data works correctly", {
  gstmar42r <- GSMAR(p=4, M=c(1, 1), params=params42gsr, model="G-StMAR", restricted=TRUE)
  gstmar42r <- add_data(data=T10Y1Y, gstmar42r)
  expect_equal(gstmar42r$data, T10Y1Y)
})


test_that("swap_parametrization works correctly", {
  gstmar42r <- GSMAR(data=T10Y1Y, p=4, M=c(1, 1), params=params42gsr, model="G-StMAR", restricted=TRUE)
  gstmar42r2 <- swap_parametrization(gstmar42r)
  expect_equal(gstmar42r2$params, c(2.20, 0.60, 1.27, -0.39, 0.24, -0.17, 0.03, 1.01, 0.30, 2.03))
})


test_that("stmar_to_gstmar works", {
  stmar12 <- GSMAR(p=1, M=2, params=params12t, model="StMAR")
  gstmar12 <- stmar_to_gstmar(stmar12)
  expect_equal(gstmar12$params, c(1.38, 0.88, 0.27, 3.8, 0.74, 3.15, 0.8, 3.6))
})


test_that("alt_gsmar does not throw errors", {
  fit11 <- suppressMessages(fitGSMAR(T10Y1Y, 1, 1, model="StMAR", ncalls=2, ncores=1, maxit=1, seeds=1:2, print_res=FALSE, ngen=1))
  alt_gsmar(fit11, which_round=1)
  expect_true(TRUE)
})

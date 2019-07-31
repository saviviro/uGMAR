library(uGMAR)
context("standard errors")

params12 <- c(1.0, 0.9, 0.25, 4.5, 0.7, 3.0, 0.8)
params12gs <- c(1.5, 0.8, 1.5, 2.9, 0.8, 1.1, 0.6, 3)

test_that("standardErrors work", {
  expect_equal(standardErrors(VIX, p=1, M=2, params=params12, model="GMAR"),
               c(0.22749623, 0.02002580, 0.03528572, 1.51607790, 0.11058642, 0.54552605, NA), tolerance=1e-3)
  expect_equal(standardErrors(VIX, p=1, M=c(1, 1), params=params12gs, model="G-StMAR"),
               c(NA, NA, NA, 0.28571566, 0.02554265, 0.2051497, NA, 0.27498723), tolerance=1e-3)
})



library(uGMAR)
context("standard errors")

params11t <- c(0.9, 0.92, 1.01, 2.89)
params12 <- c(1.7, 0.85, 0.3, 4.12, 0.73, 1.98, 0.63)
params12gs <- c(4.13, 0.73, 1.98, 1.7, 0.85, 0.3, 0.37, 9) # M1=1, M2=1

test_that("standard_errors work", {
  expect_equal(standard_errors(simudata, p=1, M=1, params=params11t, model="StMAR"),
               c(0.31627663, 0.02582825, 0.28359778, 0.36390871), tolerance=1e-2)
  expect_equal(standard_errors(simudata, p=1, M=2, params=params12, model="GMAR"),
               c(0.63767937, 0.05622560, 0.05026257, 1.48243852, 0.09497417, 0.36394897, 0.15529767), tolerance=1e-2)
  expect_equal(standard_errors(simudata, p=1, M=c(1, 1), params=params12gs, model="G-StMAR"),
               c(1.52047054, 0.09725554, 0.37288996, 0.63467556, 0.05605961, 0.05325949, 0.16981272, 3.51539538), tolerance=1e-2)
})



library(uGMAR)
context("quantileResiduals")

params12 <- c(1.0, 0.9, 0.25, 4.5, 0.7, 3.0, 0.8)
params12t <- c(1.1, 0.9, 0.3, 4.5, 0.7, 3.2, 0.8, 5, 8) # StMAR
params22 <- c(1.2, 0.8, 0.05, 0.3, 3.5, 0.8, -0.1, 2.8, 0.8)
params22t <- c(1.4, 0.8, 0.05, 0.27, 3.5, 0.9, -0.18, 3.1, 0.7, 203, 3) # StMAR
params23 <- c(2.7, 0.8, -0.06, 0.3, 3.5, 0.8, -0.07, 2.6, 7.2, 0.3, -0.01, 0.1, 0.6, 0.25)

params12r <- c(1.4, 1.8, 0.9, 0.3, 3.3, 0.8)
params12tr <- c(0.8, 0.96, 0.9, 0.4, 5.8, 0.9, 4, 272) # StMAR
params23tr <-  c(1.9, 1.6, 2.1, 0.8, -0.02, 0.4, 0.1, 3.9, 0.6, 0.3, 15, 200, 220) # StMAR

params12gs <- c(1.2, 0.8, 0.6, 1.3, 0.6, 1.1, 0.6, 3)
params23gs <- c(1, 0.1, 0.1, 1, 1.2, 0.2, 0.2, 1.2, 1.3, 0.3, -0.3, 1.3, 0.3, 0.4, 11, 12) # M1=1, M2=2
params13gsr <- c(1.3, 2.2, 1.4, 0.8, 2.4, 4.6, 0.4, 0.25, 0.15, 20) # M1=2, M2=1
params14gsr <- c(1.3, 2.2, 1.4, 1.5, 0.8, 2.4, 4.6, 0.4, 1, 0.35, 0.3, 0.2, 20, 30) # M1=2, M2=2

R1 <- matrix(c(1, 0, 0, 0, 0, 1), ncol=2)
R2 <- diag(1, ncol=3, nrow=3)
R3 <- matrix(c(0.5, 0.5), ncol=1)
R4 <- diag(1, ncol=2, nrow=2)

params21c <- c(1, 0.9, 1, 3) # R3, StMAR
params22c <- c(1, 0.1, -0.1, 1, 2, 0.2, 2, 0.8, 11, 12) # R4, R3, StMAR
params22c_2 <- c(1.2, 0.8, 0.05, 0.3, 3.5, 0.8, -0.1, 2.8, 0.8) # R4, R4 (should be same as non-constrained)

params22cr <- c(1, 2, 0.8, 1, 2, 0.7, 11, 12) # R3, StMAR
params32cr <- c(1, 2, 0.3, -0.3, 1, 2, 0.6) # R1
params23cr <- c(1.7, 1.9, 2.1, 0.8, -0.05, 0.3, 0.7, 4.5, 0.7, 0.2) # R4 (should be same as non-constrained)

params32gsc <- c(1, 0.1, 0.1, 1, 2, 0.2, 0.2, 0.2, 2, 0.6, 10) # M1=1, M2=1, R1, R2
params22gsrc <- c(1, 2, 0.5, 1, 2, 0.5, 10) # M1=1, M2=1, R3

test_that("quantileResiduals gives correct residuals for non-restricted models", {
  expect_equal(quantileResiduals_int(VIX, 1, 2, params12, StMAR=FALSE)[1], -0.2574897, tolerance=1e-3)
  expect_equal(quantileResiduals_int(VIX, 1, 2, params12t, StMAR=TRUE)[13], -0.9728816, tolerance=1e-3)
  expect_equal(quantileResiduals_int(VIX[1:10], 2, 2, params22)[1], 0.5079722, tolerance=1e-3)
  expect_equal(quantileResiduals_int(3*VIX, 2, 2, params22, epsilon=-1)[113], 5.578132, tolerance=1e-3)
  expect_equal(quantileResiduals_int(VIX, 2, 2, params22t, StMAR=TRUE)[213], 0.6504481, tolerance=1e-3)
  expect_equal(quantileResiduals_int(VIX, 2, 3, params23)[2], -0.6578804, tolerance=1e-3)
  expect_equal(quantileResiduals_int(3*VIX, 2, 3, params23)[1], 5.053051, tolerance=1e-3)
  expect_equal(quantileResiduals_int(VIX, 1, c(1, 1), params12gs, GStMAR=TRUE)[10], 0.7802472, tolerance=1e-3)
  expect_equal(quantileResiduals_int(VIX, 2, c(1, 2), params23gs, GStMAR=TRUE)[100], 1.712417, tolerance=1e-3)
})

test_that("quantileResiduals gives correct residuals for restricted models", {
  expect_equal(quantileResiduals_int(VIX, 1, 2, params12r, restricted=TRUE)[1], -0.7925692, tolerance=1e-3)
  expect_equal(quantileResiduals_int(VIX, 1, 2, params12r, restricted=TRUE, epsilon=-1)[12], 1.163688, tolerance=1e-3)
  expect_equal(quantileResiduals_int(VIX, 1, 2, params12tr, StMAR=TRUE, restricted=TRUE)[1], -0.0005443089, tolerance=1e-3)
  expect_equal(quantileResiduals_int(VIX, 2, 3, params23tr, StMAR=TRUE, restricted=TRUE)[13], 0.7850545, tolerance=1e-3)
  expect_equal(quantileResiduals_int(VIX, 1, c(2,1), params13gsr, GStMAR=TRUE, restricted=TRUE)[7], 0.143749, tolerance=1e-3)
  expect_equal(quantileResiduals_int(VIX, 1, c(2,2), params14gsr, GStMAR=TRUE, restricted=TRUE)[70], 0.6686324, tolerance=1e-3)
})

test_that("quantileResiduals gives correct residuals for constrained models", {
  expect_equal(quantileResiduals_int(VIX, 2, 1, params21c, StMAR=TRUE, constraints=TRUE, R=list(R3))[240], -0.3792682, tolerance=1e-3)
  expect_equal(quantileResiduals_int(VIX, 2, 2, params22c, StMAR=TRUE, constraints=TRUE, R=list(R4, R3))[1], 1.913685, tolerance=1e-3)
  expect_equal(quantileResiduals_int(3*VIX, 2, 2, params22c_2, constraints=TRUE, R=list(R4, R4))[113], 5.578132, tolerance=1e-3)
  expect_equal(quantileResiduals_int(VIX, 2, 2, params22cr, StMAR=TRUE, restricted=TRUE, constraints=TRUE, R=R3)[24], 0.4021366, tolerance=1e-3)
  expect_equal(quantileResiduals_int(VIX, 3, 2, params32cr, restricted=TRUE, constraints=TRUE, R=R1)[2], 8.125891, tolerance=1e-3)
  expect_equal(quantileResiduals_int(VIX, 2, 3, params23cr, restricted=TRUE, constraints=TRUE, R=R4)[130], 0.4080434, tolerance=1e-3)
  expect_equal(quantileResiduals_int(VIX, 3, c(1,1), params32gsc, GStMAR=TRUE, restricted=FALSE, constraints=TRUE, R=list(R1, R2))[13], 1.219363, tolerance=1e-3)
  expect_equal(quantileResiduals_int(VIX, 2, c(1,1), params22gsrc, GStMAR=TRUE, restricted=TRUE, constraints=TRUE, R=R3)[1], 1.269025, tolerance=1e-3)
})

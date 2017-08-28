library(uGMAR)
context("loglikelihood and mixing weights")

params12 <- c(1.0, 0.9, 0.25, 4.5, 0.7, 3.0, 0.8)
params12t <- c(1.1, 0.9, 0.3, 4.5, 0.7, 3.2, 0.8, 5, 8) # StMAR
params22 <- c(1.2, 0.8, 0.05, 0.3, 3.5, 0.8, -0.1, 2.8, 0.8)
params22t <- c(1.4, 0.8, 0.05, 0.27, 3.5, 0.9, -0.18, 3.1, 0.7, 203, 3) # StMAR
params23 <- c(2.7, 0.8, -0.06, 0.3, 3.5, 0.8, -0.07, 2.6, 7.2, 0.3, -0.01, 0.1, 0.6, 0.25)
params12bound <- c(0.4, 0.4, 0.4, 0.6, 0.6, 0.6, 0.0)
params22bound <- c(1, 1, 0.1, 1, 2, 0.2, 0.2, 0.2, 0.7)

params12r <- c(1.4, 1.8, 0.9, 0.3, 3.3, 0.8)
params12tr <- c(0.8, 0.96, 0.9, 0.4, 5.8, 0.9, 4, 272) # StMAR
params23r <- c(1.7, 1.9, 2.1, 0.8, -0.05, 0.3, 0.7, 4.5, 0.7, 0.2)
params23tr <-  c(1.9, 1.6, 2.1, 0.8, -0.02, 0.4, 0.1, 3.9, 0.6, 0.3, 15, 200, 220) # StMAR
params12rbound <- c(1, 2, 0.2, -0.1, 0.6)
params13rbound <- c(0.1, 0.2, 0.3, 0.5, 0.1, 0.2, 0.3, 0.5, 0.5)

R1 <- matrix(c(1, 0, 0, 0, 0, 1), ncol=2)
R2 <- diag(1, ncol=3, nrow=3)
R3 <- matrix(c(0.5, 0.5), ncol=1)
R4 <- diag(1, ncol=2, nrow=2)

params32c <- c(1, 0.1, -0.1, 1, 2, 0.2, -0.2, 2, 0.6, 11, 12) # R1, R1, StMAR
params33c <- c(1, 0.1, 0.1, 0.1, 1, 2, 0.2, 0.2, 0.2, 2, 3, 0.3, -0.3, 3, 0.5, 0.4) # R2, R2, R1
params21c <- c(1, 0.9, 1, 3) # R3, StMAR
params22c <- c(1, 0.1, -0.1, 1, 2, 0.2, 2, 0.8, 11, 12) # R4, R3, StMAR
params22c_2 <- c(1.2, 0.8, 0.05, 0.3, 3.5, 0.8, -0.1, 2.8, 0.8) # R4, R4 (should be same as non-constrained)

params21cr <- c(1, 1, 1) # R3 bound
params22cr <- c(1, 2, 0.8, 1, 2, 0.7, 11, 12) # R3, StMAR
params32cr <- c(1, 2, 0.3, -0.3, 1, 2, 0.6) # R1
params23cr <- c(1.7, 1.9, 2.1, 0.8, -0.05, 0.3, 0.7, 4.5, 0.7, 0.2) # R4 (should be same as non-constrained)

test_that("Loglikelihood gives correct value for non-restricted models", {
  expect_equal(loglikelihood_int(VIX, 1, 2, params12, StMAR=FALSE), -307.5011, tolerance=1e-3)
  expect_equal(loglikelihood_int(VIX, 1, 2, params12, conditional=FALSE, StMAR=FALSE), -310.944, tolerance=1e-3)
  expect_equal(loglikelihood_int(5*VIX[10:50], 1, 2, params12, conditional=FALSE, StMAR=FALSE), -1965.026, tolerance=1e-3)
  expect_equal(loglikelihood_int(VIX[100:102], 1, 2, params12, StMAR=FALSE), -1.276816, tolerance=1e-3)
  expect_equal(loglikelihood_int(VIX, 1, 2, params12t, StMAR=TRUE), -282.7787, tolerance=1e-3)
  expect_equal(loglikelihood_int(VIX, 1, 2, params12t, StMAR=TRUE, epsilon=-1), -282.7787, tolerance=1e-3)
  expect_equal(loglikelihood_int(VIX, 2, 2, params22, StMAR=FALSE), -389.3221, tolerance=1e-3)
  expect_equal(loglikelihood_int(log(VIX), 2, 2, params22, StMAR=FALSE), -684.4247, tolerance=1e-3)
  expect_equal(loglikelihood_int(VIX[10:20], 2, 2, params22, conditional=FALSE, StMAR=FALSE), -18.44737, tolerance=1e-3)
  expect_equal(loglikelihood_int(VIX, 2, 2, params22t, conditional=FALSE, StMAR=TRUE), -308.9636, tolerance=1e-3)
  expect_equal(loglikelihood_int(-10*VIX, 2, 2, params22t, StMAR=TRUE), -1297.336, tolerance=1e-3)
  expect_equal(loglikelihood_int(VIX, 2, 3, params23, StMAR=FALSE), -329.4227, tolerance=1e-3)
  expect_equal(loglikelihood_int(5*VIX, 2, 3, params23, StMAR=FALSE), -10192.5, tolerance=1e-3)
  expect_equal(loglikelihood_int(VIX[1:4], 2, 3, params23, conditional=FALSE, StMAR=FALSE), -8.210192, tolerance=1e-3)
  expect_equal(loglikelihood_int(VIX, 1, 2, params12bound, boundaries=TRUE, StMAR=FALSE), -99999)
  expect_equal(loglikelihood_int(VIX, 2, 2, params22bound, boundaries=TRUE, StMAR=FALSE), -99999)
})

test_that("Loglikelihood gives correct value for resticted models", {
  expect_equal(loglikelihood_int(VIX, 1, 2, params12r, restricted=TRUE, StMAR=FALSE), -310.2204, tolerance=1e-3)
  expect_equal(loglikelihood_int(log(VIX), 1, 2, params12r, restricted=TRUE, conditional=FALSE, StMAR=FALSE), -482.5679, tolerance=1e-3)
  expect_equal(loglikelihood_int(VIX[13:16], 1, 2, params12r, restricted=TRUE, conditional=FALSE, StMAR=FALSE), -4.759725, tolerance=1e-3)
  expect_equal(loglikelihood_int(7*VIX, 1, 2, params12r, restricted=TRUE, conditional=FALSE, StMAR=FALSE), -4061.793, tolerance=1e-3)
  expect_equal(loglikelihood_int(VIX, 1, 2, params12tr, StMAR=TRUE, restricted=TRUE), -380.6129, tolerance=1e-3)
  expect_equal(loglikelihood_int(log(VIX), 1, 2, params12tr, StMAR=TRUE, restricted=TRUE, conditional=FALSE), -402.2318, tolerance=1e-3)
  expect_equal(loglikelihood_int(VIX, 2, 3, params23r, restricted=TRUE, StMAR=FALSE), -462.2798, tolerance=1e-3)
  expect_equal(loglikelihood_int(VIX[1:5], 2, 3, params23r, restricted=TRUE, conditional=FALSE, StMAR=FALSE), -14.15429, tolerance=1e-3)
  expect_equal(loglikelihood_int(VIX, 2, 3, params23r, restricted=TRUE, conditional=FALSE, StMAR=FALSE), -470.8897, tolerance=1e-3)
  expect_equal(loglikelihood_int(10*VIX[100:150], 2, 3, params23r, restricted=TRUE, conditional=FALSE, StMAR=FALSE), -5231.426, tolerance=1e-3)
  expect_equal(loglikelihood_int(VIX, 2, 3, params23tr, StMAR=TRUE, restricted=TRUE, conditional=FALSE), -412.6474, tolerance=1e-3)
  expect_equal(loglikelihood_int(3*VIX, 2, 3, params23tr, StMAR=TRUE, restricted=TRUE, epsilon=-2), -955.0699, tolerance=1e-3)
  expect_equal(loglikelihood_int(VIX, 1, 2, params12rbound, restricted=TRUE, boundaries=TRUE, StMAR=FALSE), -99999)
  expect_equal(loglikelihood_int(VIX, 1, 3, params13rbound, restricted=TRUE, boundaries=TRUE, StMAR=FALSE), -99999)
})

test_that("Loglikelihood gives correct value for constrained models", {
  expect_equal(loglikelihood_int(VIX, 3, 2, params32c, StMAR=TRUE, constraints=TRUE, R=list(R1, R1), conditional=FALSE), -1312.314, tolerance=1e-3)
  expect_equal(loglikelihood_int(VIX, 3, 3, params33c, constraints=TRUE, R=list(R2, R2, R1)), -977.7254, tolerance=1e-3)
  expect_equal(loglikelihood_int(VIX, 2, 1, params21c, StMAR=TRUE, constraints=TRUE, R=list(R3), epsilon=-1), -315.6233, tolerance=1e-3)
  expect_equal(loglikelihood_int(VIX, 2, 2, params22c, StMAR=TRUE, constraints=TRUE, R=list(R4, R3)), -1098.969, tolerance=1e-3)
  expect_equal(loglikelihood_int(VIX, 2, 2, params22c_2, constraints=TRUE, R=list(R4, R4)), -389.3221, tolerance=1e-3)
  expect_equal(loglikelihood_int(VIX, 2, 1, params21cr, restricted=TRUE, constraints=TRUE, R=R3, boundaries=TRUE), -99999)
  expect_equal(loglikelihood_int(VIX, 2, 2, params22cr, StMAR=TRUE, restricted=TRUE, constraints=TRUE, R=R3), -381.8174, tolerance=1e-3)
  expect_equal(loglikelihood_int(VIX, 3, 2, params32cr, restricted=TRUE, constraints=TRUE, R=R1, conditional=FALSE, epsilon=-1), -7602.361, tolerance=1e-3)
  expect_equal(loglikelihood_int(VIX, 2, 3, params23cr, restricted=TRUE, constraints=TRUE, R=R4), -462.2798, tolerance=1e-3)
})

test_that("mixingWeights gives correct weights for non-restricted models", {
  expect_equal(mixingWeights_int(VIX, 1, 2, params12, StMAR=FALSE)[2, 1], 0.0007360931, tolerance=1e-6)
  expect_equal(mixingWeights_int(VIX, 1, 2, params12t, StMAR=TRUE, epsilon=-1)[20, 2], 0.14159, tolerance=1e-5)
  expect_equal(mixingWeights_int(VIX, 2, 2, params22t, StMAR=TRUE)[1, 1], 1.742628e-07, tolerance=1e-6)
  expect_equal(mixingWeights_int(VIX, 2, 3, params23, StMAR=FALSE)[13, 1], 0.2983124, tolerance=1e-6)
})

test_that("mixingWeights gives correct weights for restricted models", {
  expect_equal(mixingWeights_int(VIX, 1, 2, params12r, StMAR=FALSE, restricted=TRUE)[100, 1], 0.9360905, tolerance=1e-6)
  expect_equal(mixingWeights_int(VIX, 1, 2, params12tr, StMAR=TRUE, restricted=TRUE)[1, 2], 0.9091443, tolerance=1e-6)
  expect_equal(mixingWeights_int(VIX, 2, 3, params23r, StMAR=FALSE, restricted=TRUE)[13, 2], 0.009593899, tolerance=1e-6)
  expect_equal(mixingWeights_int(VIX, 2, 3, params23tr, StMAR=FALSE, restricted=TRUE, epsilon=-1)[111, 2], 1.090888e-22, tolerance=1e-6)
})

test_that("mixingWeights gives correct weights for constrained models", {
  expect_equal(mixingWeights_int(VIX, 3, 2, params32c, StMAR=TRUE, constraints=TRUE, R=list(R1, R1))[1, 1], 0.01096221, tolerance=1e-6)
  expect_equal(mixingWeights_int(VIX, 3, 3, params33c, StMAR=FALSE, constraints=TRUE, R=list(R2, R2, R1))[113, 3], 0.0001194389, tolerance=1e-6)
  expect_equal(mixingWeights_int(VIX, 2, 1, params21c, StMAR=TRUE, constraints=TRUE, R=list(R3))[200, 1], 1, tolerance=1e-6)
  expect_equal(mixingWeights_int(VIX, 2, 2, params22c, StMAR=TRUE, constraints=TRUE, R=list(R4, R3))[2, 2], 0.956518, tolerance=1e-6)
  expect_equal(mixingWeights_int(VIX, 2, 2, params22cr, StMAR=TRUE, restricted=TRUE, constraints=TRUE, R=R3)[1, 2], 0.9994679, tolerance=1e-6)
  expect_equal(mixingWeights_int(VIX, 3, 2, params32cr, StMAR=FALSE, restricted=TRUE, constraints=TRUE, R=R1)[100, 1], 1.07196e-29, tolerance=1e-6)
})

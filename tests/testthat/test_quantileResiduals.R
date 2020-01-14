library(uGMAR)
context("quantileResiduals")

params11 <- c(1, 0.9, 1)
params11t <- c(1, 0.8, 1, 3) # StMAR
params12 <- c(1.0, 0.9, 0.25, 4.5, 0.7, 3.0, 0.8)
params12t <- c(1.1, 0.9, 0.3, 4.5, 0.7, 3.2, 0.8, 5, 8) # StMAR
params22 <- c(1.2, 0.8, 0.05, 0.3, 3.5, 0.8, -0.1, 2.8, 0.8)
params22t <- c(1.4, 0.8, 0.05, 0.27, 3.5, 0.9, -0.18, 3.1, 0.7, 203, 3) # StMAR
params23 <- c(2.7, 0.8, -0.06, 0.3, 3.5, 0.8, -0.07, 2.6, 7.2, 0.3, -0.01, 0.1, 0.6, 0.25)

params12r <- c(1.4, 1.8, 0.9, 0.3, 3.3, 0.8)
params12tr <- c(0.8, 0.96, 0.9, 0.4, 5.8, 0.9, 4, 272) # StMAR
params23tr <-  c(1.9, 1.6, 2.1, 0.8, -0.02, 0.4, 0.1, 3.9, 0.6, 0.3, 15, 200, 220) # StMAR

params12gs <- c(1.2, 0.8, 0.6, 1.3, 0.6, 1.1, 0.6, 3)
params22gs <- c(0, -0.8, -0.199, 1.1, 0, 0.8, 0.199, 1.2, 0.5, 2.001)
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
  expect_equal(quantileResiduals_int(VIX, 1, 1, params11, model="GMAR")[4:6], c(0.660, 0.364, -0.155), tolerance=1e-3)
  expect_equal(quantileResiduals_int(VIX, 1, 1, params11t, model="StMAR")[7:9], c(0.5198048, 0.5428692, 0.5111834), tolerance=1e-4)
  expect_equal(quantileResiduals_int(VIX, 1, 2, params12, model="GMAR")[1:3], c(-0.2574897, -0.1218930, -0.9431908), tolerance=1e-4)
  expect_equal(quantileResiduals_int(VIX, 1, 2, params12t, model="StMAR")[13:17], c(-0.97288160, 1.48243917, 0.39651300, -0.18067370, -0.01133098), tolerance=1e-4)
  expect_equal(quantileResiduals_int(VIX[1:10], 2, 2, params22)[1:2], c(0.5079722, -0.3657400), tolerance=1e-4)
  expect_equal(quantileResiduals_int(3*VIX, 2, 2, params22)[113:115], c(5.578132, 4.201229, 3.801425), tolerance=1e-4)
  expect_equal(quantileResiduals_int(VIX, 2, 2, params22t, model="StMAR")[213:215], c(0.65044814, 0.95905038, -0.03750247), tolerance=1e-4)
  expect_equal(quantileResiduals_int(VIX, 2, 3, params23)[1:2], c(0.2372785, -0.6578804), tolerance=1e-4)
  expect_equal(quantileResiduals_int(3*VIX, 2, 3, params23)[1:3], c(5.053051, 2.367575, 5.296035), tolerance=1e-4)
  expect_equal(quantileResiduals_int(VIX, 1, c(1, 1), params12gs, model="G-StMAR")[10:12], c(0.7802472, 0.8500172, 1.1435514), tolerance=1e-4)
  expect_equal(quantileResiduals_int(VIX, 2, c(1, 2), params23gs, model="G-StMAR")[100:102], c(1.712417, 1.899989, 1.763852), tolerance=1e-4)
  expect_equal(quantileResiduals_int(VIX, 2, c(1, 1), params22gs, model="G-StMAR")[20:22], c(-0.5182896, -1.1170781, -2.0246423), tolerance=1e-4)

})

test_that("quantileResiduals gives correct residuals for restricted models", {
  expect_equal(quantileResiduals_int(VIX, 1, 2, params12r, restricted=TRUE)[1:3], c(-0.7925692, -0.1961077, -2.2397076), tolerance=1e-4)
  expect_equal(quantileResiduals_int(VIX, 1, 2, params12r, restricted=TRUE)[12:14], c(1.163688, -1.467829, 1.058580), tolerance=1e-4)
  expect_equal(quantileResiduals_int(VIX, 1, 2, params12tr, model="StMAR", restricted=TRUE)[2:4], c(0.1618962, -0.4549027, 0.4024535), tolerance=1e-4)
  expect_equal(quantileResiduals_int(VIX, 2, 3, params23tr, model="StMAR", restricted=TRUE)[16:18], c(0.33635368, -0.06784057, 0.62990153), tolerance=1e-4)
  expect_equal(quantileResiduals_int(VIX, 1, c(2, 1), params13gsr, model="G-StMAR", restricted=TRUE)[7:8], c(0.1437490, 0.1659238), tolerance=1e-4)
  expect_equal(quantileResiduals_int(VIX, 1, c(2, 2), params14gsr, model="G-StMAR", restricted=TRUE)[70:72], c(0.6686324, 0.8246961, 0.2662351), tolerance=1e-4)
})

test_that("quantileResiduals gives correct residuals for constrained models", {
  expect_equal(quantileResiduals_int(VIX, 2, 1, params21c, model="StMAR", constraints=list(R3))[240:241], c(-0.3792682, 1.0279554), tolerance=1e-4)
  expect_equal(quantileResiduals_int(VIX, 2, 2, params22c, model="StMAR", constraints=list(R4, R3))[1:2], c(1.913685, 1.699737), tolerance=1e-4)
  expect_equal(quantileResiduals_int(3*VIX, 2, 2, params22c_2, constraints=list(R4, R4))[113:115], c(5.578132, 4.201229, 3.801425), tolerance=1e-4)
  expect_equal(quantileResiduals_int(VIX, 2, 2, params22cr, model="StMAR", restricted=TRUE, constraints=R3)[24:25], c(0.4021366, 0.6204558), tolerance=1e-4)
  expect_equal(quantileResiduals_int(VIX, 3, 2, params32cr, restricted=TRUE, constraints=R1)[2:3], c(8.125891, 8.125891), tolerance=1e-1)
  expect_equal(quantileResiduals_int(VIX, 2, 3, params23cr, restricted=TRUE, constraints=R4)[130:131], c(0.4080434, 0.3146683), tolerance=1e-4)
  expect_equal(quantileResiduals_int(VIX, 3, c(1, 1), params32gsc, model="G-StMAR", restricted=FALSE, constraints=list(R1, R2))[13:14], c(1.219363, 1.113442), tolerance=1e-4)
  expect_equal(quantileResiduals_int(VIX, 2, c(1, 1), params22gsrc, model="G-StMAR", restricted=TRUE, constraints=R3)[1:2], c(1.2690254, 0.9798332), tolerance=1e-4)
})

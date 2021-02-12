library(uGMAR)
context("loglikelihood, mixing weights and cond moments")

params11 <- c(2.1, 0.7, 1.1)
params11t <- c(0.9, 0.92, 1.01, 2.89)
params12 <- c(1.7, 0.85, 0.3, 4.12, 0.73, 1.98, 0.63)
params12t <- c(1.1, 0.9, 0.3, 4.5, 0.7, 3.2, 0.8, 5, 8) # StMAR
params22 <- c(1.2, 0.8, 0.05, 0.3, 3.5, 0.8, -0.1, 2.8, 0.8)
params22t <- c(1.4, 0.8, 0.05, 0.27, 3.5, 0.9, -0.18, 3.1, 0.7, 203, 3) # StMAR
params23 <- c(2.7, 0.8, -0.06, 0.3, 3.5, 0.8, -0.07, 2.6, 7.2, 0.3, -0.01, 0.1, 0.6, 0.25)
params12bound <- c(0.4, 0.4, 0.4, 0.6, 0.6, 0.6, 0.0)
params22bound <- c(1, 1, 0.1, 1, 2, 0.2, 0.2, 0.2, 0.7)

params12r <- c(2.18, 2.93, 0.81, 0.29, 1.9, 0.64)
params12tr <- c(0.8, 0.96, 0.9, 0.4, 5.8, 0.9, 4, 272) # StMAR
params23r <- c(1.7, 1.9, 2.1, 0.8, -0.05, 0.3, 0.7, 4.5, 0.7, 0.2)
params23tr <-  c(1.9, 1.6, 2.1, 0.8, -0.02, 0.4, 0.1, 3.9, 0.6, 0.3, 15, 13, 10) # StMAR
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

params12gs <- c(4.13, 0.73, 1.98, 1.7, 0.85, 0.3, 0.37, 9) # M1=1, M2=1
params23gs <- c(1, 0.1, 0.1, 1, 1.2, 0.2, 0.2, 1.2, 1.3, 0.3, -0.3, 1.3, 0.3, 0.4, 11, 12) # M1=1, M2=2
params13gsr <- c(4.8, 3.31, 3.74, 0.69, 2, 0.19, 0.41, 0.34, 0.3, 9) # M1=2, M2=1

params12gs2 <- c(1.5, 0.8, 1.5, 2.9, 0.8, 1.1, 0.6, 3) # M1=1, M2=1
params13gsr2 <- c(1.3, 1, 1.4, 0.8, 0.4, 2, 0.2, 0.25, 0.15, 11) # M1=2, M2=1

params32gsc <- c(1, 0.1, 0.1, 1, 2, 0.2, 0.2, 0.2, 2, 0.6, 10) # M1=1, M2=1, R1, R2
params22gsrc <- c(1, 2, 0.5, 1, 2, 0.5, 10) # M1=1, M2=1, R3

test_that("Loglikelihood gives correct value for non-restricted models", {
  expect_equal(loglikelihood_int(simudata, p=1, M=1, params=params11, model="GMAR", constraints=NULL), -549.1883, tolerance=1e-3)
  expect_equal(loglikelihood_int(simudata, p=1, M=1, params=params11t, model="StMAR", constraints=NULL), -249.2031, tolerance=1e-3)
  expect_equal(loglikelihood_int(simudata, p=1, M=2, params=params12, model="GMAR", constraints=NULL), -238.8676, tolerance=1e-3)
  expect_equal(loglikelihood_int(simudata, p=1, M=2, params=params12, conditional=FALSE, model="GMAR"), -240.85, tolerance=1e-3)
  expect_equal(loglikelihood_int(5*simudata[10:50], p=1, M=2, params=params12, conditional=FALSE, model="GMAR"), -1829.738, tolerance=1e-3)
  expect_equal(loglikelihood_int(simudata[100:102], p=1, M=2, params=params12, model="GMAR"), -2.528468, tolerance=1e-3)
  expect_equal(loglikelihood_int(simudata, p=1, M=2, params=params12t, model="StMAR"), -245.3812, tolerance=1e-3)
  expect_equal(loglikelihood_int(simudata, p=1, M=2, params=params12t, model="StMAR"), -245.3812, tolerance=1e-3)
  expect_equal(loglikelihood_int(simudata, p=2, M=2, params=params22, model="GMAR"), -309.5868, tolerance=1e-3)
  expect_equal(loglikelihood_int(log(simudata), p=2, M=2, params=params22, model="GMAR"), -542.2795, tolerance=1e-3)
  expect_equal(loglikelihood_int(simudata[10:20], p=2, M=2, params=params22, conditional=FALSE, model="GMAR"), -18.19554, tolerance=1e-3)
  expect_equal(loglikelihood_int(simudata, p=2, M=2, params=params22t, conditional=FALSE, model="StMAR"), -261.7327, tolerance=1e-3)
  expect_equal(loglikelihood_int(-10*simudata, p=2, M=2, params=params22t, model="StMAR"), -1027.096, tolerance=1e-3)
  expect_equal(loglikelihood_int(simudata, p=2, M=3, params=params23, model="GMAR"), -255.7899, tolerance=1e-3)
  expect_equal(loglikelihood_int(5*simudata, p=2, M=3, params=params23, model="GMAR"), -8197.359, tolerance=1e-3)
  expect_equal(loglikelihood_int(simudata[1:4], p=2, M=3, params=params23, conditional=FALSE, model="GMAR"), -7.533156, tolerance=1e-3)
  expect_equal(loglikelihood_int(simudata, p=1, M=2, params=params12bound, boundaries=TRUE, model="GMAR", minval=-9999), -9999)
  expect_equal(loglikelihood_int(simudata, p=2, M=2, params=params22bound, boundaries=TRUE, model="GMAR", minval=-9999), -9999)
  expect_equal(loglikelihood_int(simudata, p=1, M=c(1, 1), params=params12gs, model="G-StMAR", conditional=FALSE), -242.5499, tolerance=1e-3)
  expect_equal(loglikelihood_int(simudata, p=2, M=c(1, 2), params=params23gs, model="G-StMAR"), -778.8025, tolerance=1e-3)
  expect_equal(loglikelihood_int(simudata, p=1, M=c(1, 1), params=params12gs2, model="G-StMAR"), -315.1921, tolerance=1e-3)

  expect_equal(loglikelihood_int(-10*simudata, p=2, M=2, params=params22t, model="StMAR", parametrization="mean"), -1014.057, tolerance=1e-3)
  expect_equal(loglikelihood_int(simudata, p=2, M=3, params=params23, model="GMAR", parametrization="mean"), -545.4563, tolerance=1e-3)
  expect_equal(loglikelihood_int(-0.7*simudata, p=1, M=c(1, 1), params=params12gs, model="G-StMAR", conditional=FALSE), -538.9894, tolerance=1e-3)
  expect_equal(loglikelihood_int(simudata, p=2, M=c(1, 2), params=params23gs, model="G-StMAR", parametrization="mean"), -842.6103, tolerance=1e-3)
  expect_equal(loglikelihood_int(-0.3*simudata, p=2, M=c(1, 2), params=params23gs, model="G-StMAR"), -795.1153, tolerance=1e-3)
})


test_that("Loglikelihood gives correct value for resticted models", {
  expect_equal(loglikelihood_int(simudata, p=1, M=2, params=params12r, restricted=TRUE, model="GMAR"), -240.2778, tolerance=1e-3)
  expect_equal(loglikelihood_int(log(simudata), p=1, M=2, params=params12r, restricted=TRUE, conditional=FALSE, model="GMAR"), -579.5836, tolerance=1e-3)
  expect_equal(loglikelihood_int(simudata[13:16], p=1, M=2, params=params12r, restricted=TRUE, conditional=FALSE, model="GMAR"), -3.675981, tolerance=1e-3)
  expect_equal(loglikelihood_int(7*simudata, p=1, M=2, params=params12r, restricted=TRUE, conditional=FALSE, model="GMAR"), -13265.4, tolerance=1e-3)
  expect_equal(loglikelihood_int(simudata, p=1, M=2, params=params12tr, model="StMAR", restricted=TRUE), -299.5897, tolerance=1e-3)
  expect_equal(loglikelihood_int(log(simudata), p=1, M=2, params=params12tr, model="StMAR", restricted=TRUE, conditional=FALSE), -320.248, tolerance=1e-3)
  expect_equal(loglikelihood_int(simudata, p=2, M=3, params=params23r, restricted=TRUE, model="GMAR"), -365.42, tolerance=1e-3)
  expect_equal(loglikelihood_int(simudata[1:5], p=2, M=3, params=params23r, restricted=TRUE, conditional=FALSE, model="GMAR"), -12.76254, tolerance=1e-3)
  expect_equal(loglikelihood_int(simudata, p=2, M=3, params=params23r, restricted=TRUE, conditional=FALSE, model="GMAR"), -372.5312, tolerance=1e-3)
  expect_equal(loglikelihood_int(10*simudata[100:150], p=2, M=3, params=params23r, restricted=TRUE, conditional=FALSE, model="GMAR"), -5005.734, tolerance=1e-3)
  expect_equal(loglikelihood_int(simudata, p=2, M=3, params=params23tr, model="StMAR", restricted=TRUE, conditional=FALSE), -319.8685, tolerance=1e-3)
  expect_equal(loglikelihood_int(3*simudata, p=2, M=3, params=params23tr, model="StMAR", restricted=TRUE), -669.3697, tolerance=1e-3)
  expect_equal(loglikelihood_int(simudata, p=1, M=2, params=params12rbound, restricted=TRUE, boundaries=TRUE, model="GMAR", minval=-9999), -9999)
  expect_equal(loglikelihood_int(simudata, p=1, M=3, params=params13rbound, restricted=TRUE, boundaries=TRUE, model="GMAR", minval=-9999), -9999)
  expect_equal(loglikelihood_int(simudata, p=1, M=c(2, 1), params=params13gsr, model="G-StMAR", restricted=TRUE, conditional=FALSE), -239.2833, tolerance=1e-3)
  expect_equal(loglikelihood_int(simudata, p=1, M=c(2, 1), params=params13gsr2, model="G-StMAR", restricted=TRUE), -409.1438, tolerance=1e-3)

  expect_equal(loglikelihood_int(3*simudata, p=2, M=3, params=params23tr, model="StMAR", restricted=TRUE, parametrization="mean"), -711.2051, tolerance=1e-3)
  expect_equal(loglikelihood_int(-0.6*simudata, p=1, M=c(2, 1), params=params13gsr2, model="G-StMAR", restricted=TRUE), -573.1651, tolerance=1e-3)
})


test_that("Loglikelihood gives correct value for constrained models", {
  expect_equal(loglikelihood_int(simudata, p=3, M=c(1, 1), params=params32gsc, model="G-StMAR", constraints=list(R1, R2), conditional=FALSE), -532.7867, tolerance=1e-3)
  expect_equal(loglikelihood_int(simudata, p=2, M=c(1, 1), params=params22gsrc, model="G-StMAR", restricted=TRUE, constraints=R3, conditional=FALSE), -623.5345, tolerance=1e-3)
  expect_equal(loglikelihood_int(simudata, p=3, M=2, params=params32c, model="StMAR", constraints=list(R1, R1), conditional=FALSE), -1041.463, tolerance=1e-3)
  expect_equal(loglikelihood_int(simudata, p=3, M=3, params=params33c, constraints=list(R2, R2, R1)), -793.5932, tolerance=1e-3)
  expect_equal(loglikelihood_int(simudata, p=2, M=1, params=params21c, model="StMAR", constraints=list(R3)), -274.3698, tolerance=1e-3)
  expect_equal(loglikelihood_int(simudata, p=2, M=2, params=params22c, model="StMAR", constraints=list(R4, R3)), -870.892, tolerance=1e-3)
  expect_equal(loglikelihood_int(simudata, p=2, M=2, params=params22c_2, constraints=list(R4, R4)), -309.5868, tolerance=1e-3)
  expect_equal(loglikelihood_int(simudata, p=2, M=1, params=params21cr, restricted=TRUE, constraints=R3, boundaries=TRUE, minval=-9999), -9999)
  expect_equal(loglikelihood_int(simudata, p=2, M=2, params=params22cr, model="StMAR", restricted=TRUE, constraints=R3), -314.7658, tolerance=1e-3)
  expect_equal(loglikelihood_int(simudata, p=3, M=2, params=params32cr, restricted=TRUE, constraints=R1, conditional=FALSE), -5990.307, tolerance=1e-3)
  expect_equal(loglikelihood_int(simudata, p=2, M=3, params=params23cr, restricted=TRUE, constraints=R4), -365.42, tolerance=1e-3)

  expect_equal(loglikelihood_int(3*simudata, p=3, M=2, params=params32c, model="StMAR", constraints=list(R1, R1), conditional=FALSE), -1325.383, tolerance=1e-3)
  expect_equal(loglikelihood_int(7*simudata, p=2, M=2, params=params22cr, model="StMAR", restricted=TRUE, constraints=R3), -851.2123, tolerance=1e-3)
  expect_equal(loglikelihood_int(0.7*simudata, p=3, M=2, params=params32cr, restricted=TRUE, constraints=R1, conditional=FALSE), -2643.472, tolerance=1e-3)
})


test_that("mixing_weights gives correct values for non-restricted models", {
  expect_equal(mixing_weights_int(simudata, p=1, M=c(1, 1), params=params12gs, model="G-StMAR")[10, 2], 0.9631588, tolerance=1e-3)
  expect_equal(mixing_weights_int(simudata, p=2, M=c(1, 2), params=params23gs, model="G-StMAR")[20, 3], 0.2213917, tolerance=1e-3)
  expect_equal(mixing_weights_int(simudata, p=1, M=2, params=params12, model="GMAR")[2, 1], 0.656533, tolerance=1e-3)
  expect_equal(mixing_weights_int(simudata, p=2, M=2, params=params22t, model="StMAR")[1, 1], 0.007149177, tolerance=1e-4)
  expect_equal(mixing_weights_int(simudata, p=2, M=3, params=params23, model="GMAR")[13, 1], 0.8337944, tolerance=1e-3)
})

test_that("mixing_weights gives correct values for restricted models", {
  expect_equal(mixing_weights_int(simudata, p=1, M=c(2, 1), params=params13gsr2, model="G-StMAR", restricted=TRUE)[20, 2], 0.9781194, tolerance=1e-3)
  expect_equal(mixing_weights_int(simudata, p=1, M=c(2, 1), params=params13gsr, model="G-StMAR", restricted=TRUE)[40, 1], 0.01115213, tolerance=1e-3)
  expect_equal(mixing_weights_int(simudata, p=2, M=3, params=params23r, model="GMAR", restricted=TRUE)[13, 2], 0.1214211, tolerance=1e-6)
  expect_equal(mixing_weights_int(simudata, p=2, M=3, params=params23tr, model="StMAR", restricted=TRUE)[5, 3], 0.9339346, tolerance=1e-6)
})

test_that("mixing_weights gives correct values for constrained models", {
  expect_equal(mixing_weights_int(simudata, p=3, M=c(1, 1), params=params32gsc, model="G-StMAR", constraints=list(R1, R2))[1, 2], 1, tolerance=1e-3)
  expect_equal(mixing_weights_int(simudata, p=2, M=c(1, 1), params=params22gsrc, model="G-StMAR", restricted=TRUE, constraints=R3)[1, 2], 1, tolerance=1e-3)
  expect_equal(mixing_weights_int(simudata, p=3, M=2, params=params32c, model="StMAR", constraints=list(R1, R1))[1, 1], 0.008782872, tolerance=1e-6)
  expect_equal(mixing_weights_int(simudata, p=3, M=3, params=params33c, model="GMAR", constraints=list(R2, R2, R1))[113, 3], 0.0001672733, tolerance=1e-6)
  expect_equal(mixing_weights_int(simudata, p=2, M=1, params=params21c, model="StMAR", constraints=list(R3))[100, 1], 1, tolerance=1e-3)
  expect_equal(mixing_weights_int(simudata, p=3, M=2, params=params32cr, model="GMAR", restricted=TRUE, constraints=R1)[150, 2], 1, tolerance=1e-3)
})

test_that("cond_moments gives correct values for non-restricted models", {
  expect_equal(cond_moments(simudata, p=1, M=c(1, 1), params=params12gs2, model="G-StMAR", to_return="regime_cmeans")[10,], c(10.47718, 11.87718), tolerance=1e-3)
  expect_equal(cond_moments(simudata, p=1, M=c(1, 1), params=params12gs2, model="G-StMAR", to_return="regime_cvars")[10,], c(1.500000, 2.484769), tolerance=1e-3)
  expect_equal(cond_moments(simudata, p=1, M=c(1, 1), params=params12gs, model="G-StMAR", to_return="total_cmeans")[1:2], c(12.70827, 13.03498), tolerance=1e-3)
  expect_equal(cond_moments(simudata, p=1, M=c(1, 1), params=params12gs, model="G-StMAR", to_return="total_cvars")[111:112], c(1.2127486, 0.8788069), tolerance=1e-3)
  expect_equal(cond_moments(simudata, p=2, M=c(1, 2), params=params23gs, model="G-StMAR", to_return="regime_cmeans")[20,], c(3.396052, 5.992103, 1.216248), tolerance=1e-3)
  expect_equal(cond_moments(simudata, p=2, M=c(1, 2), params=params23gs, model="G-StMAR", to_return="regime_cvars")[20,], c(1.00000, 14.02534, 14.39489), tolerance=1e-3)
  expect_equal(cond_moments(simudata, p=2, M=c(1, 2), params=params23gs, model="G-StMAR", to_return="total_cvars")[20], 18.03888, tolerance=1e-3)
  expect_equal(cond_moments(simudata, p=1, M=2, params=params12, model="GMAR", to_return="regime_cmeans")[2,], c(12.70238, 13.56910), tolerance=1e-3)
  expect_equal(cond_moments(simudata, p=1, M=2, params=params12, model="GMAR", to_return="regime_cvars")[21,], c(0.30, 1.98), tolerance=1e-3)
  expect_equal(cond_moments(simudata, p=2, M=3, params=params23, model="GMAR", to_return="total_cmeans")[13], 10.74125, tolerance=1e-4)
  expect_equal(cond_moments(simudata, p=2, M=3, params=params23, model="GMAR", to_return="total_cvars")[100], 2.480149, tolerance=1e-4)
  expect_equal(cond_moments(simudata, p=1, M=2, params=params12t, model="StMAR", to_return="total_cmeans")[1], 12.6818, tolerance=1e-4)
  expect_equal(cond_moments(simudata, p=1, M=2, params=params12t, model="StMAR", to_return="total_cvars")[2], 1.311637, tolerance=1e-4)
  expect_equal(cond_moments(simudata, p=2, M=2, params=params22t, model="StMAR", to_return="regime_cmeans")[1,], c(12.38894, 12.86807), tolerance=1e-4)
  expect_equal(cond_moments(simudata, p=2, M=2, params=params22t, model="StMAR", to_return="regime_cvars")[2,], c(0.2964307, 1.3786931), tolerance=1e-4)
})

test_that("cond_moments gives correct values for restricted models", {
  expect_equal(cond_moments(simudata, p=1, M=c(2, 1), params=params13gsr2, model="G-StMAR", restricted=TRUE, to_return="regime_cmeans")[30,], c(10.61347, 10.31347, 10.71347), tolerance=1e-3)
  expect_equal(cond_moments(simudata, p=1, M=c(2, 1), params=params13gsr2, model="G-StMAR", restricted=TRUE, to_return="regime_cvars")[30,], c(0.4000, 2.000, 0.9556815), tolerance=1e-3)
  expect_equal(cond_moments(simudata, p=1, M=c(2, 1), params=params13gsr, model="G-StMAR", restricted=TRUE, to_return="total_cmeans")[40], 10.33982, tolerance=1e-4)
  expect_equal(cond_moments(simudata, p=1, M=c(2, 1), params=params13gsr, model="G-StMAR", restricted=TRUE, to_return="total_cvars")[40], 0.2860245, tolerance=1e-4)
  expect_equal(cond_moments(simudata, p=1, M=2, params=params12r, model="GMAR", restricted=TRUE, to_return="total_cmeans")[100:101], c(12.79224, 13.15661), tolerance=1e-4)
  expect_equal(cond_moments(simudata, p=1, M=2, params=params12tr, model="StMAR", restricted=TRUE, to_return="total_cvars")[11:12], c(3.319071, 4.097902), tolerance=1e-4)
  expect_equal(cond_moments(simudata, p=2, M=3, params=params23r, model="GMAR", restricted=TRUE, to_return="regime_cmeans")[13,], c(9.744875, 9.944875, 10.144875), tolerance=1e-4)
  expect_equal(cond_moments(simudata, p=2, M=3, params=params23tr, model="StMAR", restricted=TRUE, to_return="regime_cvars")[111,], c(0.8589741, 1.1162082, 3.5958124), tolerance=1e-4)
})

test_that("cond_moments gives correct values for constrained models", {
  expect_equal(cond_moments(simudata, p=3, M=c(1, 1), params=params32gsc, model="G-StMAR", constraints=list(R1, R2), to_return="total_cmeans")[1], 9.890363, tolerance=1e-2)
  expect_equal(cond_moments(simudata, p=2, M=c(1, 1), params=params22gsrc, model="G-StMAR", restricted=TRUE, constraints=R3, to_return="total_cvars")[1], 11.30548, tolerance=1e-4)
  expect_equal(cond_moments(simudata, p=3, M=2, params=params32c, model="StMAR", constraints=list(R1, R1), to_return="regime_cmeans")[1,], c(1.115771, 2.231541), tolerance=1e-3)
  expect_equal(cond_moments(simudata, p=3, M=3, params=params33c, model="GMAR", constraints=list(R2, R2, R1), to_return="total_cvars")[113], 2.006294, tolerance=1e-4)
  expect_equal(cond_moments(simudata, p=2, M=1, params=params21c, model="StMAR", constraints=list(R3), to_return="regime_cvars")[100,], 1.228252, tolerance=1e-4)
  expect_equal(cond_moments(simudata, p=2, M=2, params=params22c, model="StMAR", constraints=list(R4, R3), to_return="regime_cmeans")[2, ], c(1.088879, 4.677675), tolerance=1e-3)
  expect_equal(cond_moments(simudata, p=2, M=2, params=params22cr, model="StMAR", restricted=TRUE, constraints=R3, to_return="total_cmeans")[10:11], c(11.21663, 11.71482), tolerance=1e-4)
  expect_equal(cond_moments(simudata, p=3, M=2, params=params32cr, model="GMAR", restricted=TRUE, constraints=R1, to_return="total_cvars")[100], 2, tolerance=1e-4)
})

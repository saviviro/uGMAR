library(uGMAR)
context("quantile_residuals")

params11 <- c(1, 0.9, 1)
params11t <- c(0.9, 0.92, 1.01, 2.89)
params12 <- c(1.7, 0.85, 0.3, 4.12, 0.73, 1.98, 0.63)
params12t <- c(1.1, 0.9, 0.3, 4.5, 0.7, 3.2, 0.8, 5, 8) # StMAR
params22 <- c(1.2, 0.8, 0.05, 0.3, 3.5, 0.8, -0.1, 2.8, 0.8)
params22t <- c(1.4, 0.8, 0.05, 0.27, 3.5, 0.9, -0.18, 3.1, 0.7, 13, 3) # StMAR
params23 <- c(2.7, 0.8, -0.06, 0.3, 3.5, 0.8, -0.07, 2.6, 7.2, 0.3, -0.01, 0.1, 0.6, 0.25)

params12r <- c(1.4, 1.8, 0.9, 0.3, 3.3, 0.8)
params12tr <- c(0.8, 0.96, 0.9, 0.4, 5.8, 0.9, 4, 15) # StMAR
params23tr <-  c(1.9, 1.6, 2.1, 0.8, -0.02, 0.4, 0.1, 3.9, 0.6, 0.3, 15, 12, 13) # StMAR

params12gs <- c(4.13, 0.73, 1.98, 1.7, 0.85, 0.3, 0.37, 9) # M1=1, M2=1
params22gs <- c(0, -0.8, -0.19, 1.1, 0, 0.8, 0.19, 1.2, 0.5, 2.01)
params23gs <- c(1, 0.1, 0.1, 1, 1.2, 0.2, 0.2, 1.2, 1.3, 0.3, -0.3, 1.3, 0.3, 0.4, 11, 12) # M1=1, M2=2
params13gsr <- c(4.8, 3.31, 3.74, 0.69, 2, 0.19, 0.41, 0.34, 0.3, 9) # M1=2, M2=1
params14gsr <- c(1.3, 2.2, 1.4, 1.5, 0.8, 2.4, 4.6, 0.4, 1, 0.35, 0.3, 0.2, 12, 13) # M1=2, M2=2

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

test_that("quantile_residuals gives correct residuals for non-restricted models", {
  expect_equal(quantile_residuals_int(simudata, 1, 1, params11, model="GMAR")[4:6], c(-0.2490894, -0.5921886, -0.2115230), tolerance=1e-3)
  expect_equal(quantile_residuals_int(simudata, 1, 1, params11t, model="StMAR")[7:9], c(0.1197390, -1.5075963, 0.8017146), tolerance=1e-4)
  expect_equal(quantile_residuals_int(simudata, 1, 2, params12, model="GMAR")[1:3], c(0.5364783, 0.9984865, -0.2941690), tolerance=1e-4)
  expect_equal(quantile_residuals_int(simudata, 1, 2, params12t, model="StMAR")[13:14], c(-1.6866567, 0.8906252), tolerance=1e-4)
  expect_equal(quantile_residuals_int(simudata[1:10], 2, 2, params22)[1:2], c(0.7442086, 0.1990116), tolerance=1e-4)
  expect_equal(quantile_residuals_int(3*simudata, 2, 2, params22)[113:115], c(4.395782, 3.864548, 3.409263), tolerance=1e-4)
  expect_equal(quantile_residuals_int(simudata, 2, 2, params22t, model="StMAR")[113:115], c(0.13373040, 0.05373437, -0.35308351), tolerance=1e-4)
  expect_equal(quantile_residuals_int(simudata, 2, 3, params23)[1:2], c(0.6549603, -0.0296626), tolerance=1e-4)
  expect_equal(quantile_residuals_int(3*simudata, 2, 3, params23)[1:3], c(5.950277, 4.238304, 3.559215), tolerance=1e-4)
  expect_equal(quantile_residuals_int(simudata, 1, c(1, 1), params12gs, model="G-StMAR")[10:12], c(1.106072, 1.208008, -1.113996), tolerance=1e-4)
  expect_equal(quantile_residuals_int(simudata, 2, c(1, 2), params23gs, model="G-StMAR")[100:102], c(1.920744, 1.949446, 1.834405), tolerance=1e-4)
  expect_equal(quantile_residuals_int(simudata, 2, c(1, 1), params22gs, model="G-StMAR")[20:22], c(2.0461076, 0.7647510, -0.6827557), tolerance=1e-4)
})

test_that("quantile_residuals gives correct residuals for restricted models", {
  expect_equal(quantile_residuals_int(simudata, 1, 2, params12r, restricted=TRUE)[1:3], c(0.2269731, 1.3225061, -0.4425795), tolerance=1e-4)
  expect_equal(quantile_residuals_int(simudata, 1, 2, params12r, restricted=TRUE)[12:14], c(-1.55925652, -1.77122955, 0.06400736), tolerance=1e-4)
  expect_equal(quantile_residuals_int(simudata, 1, 2, params12tr, model="StMAR", restricted=TRUE)[2:3], c(0.7061444, 0.1243444), tolerance=1e-4)
  expect_equal(quantile_residuals_int(simudata, 2, 3, params23tr, model="StMAR", restricted=TRUE)[16:18], c(1.175297890, 0.007549432, 0.533408786), tolerance=1e-4)
  expect_equal(quantile_residuals_int(simudata, 1, c(2, 1), params13gsr, model="G-StMAR", restricted=TRUE)[7:8], c(0.008349508, -1.692564492), tolerance=1e-4)
  expect_equal(quantile_residuals_int(simudata, 1, c(2, 2), params14gsr, model="G-StMAR", restricted=TRUE)[70:71], c(0.2219573, 0.5713081), tolerance=1e-4)
})

test_that("quantile_residuals gives correct residuals for constrained models", {
  expect_equal(quantile_residuals_int(simudata, 2, 1, params21c, model="StMAR", constraints=list(R3))[140:141], c(0.8115135, -0.2643541), tolerance=1e-4)
  expect_equal(quantile_residuals_int(simudata, 2, 2, params22c, model="StMAR", constraints=list(R4, R3))[1:2], c(2.173890, 2.017603), tolerance=1e-4)
  expect_equal(quantile_residuals_int(3*simudata, 2, 2, params22c_2, constraints=list(R4, R4))[113:115], c(4.395782, 3.864548, 3.409263), tolerance=1e-4)
  expect_equal(quantile_residuals_int(simudata, 2, 2, params22cr, model="StMAR", restricted=TRUE, constraints=R3)[24:25], c(0.1822218, 0.3833481), tolerance=1e-4)
  expect_equal(quantile_residuals_int(simudata, 3, 2, params32cr, restricted=TRUE, constraints=R1)[2:3], c(7.633364, 7.321785), tolerance=1e-1)
  expect_equal(quantile_residuals_int(simudata, 2, 3, params23cr, restricted=TRUE, constraints=R4)[130:131], c(0.49114035, 0.07224426), tolerance=1e-4)
  expect_equal(quantile_residuals_int(simudata, 3, c(1, 1), params32gsc, model="G-StMAR", restricted=FALSE, constraints=list(R1, R2))[13:14], c(1.016407, 1.215752), tolerance=1e-4)
  expect_equal(quantile_residuals_int(simudata, 2, c(1, 1), params22gsrc, model="G-StMAR", restricted=TRUE, constraints=R3)[1:2], c(1.633429, 1.416753), tolerance=1e-4)
})

library(uGMAR)
context("unconditional moments")

params11t <- c(-2, 0.8, 1, 12) # StMAR
params23 <- c(2.7, 0.8, -0.06, 0.3, 3.5, 0.8, -0.07, 2.6, 7.2, 0.3, -0.01, 0.1, 0.6, 0.25)
params23t <- c(1.8, 0.9, -0.06, 0.4, 7.2, 0.3, -0.009, 0.1, 3.1, 0.83, -0.05, 3.7, 0.7, 0.2, 11, 339, 198) # StMAR
params13gsr <- c(1.3, 1, 1.4, 0.8, 0.4, 2, 0.2, 0.25, 0.15, 20) # M1=2, M2=1

R1 <- matrix(c(1, 0, 0, 0, 0, 1), ncol=2)
R2 <- diag(1, ncol=3, nrow=3)
R3 <- matrix(c(0.5, 0.5), ncol=1)
params32c <- c(1, 0.1, -0.1, 1, 2, 0.2, -0.2, 2, 0.6, 11, 12) # R1, R1, StMAR
params33c <- c(1, 0.1, 0.1, 0.1, 1, 2, 0.2, 0.2, 0.2, 2, 3, 0.3, -0.3, 3, 0.5, 0.4) # R2, R2, R1
params22gsrc <- c(1, 2, 0.5, 1, 2, 0.5, 10) # M1=1, M2=1, R3


stmar11 <- GSMAR(p=1, M=1, params=params11t, model="StMAR")
gmar23 <- GSMAR(p=2, M=3, params=params23, model="GMAR")
stmar23 <- GSMAR(p=2, M=3, params=params23t, model="StMAR")
gstmar13r <- GSMAR(p=1, M=c(2, 1), params=params13gsr, model="G-StMAR", restricted=TRUE)

stmar32c <- GSMAR(p=3, M=2, params=params32c, model="StMAR", constraints=list(R1, R1))
gmar33c <- GSMAR(p=3, M=3, params=params33c, model="GMAR", constraints=list(R2, R2, R1))
gstmar22cr <- GSMAR(p=2, M=c(1, 1), params=params22gsrc, model="G-StMAR", restricted=TRUE, constraints=R3)

test_that("get_regime_means gives correct values", {
  expect_equal(get_regime_means(stmar11), -10, tolerance=1e-5)
  expect_equal(get_regime_means(gmar23), c(10.38462, 12.96296, 10.14085), tolerance=1e-5)
  expect_equal(get_regime_means(stmar23), c(11.25000, 10.15515, 14.09091), tolerance=1e-5)
  expect_equal(get_regime_means(gstmar13r), c(6.5, 5.0, 7.0), tolerance=1e-5)
  expect_equal(get_regime_means(stmar32c), c(1, 2), tolerance=1e-5)
  expect_equal(get_regime_means(gmar33c), c(1.428571, 5.000000, 3.000000), tolerance=1e-5)
  expect_equal(get_regime_means(gstmar22cr), c(2, 4), tolerance=1e-5)
})

test_that("get_regime_autocovs gives correct values", {
  expect_equal(get_regime_autocovs(stmar11)[1,], 2.222222, tolerance=1e-5)
  expect_equal(get_regime_autocovs(gmar23)[2,], c(0.380391743, 3.129034388, 0.008677247), tolerance=1e-5)
  expect_equal(get_regime_autocovs(stmar23)[1,], c(1.22123317, 0.03261833, 7.81584037), tolerance=1e-5)
  expect_equal(get_regime_autocovs(gstmar13r)[1,], c(0.8888889, 4.444444, 0.4444444), tolerance=1e-5)
  expect_equal(get_regime_autocovs(stmar32c)[3,], c(-0.1020408, -0.4347826), tolerance=1e-5)
  expect_equal(get_regime_autocovs(gmar33c)[3,], c(0.1298701, 0.8333333, -1.0975610), tolerance=1e-5)
  expect_equal(get_regime_autocovs(gstmar22cr)[2,], c(0.4, 0.8), tolerance=1e-5)
})

test_that("get_regime_vars gives correct values", {
  expect_equal(get_regime_vars(stmar11), 2.777778, tolerance=1e-5)
  expect_equal(get_regime_vars(gmar23), c(0.6995407, 5.9247333, 0.1096873), tolerance=1e-5)
  expect_equal(get_regime_vars(stmar23), c(1.4383413, 0.1097063, 9.8875089), tolerance=1e-5)
  expect_equal(get_regime_vars(gstmar13r), c(1.1111111, 5.5555556, 0.5555556), tolerance=1e-5)
  expect_equal(get_regime_vars(stmar32c), c(1.020408, 2.173913), tolerance=1e-5)
  expect_equal(get_regime_vars(gmar33c), c(1.038961, 2.500000, 3.658537), tolerance=1e-5)
  expect_equal(get_regime_vars(gstmar22cr), c(1.2, 2.4), tolerance=1e-5)
})
